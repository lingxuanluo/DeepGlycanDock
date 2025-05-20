import copy

import numpy as np
from torch_geometric.loader import DataLoader
from tqdm import tqdm
from utils import so3, torus
from utils_glycan.sampling import randomize_position, sampling
import torch
from utils_glycan.diffusion_utils import get_t_schedule
from torch_geometric.data import Dataset,Data
from loguru import logger

class ListDataset(Dataset):
    def __init__(self, list):
        super().__init__()
        self.data_list = list
    def len(self) -> int:
        return len(self.data_list)
    def get(self, idx: int) -> Data:
        return self.data_list[idx]
import gc
def loss_function(tr_pred, rot_pred, tor_pred, data, t_to_sigma, device, ligand_pred, atom_close_labels, residue_pred, residue_close_labels, tr_weight=1, rot_weight=1,
                  tor_weight=1, atom_type_weight=5, residue_type_weight=5, apply_mean=True, no_torsion=False):
    tr_sigma, rot_sigma, tor_sigma = t_to_sigma(
        *[data.complex_t[noise_type]
          for noise_type in ['tr', 'rot', 'tor']])
    mean_dims = (0, 1) if apply_mean else 1
    bce_loss_fn = torch.nn.BCEWithLogitsLoss(reduction='none')
    BCE_loss = bce_loss_fn(ligand_pred, atom_close_labels.float())  # Ensure labels are float   
    atom_close_labels = atom_close_labels.type(torch.long)        
    at = torch.tensor([0.25, 0.75]).to(device).gather(0, atom_close_labels.data.view(-1))      
    pt = torch.exp(-BCE_loss)        
    F_loss = at*(1-pt)**2 * BCE_loss        
    atom_type_loss = F_loss.mean()

    BCE_loss_residue = bce_loss_fn(residue_pred, residue_close_labels.float())
    residue_close_labels = residue_close_labels.type(torch.long)
    at = torch.tensor([0.25, 0.75]).to(device).gather(0, residue_close_labels.data.view(-1))
    pt = torch.exp(-BCE_loss_residue)
    residue_type_loss = at*(1-pt)**2 * BCE_loss_residue
    residue_type_loss = residue_type_loss.mean()

    tr_score = data.tr_score
    tr_sigma = tr_sigma.unsqueeze(-1)

    tr_loss = ((tr_pred - tr_score) ** 2 * tr_sigma ** 2).mean(dim=mean_dims)
    tr_base_loss = (tr_score ** 2 * tr_sigma ** 2).mean(dim=mean_dims)
    # rotation component
    rot_score = data.rot_score
    rot_score_norm = so3.score_norm(rot_sigma.cpu()).unsqueeze(-1).to(device)
    rot_loss = (((rot_pred - rot_score) / rot_score_norm) ** 2).mean(dim=mean_dims)
    rot_base_loss = ((rot_score / rot_score_norm) ** 2).mean(dim=mean_dims)
    # torsion component
    if not no_torsion:

        edge_tor_sigma = torch.from_numpy(
            np.concatenate(data.tor_sigma_edge))
        
        tor_score = data.tor_score
        tor_score_norm2 = torch.tensor(torus.score_norm(edge_tor_sigma.cpu().numpy())).float().to(device)

        tor_loss = ((tor_pred - tor_score) ** 2 / tor_score_norm2)
        tor_base_loss = ((tor_score ** 2 / tor_score_norm2))
        if apply_mean:
            tor_loss, tor_base_loss = tor_loss.mean() * torch.ones(1, dtype=torch.float,device = device), tor_base_loss.mean() * torch.ones(1, dtype=torch.float,device = device)
        else:
            index = data['ligand'].batch[
                data['ligand', 'ligand'].edge_index[0][data['ligand'].edge_mask]]
            num_graphs = data.num_graphs
            t_l, t_b_l, c = torch.zeros(num_graphs,device = device), torch.zeros(num_graphs,device = device), torch.zeros(num_graphs,device = device)

            c.index_add_(0, index, torch.ones(tor_loss.shape,device = device))
            c = c + 0.0001
            t_l.index_add_(0, index, tor_loss)
            t_b_l.index_add_(0, index, tor_base_loss)
            tor_loss, tor_base_loss = t_l / c, t_b_l / c
    else:
        if apply_mean:
            tor_loss, tor_base_loss = torch.zeros(1, dtype=torch.float,device = device), torch.zeros(1, dtype=torch.float,device = device)
        else:
            tor_loss, tor_base_loss = torch.zeros(len(rot_loss), dtype=torch.float,device = device), torch.zeros(len(rot_loss), dtype=torch.float,device = device)

    loss = tr_loss * tr_weight + rot_loss * rot_weight + tor_loss * tor_weight + atom_type_loss * atom_type_weight + residue_type_loss * residue_type_weight
  
    return loss, tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss, atom_type_loss, residue_type_loss

class AverageMeter():
    def __init__(self, types, unpooled_metrics=False, intervals=1):
        self.types = types
        self.intervals = intervals
        self.count = 0 if intervals == 1 else torch.zeros(len(types), intervals)
        self.acc = {t: torch.zeros(intervals) for t in types}
        self.unpooled_metrics = unpooled_metrics

    def add(self, vals, interval_idx=None):
        if self.intervals == 1:
            self.count += 1 if vals[0].dim() == 0 else len(vals[0])
            for type_idx, v in enumerate(vals):
                self.acc[self.types[type_idx]] += v.sum() if self.unpooled_metrics else v
        else:
            for type_idx, v in enumerate(vals):

                self.count[type_idx].index_add_(0, interval_idx[type_idx], torch.ones(len(v)))
                if not torch.allclose(v, torch.tensor(0.0)):
                    self.acc[self.types[type_idx]].index_add_(0, interval_idx[type_idx], v)
    def summary(self):
        if self.intervals == 1:
            out = {k: v.item() / self.count for k, v in self.acc.items()}
            return out
        else:
            out = {}
            for i in range(self.intervals):
                for type_idx, k in enumerate(self.types):
                    out['int' + str(i) + '_' + k] = (
                            list(self.acc.values())[type_idx][i] / self.count[type_idx][i]).item()
            return out


def train_epoch(model, loader, optimizer, device, t_to_sigma, loss_fn, accelerator, ema_weights):
    model.train()
    # if mdn_mode:
        # 
    meter = AverageMeter(['loss', 'tr_loss', 'rot_loss', 'tor_loss', 'tr_base_loss', 'rot_base_loss', 'tor_base_loss', 'atom_type_loss', 'residue_type_loss'])
    pbar = tqdm(loader, total=len(loader),disable=not accelerator.is_local_main_process)
    for data in pbar:
        if device.type == 'cuda' and len(data) == 1 or device.type == 'cpu' and data.num_graphs == 1:
            logger.info("Skipping batch of size 1 since otherwise batchnorm would not work.")
        optimizer.zero_grad()
        try:
            tr_pred, rot_pred, tor_pred, ligand_pred, atom_close_labels, residue_pred, residue_close_labels  = model(data)
            with accelerator.autocast():
                loss, tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss, atom_type_loss, residue_type_loss = \
                    loss_fn(tr_pred, rot_pred, tor_pred, data, t_to_sigma, device, ligand_pred, atom_close_labels, residue_pred, residue_close_labels)
            if not torch.isnan(loss.mean()):
                accelerator.backward(loss)
                optimizer.step()
            else:
                logger.info(f'loss is nan in these data samples: {data.name}')
                # continue
                loss = torch.nan_to_num(loss)
            # gather all loss for plot
            loss, tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss, atom_type_loss, residue_type_loss = \
                accelerator.gather(loss),accelerator.gather(tr_loss), accelerator.gather(rot_loss), \
                    accelerator.gather(tor_loss), accelerator.gather(tr_base_loss), accelerator.gather(rot_base_loss), accelerator.gather(tor_base_loss), accelerator.gather(atom_type_loss), accelerator.gather(residue_type_loss)

            ema_weights.update(model.parameters())
            meter.add([loss.mean().cpu().detach(), tr_loss.mean().cpu().detach(), rot_loss.mean().cpu().detach(), tor_loss.mean().cpu().detach(), tr_base_loss.mean().cpu().detach(), 
                       rot_base_loss.mean().cpu().detach(), tor_base_loss.mean().cpu().detach(), atom_type_loss.mean().cpu().detach(), residue_type_loss.mean().cpu().detach()])
        except RuntimeError as e:
            if 'out of memory' in str(e):
                logger.info('| WARNING: ran out of memory, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                optimizer.zero_grad()
                del data

                gc.collect()
                torch.cuda.empty_cache()
                continue
            elif 'Input mismatch' in str(e):
                logger.info('| WARNING: weird torch_cluster error, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                optimizer.zero_grad()
                del data
                gc.collect()
                torch.cuda.empty_cache()
                continue
            else:
                raise e
    logger.info('clear last train batch data and model grad')
    for p in model.parameters():
        if p.grad is not None:
            del p.grad  # free some memory
    del data
    gc.collect()
    torch.cuda.empty_cache()
    return meter.summary()


def test_epoch(model, loader, device, t_to_sigma, loss_fn, accelerator, test_sigma_intervals=False,model_type = 'energy_score_model'):
    if not model_type == 'energy_score_model':
        model.eval()
    meter = AverageMeter(['loss', 'tr_loss', 'rot_loss', 'tor_loss', 'tr_base_loss', 'rot_base_loss', 'tor_base_loss', 'atom_type_loss', 'residue_type_loss'],
                         unpooled_metrics=True)

    if test_sigma_intervals:
        meter_all = AverageMeter(
            ['loss', 'tr_loss', 'rot_loss', 'tor_loss', 'tr_base_loss', 'rot_base_loss', 'tor_base_loss', 'atom_type_loss', 'residue_type_loss'],
            unpooled_metrics=True, intervals=10)
    
    for data in tqdm(loader, total=len(loader), disable=not accelerator.is_local_main_process):
        try:
            if not model_type == 'energy_score_model':
                with torch.no_grad():
                    tr_pred, rot_pred, tor_pred, ligand_pred, atom_close_labels, residue_pred, residue_close_labels = model(data)
            else:
                tr_pred, rot_pred, tor_pred, ligand_pred, atom_close_labels, residue_pred, residue_close_labels = model(data)
            with accelerator.autocast():
                loss, tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss, atom_type_loss, residue_type_loss = \
                    loss_fn(tr_pred, rot_pred, tor_pred, data, t_to_sigma, device, ligand_pred, atom_close_labels, residue_pred, residue_close_labels, apply_mean=False)


            loss, tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss, atom_type_loss, residue_type_loss = \
                accelerator.gather(loss),accelerator.gather(tr_loss), accelerator.gather(rot_loss), \
                    accelerator.gather(tor_loss), accelerator.gather(tr_base_loss), accelerator.gather(rot_base_loss), \
            accelerator.gather(tor_base_loss), accelerator.gather(atom_type_loss), accelerator.gather(residue_type_loss)

            metrics = [loss.mean().cpu().detach(), tr_loss.mean().cpu().detach(), \
                       rot_loss.mean().cpu().detach(), tor_loss.mean().cpu().detach(), \
                        tr_base_loss.mean().cpu().detach(), rot_base_loss.mean().cpu().detach(), tor_base_loss.mean().cpu().detach(), \
                        atom_type_loss.mean().cpu().detach(), residue_type_loss.mean().cpu().detach()]
            meter.add(metrics)


            if test_sigma_intervals > 0:
                complex_t_tr, complex_t_rot, complex_t_tor = [data.complex_t[noise_type] for
                                                              noise_type in ['tr', 'rot', 'tor']]
                sigma_index_tr = torch.round(complex_t_tr.cpu() * (10 - 1)).long().to(device)
                sigma_index_rot = torch.round(complex_t_rot.cpu() * (10 - 1)).long().to(device)
                sigma_index_tor = torch.round(complex_t_tor.cpu() * (10 - 1)).long().to(device)
                sigma_index_tr = accelerator.gather(sigma_index_tr).cpu().detach()
                sigma_index_rot = accelerator.gather(sigma_index_rot).cpu().detach()
                sigma_index_tor = accelerator.gather(sigma_index_tor).cpu().detach()
                meter_all.add(
                    metrics,
                    [sigma_index_tr, sigma_index_tr, sigma_index_rot, sigma_index_tor, sigma_index_tr, sigma_index_rot,
                     sigma_index_tor, sigma_index_tr])
        except RuntimeError as e:
            if 'out of memory' in str(e):
                logger.info('| WARNING: ran out of memory, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                del data
                gc.collect()
                torch.cuda.empty_cache()
                continue
            elif 'Input mismatch' in str(e):
                logger.info('| WARNING: weird torch_cluster error, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                del data
                gc.collect()
                torch.cuda.empty_cache()
                continue
            else:
                raise e
    logger.info('clear val batch data and model grad')
    for p in model.parameters():
        if p.grad is not None:
            del p.grad  # free some memory
    del data
    gc.collect()
    torch.cuda.empty_cache()
    out = meter.summary()
    if test_sigma_intervals > 0: out.update(meter_all.summary())
    return out


def inference_epoch(model, complex_graphs, device, t_to_sigma, args,accelerator):
    
    t_schedule = get_t_schedule(inference_steps=args.inference_steps)
    tr_schedule, rot_schedule, tor_schedule = t_schedule, t_schedule, t_schedule

    dataset = ListDataset(complex_graphs)
    loader = DataLoader(dataset=dataset, batch_size=1, shuffle=False)
    loader = accelerator.prepare(loader)
    rmsds = []
    logger.info(f'dataset size {len(dataset)}')
    for orig_complex_graph in tqdm(loader,disable=not accelerator.is_local_main_process):

        data_list = [copy.deepcopy(orig_complex_graph)]
        randomize_position(data_list, args.no_torsion, False, args.tr_sigma_max)

        predictions_list = None
        confidences = None
        failed_convergence_counter = 0
        while predictions_list == None:
            try:
                predictions_list, confidences, ligand_pred_list, residue_pred_list = sampling(input_data_list=data_list, model=model.module if device.type=='cuda' else model,
                                                         inference_steps=args.inference_steps,
                                                         tr_schedule=tr_schedule, rot_schedule=rot_schedule,
                                                         tor_schedule=tor_schedule,
                                                         device=device, t_to_sigma=t_to_sigma, model_args=args)
            except Exception as e:
                if 'failed to converge' in str(e):
                    failed_convergence_counter += 1
                    if failed_convergence_counter > 5:
                        logger.info('| WARNING: SVD failed to converge 5 times - skipping the complex')
                        break
                    logger.info('| WARNING: SVD failed to converge - trying again with a new sample')
                else:
                    raise e
        if failed_convergence_counter > 5: continue
        if args.no_torsion:
            orig_complex_graph['ligand'].orig_pos = (orig_complex_graph['ligand'].pos.cpu().numpy() +
                                                     orig_complex_graph.original_center.cpu().numpy())

        filterHs = torch.not_equal(predictions_list[0]['ligand'].x[:, 0], 0).cpu().numpy()

        if isinstance(orig_complex_graph['ligand'].orig_pos, list):
            orig_complex_graph['ligand'].orig_pos = orig_complex_graph['ligand'].orig_pos[0]

        ligand_pos = np.asarray(
            [complex_graph['ligand'].pos.cpu().numpy()[filterHs] for complex_graph in predictions_list])
        orig_ligand_pos = np.expand_dims(
            orig_complex_graph['ligand'].orig_pos[filterHs] - orig_complex_graph.original_center.cpu().numpy(), axis=0)
        rmsd = np.sqrt(((ligand_pos - orig_ligand_pos) ** 2).sum(axis=2).mean(axis=1))
        rmsds.append(rmsd)
    rmsds = np.array(rmsds)
    logger.info(f'rmsd: {rmsds}')
    losses = {'rmsds_lt2': (100 * (rmsds < 2).sum() / len(rmsds)),
              'rmsds_lt5': (100 * (rmsds < 5).sum() / len(rmsds))}
    del dataset, loader,predictions_list, confidences,ligand_pos, orig_ligand_pos, rmsd,filterHs
    gc.collect()
    torch.cuda.empty_cache()
    return losses

def inference_epoch_parallel(model, complex_graphs, device, t_to_sigma, args,accelerator):
    t_schedule = get_t_schedule(inference_steps=args.inference_steps)
    tr_schedule, rot_schedule, tor_schedule = t_schedule, t_schedule, t_schedule

    dataset = ListDataset(complex_graphs)
    loader = DataLoader(dataset=dataset, batch_size=args.batch_size, shuffle=False)
    loader = accelerator.prepare(loader)
    rmsds = []
    for orig_complex_graph in tqdm(loader,disable=not accelerator.is_local_main_process):
        orig_complex_graph_list = orig_complex_graph.to_data_list()
        data_list = [copy.deepcopy(graph) for graph in orig_complex_graph_list ]
        randomize_position(data_list, args.no_torsion, False, args.tr_sigma_max)

        predictions_list = None
        confidences = None
        failed_convergence_counter = 0
        while predictions_list == None:
            try:
                predictions_list, confidences, ligand_pred_list, residue_pred_list = sampling(input_data_list=data_list, model=model if device.type=='cuda' else model,
                                                         inference_steps=args.inference_steps,
                                                         tr_schedule=tr_schedule, rot_schedule=rot_schedule,
                                                         tor_schedule=tor_schedule,
                                                         device=device, t_to_sigma=t_to_sigma, model_args=args)
            except Exception as e:
                if 'failed to converge' in str(e):
                    failed_convergence_counter += 1
                    if failed_convergence_counter > 5:
                        logger.info('| WARNING: SVD failed to converge 5 times - skipping the complex')
                        break
                    logger.info('| WARNING: SVD failed to converge - trying again with a new sample')
                else:
                    raise e
        if failed_convergence_counter > 5: continue
        for pos_idx ,(predict_graph,orig_graph) in enumerate(zip(predictions_list,orig_complex_graph_list)):

            filterHs = torch.not_equal(predict_graph['ligand'].x[:, 0], 0)
            ligand_pos = predict_graph['ligand'].pos[filterHs].to(model.device)
            orig_ligand_pos = orig_graph['ligand'].pos[filterHs].to(model.device)
            # logger.info(f'{ligand_pos} ,{orig_ligand_pos}')
            rmsd = torch.sqrt(((ligand_pos - orig_ligand_pos) ** 2).sum())
            rmsds.append(rmsd)
    rmsds = torch.stack(rmsds)

    rmsds = accelerator.gather(rmsds)

    losses = {'rmsds_lt2': (100 * (rmsds < 2).sum() / len(rmsds)),
              'rmsds_lt5': (100 * (rmsds < 5).sum() / len(rmsds))}
    del dataset, loader,predictions_list, confidences,ligand_pos, orig_ligand_pos, rmsd,filterHs,complex_graphs
    gc.collect()
    torch.cuda.empty_cache()
    return losses