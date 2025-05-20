import math

from e3nn import o3
import torch
from torch import nn
from torch.nn import functional as F
from torch_cluster import radius, radius_graph
from torch_scatter import scatter, scatter_mean
import numpy as np
from e3nn.nn import BatchNorm
from utils import so3, torus
from datasets.process_mols import lig_feature_dims, rec_residue_feature_dims

# to use a graph to make massage passing ,between surface and residue,and then just use surface to update!
# compare with version 2 , this version use more layers to update surface nodes

class AtomEncoder(torch.nn.Module):
    def __init__(self, emb_dim, feature_dims, sigma_embed_dim, lm_embedding_type= None):
        # first element of feature_dims tuple is a list with the lenght of each categorical feature and the second is the number of scalar features
        super(AtomEncoder, self).__init__()
        self.atom_embedding_list = torch.nn.ModuleList()
        self.num_categorical_features = len(feature_dims[0])
        self.num_scalar_features = feature_dims[1] + sigma_embed_dim
        self.lm_embedding_type = lm_embedding_type
        for i, dim in enumerate(feature_dims[0]):
            emb = torch.nn.Embedding(dim, emb_dim)
            torch.nn.init.xavier_uniform_(emb.weight.data)
            self.atom_embedding_list.append(emb)

        if self.num_scalar_features > 0:
            self.linear = torch.nn.Linear(self.num_scalar_features, emb_dim)
        if self.lm_embedding_type is not None:
            if self.lm_embedding_type == 'esm':
                self.lm_embedding_dim = 1280
            else: raise ValueError('LM Embedding type was not correctly determined. LM embedding type: ', self.lm_embedding_type)
            self.lm_embedding_layer = torch.nn.Linear(self.lm_embedding_dim + emb_dim, emb_dim)
    def forward(self, x):
        x_embedding = 0
        if self.lm_embedding_type is not None:
            assert x.shape[1] == self.num_categorical_features + self.num_scalar_features + self.lm_embedding_dim
        else:
            assert x.shape[1] == self.num_categorical_features + self.num_scalar_features
        for i in range(self.num_categorical_features):
            x_embedding += self.atom_embedding_list[i](x[:, i].long())

        if self.num_scalar_features > 0:
            x_embedding += self.linear(x[:, self.num_categorical_features:self.num_categorical_features + self.num_scalar_features])
        if self.lm_embedding_type is not None:
            x_embedding = self.lm_embedding_layer(torch.cat([x_embedding, x[:, -self.lm_embedding_dim:]], axis=1))
        return x_embedding


class TensorProductConvLayer(torch.nn.Module):
    def __init__(self, in_irreps, sh_irreps, out_irreps, n_edge_features, residual=True, batch_norm=True, dropout=0.0,
                 hidden_features=None):
        super(TensorProductConvLayer, self).__init__()
        self.in_irreps = in_irreps
        self.out_irreps = out_irreps
        self.sh_irreps = sh_irreps
        self.residual = residual
        if hidden_features is None:
            hidden_features = n_edge_features

        self.tp = tp = o3.FullyConnectedTensorProduct(in_irreps, sh_irreps, out_irreps, shared_weights=False)

        self.fc = nn.Sequential(
            nn.Linear(n_edge_features, hidden_features),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_features, tp.weight_numel)
        )
        self.batch_norm = BatchNorm(out_irreps) if batch_norm else None

    def forward(self, node_attr, edge_index, edge_attr, edge_sh, out_nodes=None, reduce='mean'):

        edge_src, edge_dst = edge_index
        tp = self.tp(node_attr[edge_dst], edge_sh, self.fc(edge_attr))

        out_nodes = out_nodes or node_attr.shape[0]
        out = scatter(tp, edge_src, dim=0, dim_size=out_nodes, reduce=reduce)

        if self.residual:
            padded = F.pad(node_attr, (0, out.shape[-1] - node_attr.shape[-1]))
            out = out + padded

        if self.batch_norm:

            out = self.batch_norm(out)
        return out

class TensorProductLigConvLayer(torch.nn.Module):
    def __init__(self, in_irreps, sh_irreps, out_irreps, n_edge_features, residual=True, batch_norm=True, dropout=0.0,
                 hidden_features=None):
        super(TensorProductLigConvLayer, self).__init__()
        self.in_irreps = in_irreps
        self.out_irreps = out_irreps
        self.sh_irreps = sh_irreps
        self.residual = residual
        if hidden_features is None:
            hidden_features = n_edge_features

        self.tp = tp = o3.FullyConnectedTensorProduct(in_irreps, sh_irreps, out_irreps, shared_weights=False)

        self.fc = nn.Sequential(
            nn.Linear(n_edge_features, hidden_features),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_features, tp.weight_numel)
        )
        self.batch_norm = BatchNorm(out_irreps) if batch_norm else None
        self.global_projection = nn.Sequential(
            nn.Linear(1024, hidden_features),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_features, o3.Irreps(self.out_irreps).dim)
        )

    def forward(self, node_attr, edge_index, edge_attr, edge_sh, global_graph_embedding, ptr, add_in = False, out_nodes=None, reduce='mean'):

        edge_src, edge_dst = edge_index
        tp = self.tp(node_attr[edge_dst], edge_sh, self.fc(edge_attr))

        out_nodes = out_nodes or node_attr.shape[0]
        out = scatter(tp, edge_src, dim=0, dim_size=out_nodes, reduce=reduce)

        if global_graph_embedding is not None and add_in == True:
            global_graph_embedding = self.global_projection(global_graph_embedding)
            num_nodes = out.shape[0]
            num_nodes_per_ligand = ptr[1:] - ptr[:-1]
            ligand_global_attr_expanded = torch.repeat_interleave(global_graph_embedding, num_nodes_per_ligand, dim=0)
            out = out + ligand_global_attr_expanded

        if self.residual:
            padded = F.pad(node_attr, (0, out.shape[-1] - node_attr.shape[-1]))
            out = out + padded

        if self.batch_norm:
            out = self.batch_norm(out)
        return out

class TensorProductScoreModel(torch.nn.Module):
    def __init__(self, t_to_sigma, device, timestep_emb_func, in_lig_edge_features=10, in_rec_edge_features = 5, sigma_embed_dim=32, sh_lmax=2,
                 ns=16, nv=4, num_conv_layers=2, lig_max_radius=5, rec_max_radius=30, cross_max_distance=250,
                 center_max_distance=30, distance_embed_dim=32, cross_distance_embed_dim=32, no_torsion=False,
                 scale_by_sigma=True, use_second_order_repr=False, batch_norm=True,
                 dynamic_max_cross=False, dropout=0.0, lm_embedding_type=None, confidence_mode=False,
                 confidence_dropout=0, confidence_no_batchnorm=False, num_confidence_outputs=1):
        super(TensorProductScoreModel, self).__init__()
        self.t_to_sigma = t_to_sigma
        self.in_lig_edge_features = in_lig_edge_features
        self.sigma_embed_dim = sigma_embed_dim
        self.lig_max_radius = lig_max_radius
        self.rec_max_radius = rec_max_radius
        self.cross_max_distance = cross_max_distance
        self.dynamic_max_cross = dynamic_max_cross
        self.center_max_distance = center_max_distance
        self.distance_embed_dim = distance_embed_dim
        self.cross_distance_embed_dim = cross_distance_embed_dim
        self.sh_irreps = o3.Irreps.spherical_harmonics(lmax=sh_lmax)
        self.ns, self.nv = ns, nv
        self.scale_by_sigma = scale_by_sigma
        self.device = device
        self.no_torsion = no_torsion
        self.timestep_emb_func = timestep_emb_func
        self.confidence_mode = confidence_mode
        self.num_conv_layers = num_conv_layers

        self.lig_node_embedding = AtomEncoder(emb_dim=ns, feature_dims=lig_feature_dims, sigma_embed_dim=sigma_embed_dim)
        self.lig_global_attr = nn.Sequential(nn.Linear(1024, ns//2), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns//2, ns//2))
        self.lig_node_attr_final = nn.Sequential(nn.Linear(int(ns*1.5), ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.lig_edge_embedding = nn.Sequential(nn.Linear(in_lig_edge_features + sigma_embed_dim + distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.rec_node_embedding = AtomEncoder(emb_dim=ns, feature_dims=rec_residue_feature_dims, sigma_embed_dim=0, lm_embedding_type=lm_embedding_type)
        self.rec_edge_embedding = nn.Sequential(nn.Linear(in_rec_edge_features + distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        # surface embeddings
        self.surface_node_embedding = AtomEncoder(emb_dim=ns, feature_dims=[[],4], sigma_embed_dim=sigma_embed_dim)
        self.surface_edge_embedding = nn.Sequential(nn.Linear(3 + sigma_embed_dim + distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.cross_edge_embedding = nn.Sequential(nn.Linear(sigma_embed_dim + cross_distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.surface_rec_cross_edge_embedding = nn.Sequential(nn.Linear(cross_distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.lig_distance_expansion = GaussianSmearing(0.0, lig_max_radius, distance_embed_dim)
        self.rec_distance_expansion = GaussianSmearing(0.0, rec_max_radius, distance_embed_dim)
        self.surface_distance_expansion = GaussianSmearing(0.0, rec_max_radius, distance_embed_dim)

        self.cross_distance_expansion = GaussianSmearing(0.0, cross_max_distance, cross_distance_embed_dim)

        if use_second_order_repr:
            irrep_seq = [
                f'{ns}x0e',
                f'{ns}x0e + {nv}x1o + {nv}x2e',
                f'{ns}x0e + {nv}x1o + {nv}x2e + {nv}x1e + {nv}x2o',
                f'{ns}x0e + {nv}x1o + {nv}x2e + {nv}x1e + {nv}x2o + {ns}x0o'
            ]
        else:
            irrep_seq = [
                f'{ns}x0e',
                f'{ns}x0e + {nv}x1o',
                f'{ns}x0e + {nv}x1o + {nv}x1e',
                f'{ns}x0e + {nv}x1o + {nv}x1e + {ns}x0o'
            ]
        lig_conv_layers= []
        lig_conv_layers1= []
        # surface modules
        surface_conv_layers,lig_to_surface_conv_layers, surface_to_lig_conv_layers = [], [],[]
        residue_to_surface_conv_layers = []
        rec_conv_layers = []
        for i in range(num_conv_layers):
            in_irreps = irrep_seq[min(i, len(irrep_seq) - 1)]
            out_irreps = irrep_seq[min(i + 1, len(irrep_seq) - 1)]
            parameters = {
                'in_irreps': in_irreps,
                'sh_irreps': self.sh_irreps,
                'out_irreps': out_irreps,
                'n_edge_features': 3 * ns,
                'hidden_features': 3 * ns,
                'residual': False,
                'batch_norm': batch_norm,
                'dropout': dropout
            }
            if i ==0:
                residue_to_surface_conv_layers.append(TensorProductConvLayer(** {
                'in_irreps': f'{ns}x0e + {nv}x1o + {nv}x1e + {ns}x0o',
                'sh_irreps': self.sh_irreps,
                'out_irreps':  in_irreps,
                'n_edge_features': 3 * ns,
                'hidden_features': 3 * ns,
                'residual': False,
                'batch_norm': batch_norm,
                'dropout': dropout
            }))
                rec_conv_layers.append(TensorProductConvLayer(** {
                'in_irreps': in_irreps,
                'sh_irreps': self.sh_irreps,
                'out_irreps': f'{ns}x0e + {nv}x1o + {nv}x1e + {ns}x0o',
                'n_edge_features': 3 * ns,
                'hidden_features': 3 * ns,
                'residual': False,
                'batch_norm': batch_norm,
                'dropout': dropout
            }))

            lig_layer = TensorProductConvLayer(**parameters)
            
            lig_conv_layers.append(lig_layer)
            if i == 0:
                lig_layer1 = TensorProductLigConvLayer(**parameters)
                lig_conv_layers1.append(lig_layer1)

            if i != num_conv_layers - 1:

                # surface layers
                surface_layer = TensorProductConvLayer(**parameters)
                surface_conv_layers.append(surface_layer)
                lig_to_surface_layer = TensorProductConvLayer(**parameters)
                lig_to_surface_conv_layers.append(lig_to_surface_layer)

   
            # surface layers
            surface_to_lig_layer = TensorProductConvLayer(**parameters)
            surface_to_lig_conv_layers.append(surface_to_lig_layer)

        self.lig_conv_layers = nn.ModuleList(lig_conv_layers)
        self.lig_conv_layers1 = nn.ModuleList(lig_conv_layers1)
        self.rec_conv_layers = nn.ModuleList(rec_conv_layers)

        # surface cross residue layer
        self.residue_to_surface_conv_layers = nn.ModuleList(residue_to_surface_conv_layers)
        # surface layers
        self.surface_conv_layers = nn.ModuleList(surface_conv_layers)
        self.lig_to_surface_conv_layers = nn.ModuleList(lig_to_surface_conv_layers)
        self.surface_to_lig_conv_layers = nn.ModuleList(surface_to_lig_conv_layers)

        if self.confidence_mode:
            self.confidence_predictor = nn.Sequential(
                nn.Linear(2*self.ns if num_conv_layers >= 3 else self.ns,ns),
                nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                nn.ReLU(),
                nn.Dropout(confidence_dropout),
                nn.Linear(ns, ns),
                nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                nn.ReLU(),
                nn.Dropout(confidence_dropout),
                nn.Linear(ns, num_confidence_outputs)
            )
        else:
            # center of mass translation and rotation components

            self.bind_or_not_predictor = nn.Sequential(
                nn.Linear(3*self.ns+12, ns),
                nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                nn.ReLU(),
                nn.Linear(ns, ns),
                nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                nn.ReLU(),
                nn.Linear(ns, 1)
            )

            self.residue_bind_or_not_predictor = nn.Sequential(
                nn.Linear(3*self.ns+12, ns),
                nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                nn.ReLU(),
                nn.Linear(ns, ns),
                nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                nn.ReLU(),
                nn.Linear(ns, 1)
            )

            self.center_distance_expansion = GaussianSmearing(0.0, center_max_distance, distance_embed_dim)
            self.center_edge_embedding = nn.Sequential(
                nn.Linear(distance_embed_dim + sigma_embed_dim, ns),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(ns, ns)
            )
            self.final_conv = TensorProductConvLayer(
                in_irreps=self.lig_conv_layers[-1].out_irreps,
                sh_irreps=self.sh_irreps,
                out_irreps=f'2x1o + 2x1e',
                n_edge_features=2 * ns,
                residual=False,
                dropout=dropout,
                batch_norm=batch_norm
            )
            self.tr_final_layer = nn.Sequential(nn.Linear(1 + sigma_embed_dim, ns),nn.Dropout(dropout), nn.ReLU(), nn.Linear(ns, 1))
            self.rot_final_layer = nn.Sequential(nn.Linear(1 + sigma_embed_dim, ns),nn.Dropout(dropout), nn.ReLU(), nn.Linear(ns, 1))

            if not no_torsion:
                # torsion angles components
                self.final_edge_embedding = nn.Sequential(
                    nn.Linear(distance_embed_dim, ns),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                    nn.Linear(ns, ns)
                )
                self.final_tp_tor = o3.FullTensorProduct(self.sh_irreps, "2e")
                self.tor_bond_conv = TensorProductConvLayer(
                    in_irreps=self.lig_conv_layers[-1].out_irreps,
                    sh_irreps=self.final_tp_tor.irreps_out,
                    out_irreps=f'{ns}x0o + {ns}x0e',
                    n_edge_features=3 * ns,
                    residual=False,
                    dropout=dropout,
                    batch_norm=batch_norm
                )
                self.tor_final_layer = nn.Sequential(
                    nn.Linear(2 * ns, ns, bias=False),
                    nn.Tanh(),
                    nn.Dropout(dropout),
                    nn.Linear(ns, 1, bias=False)
                )
    def forward(self, data):
        if not self.confidence_mode:
            tr_sigma, rot_sigma, tor_sigma = self.t_to_sigma(*[data.complex_t[noise_type] for noise_type in ['tr', 'rot', 'tor']])
        else:
            tr_sigma, rot_sigma, tor_sigma = [data.complex_t[noise_type] for noise_type in ['tr', 'rot', 'tor']]

        # build ligand graph
        lig_node_attr, lig_edge_index, lig_edge_attr, lig_edge_sh = self.build_lig_conv_graph(data)
        lig_src, lig_dst = lig_edge_index
        lig_node_attr = self.lig_node_embedding(lig_node_attr)
        lig_edge_attr = self.lig_edge_embedding(lig_edge_attr)
        # build receptor graph
        rec_node_attr, rec_edge_index, rec_edge_attr, rec_edge_sh = self.build_rec_conv_graph(data)
        rec_src, rec_dst = rec_edge_index
        rec_node_attr = self.rec_node_embedding(rec_node_attr)
        rec_edge_attr = self.rec_edge_embedding(rec_edge_attr)
        # build surface graph
        surface_node_attr, surface_edge_index, surface_edge_attr, surface_edge_sh = self.build_surface_conv_graph(data)
        surface_src, surface_dst = surface_edge_index
        surface_node_attr = self.surface_node_embedding(surface_node_attr)
        surface_edge_attr = self.surface_edge_embedding(surface_edge_attr)
        # use a layer to get residue feature to surface nodes
        # then  drop residus nodes
        # build cross graph
        if self.dynamic_max_cross:
            # this distance may can be changed for given pocket
            cross_cutoff = (tr_sigma * 3 + 10).unsqueeze(1)
        else:
            cross_cutoff = self.cross_max_distance

        # surface cross graph build
        surface_cross_edge_index, surface_cross_edge_attr, surface_cross_edge_sh = self.build_surface_cross_conv_graph(data, cross_cutoff)
        surface_cross_lig, surface_cross_rec = surface_cross_edge_index
        surface_cross_edge_attr = self.cross_edge_embedding(surface_cross_edge_attr)

        # surface ,residue cross graph builld this info will use one shot
        surface_rec_cross_edge_index, surface_rec_cross_edge_attr, surface_rec_cross_edge_sh = self.build_surface_rec_cross_conv_graph(data)
        surface_rec_cross_rec, surface_rec_cross_surface = surface_rec_cross_edge_index
        surface_rec_cross_edge_attr = self.surface_rec_cross_edge_embedding( surface_rec_cross_edge_attr)

        residue_to_surface_edge_attr_ = torch.cat([surface_rec_cross_edge_attr, rec_node_attr[ surface_rec_cross_rec, :self.ns], surface_node_attr[surface_rec_cross_surface, :self.ns]], -1)
        
        # update receptor embedding and then update embedding to surface
        rec_edge_attr_ = torch.cat([rec_edge_attr, rec_node_attr[rec_src, :self.ns], rec_node_attr[rec_dst, :self.ns]], -1)
        rec_intra_update = self.rec_conv_layers[0](rec_node_attr, rec_edge_index, rec_edge_attr_, rec_edge_sh)
        rec_node_attr = F.pad(rec_node_attr, (0, rec_intra_update.shape[-1] - rec_node_attr.shape[-1]))
        rec_node_attr = rec_node_attr + rec_intra_update

        # just one layer for feature update ,maybe can add more layers?
        surface_inter_residue_update = self.residue_to_surface_conv_layers[0](rec_node_attr, torch.flip(surface_rec_cross_edge_index,dims = [0]), residue_to_surface_edge_attr_, surface_rec_cross_edge_sh,
                                                            out_nodes=surface_node_attr.shape[0])
        
        surface_node_attr = F.pad(surface_node_attr, (0, surface_inter_residue_update.shape[-1] - surface_node_attr.shape[-1]))
        surface_node_attr = surface_node_attr + surface_inter_residue_update

        for l in range(len(self.lig_conv_layers)):

            # intra graph message passing
            lig_edge_attr_ = torch.cat([lig_edge_attr, lig_node_attr[lig_src, :self.ns], lig_node_attr[lig_dst, :self.ns]], -1)
            if l == 0:
                lig_intra_update = self.lig_conv_layers1[l](lig_node_attr, lig_edge_index, lig_edge_attr_, lig_edge_sh, data['ligand']['ligand_embedding'].reshape(-1, 1024), data['ligand']['ptr'], add_in = True)
            else:
                lig_intra_update = self.lig_conv_layers[l](lig_node_attr, lig_edge_index, lig_edge_attr_, lig_edge_sh)


            # surface inter graph message passing
            surface_to_lig_edge_attr_ = torch.cat([surface_cross_edge_attr, lig_node_attr[surface_cross_lig, :self.ns], surface_node_attr[surface_cross_rec, :self.ns]], -1)
            surface_lig_inter_update = self.surface_to_lig_conv_layers[l](surface_node_attr, surface_cross_edge_index, surface_to_lig_edge_attr_, surface_cross_edge_sh,
                                                              out_nodes=lig_node_attr.shape[0])
            
            if l != len(self.lig_conv_layers) - 1:

                # surface intra graph message passing
                surface_edge_attr_ = torch.cat([surface_edge_attr, surface_node_attr[surface_src, :self.ns], surface_node_attr[surface_dst, :self.ns]], -1)
                surface_intra_update = self.surface_conv_layers[l](surface_node_attr, surface_edge_index, surface_edge_attr_, surface_edge_sh)

                # lig to surface inter graph message passing
                lig_to_surface_edge_attr_ = torch.cat([surface_cross_edge_attr, lig_node_attr[surface_cross_lig, :self.ns], surface_node_attr[surface_cross_rec, :self.ns]], -1)
                
                surface_inter_update = self.lig_to_surface_conv_layers[l](lig_node_attr, torch.flip(surface_cross_edge_index, dims=[0]), lig_to_surface_edge_attr_, surface_cross_edge_sh,
                                                              out_nodes=surface_node_attr.shape[0])

            # padding original features
            lig_node_attr = F.pad(lig_node_attr, (0, lig_intra_update.shape[-1] - lig_node_attr.shape[-1]))
            # update features with residual updates
            lig_node_attr = lig_node_attr + lig_intra_update + surface_lig_inter_update
            if l != len(self.lig_conv_layers) - 1:

                # surface update
                surface_node_attr = F.pad(surface_node_attr, (0, surface_intra_update.shape[-1] - surface_node_attr.shape[-1]))
                surface_node_attr = surface_node_attr + surface_intra_update + surface_inter_update

        # compute confidence score
        if self.confidence_mode:
            scalar_lig_attr = torch.cat([lig_node_attr[:,:self.ns],lig_node_attr[:,-self.ns:] ], dim=1) if self.num_conv_layers >= 3 else lig_node_attr[:,:self.ns]
            confidence = self.confidence_predictor(scatter_mean(scalar_lig_attr, data['ligand'].batch, dim=0)).squeeze(dim=-1)
            return confidence

        # compute translational and rotational score vectors
        center_edge_index, center_edge_attr, center_edge_sh = self.build_center_conv_graph(data)
        center_edge_attr = self.center_edge_embedding(center_edge_attr)
        center_edge_attr = torch.cat([center_edge_attr, lig_node_attr[center_edge_index[1], :self.ns]], -1)
        # print(lig_node_attr, center_edge_index, center_edge_attr, center_edge_sh,data.num_graphs)
        global_pred = self.final_conv(lig_node_attr, center_edge_index, center_edge_attr, center_edge_sh, out_nodes=data.num_graphs)

        tr_pred = global_pred[:, :3] + global_pred[:, 6:9]
        rot_pred = global_pred[:, 3:6] + global_pred[:, 9:]
        data.graph_sigma_emb = self.timestep_emb_func(data.complex_t['tr'])

        # fix the magnitude of translational and rotational score vectors
        tr_norm = torch.linalg.vector_norm(tr_pred, dim=1).unsqueeze(1)
        tr_pred = tr_pred / tr_norm * self.tr_final_layer(torch.cat([tr_norm, data.graph_sigma_emb], dim=1))
        rot_norm = torch.linalg.vector_norm(rot_pred, dim=1).unsqueeze(1)
        rot_pred = rot_pred / rot_norm * self.rot_final_layer(torch.cat([rot_norm, data.graph_sigma_emb], dim=1))

        if self.scale_by_sigma:
            tr_pred = tr_pred / tr_sigma.unsqueeze(1)
            rot_pred = rot_pred * so3.score_norm(rot_sigma.cpu()).unsqueeze(1).to(data['ligand'].x.device)

        if self.no_torsion or data['ligand'].edge_mask.sum() == 0: return tr_pred, rot_pred, torch.empty(0, device=self.device)
        # torsional components
        tor_bonds, tor_edge_index, tor_edge_attr, tor_edge_sh = self.build_bond_conv_graph(data)
        tor_bond_vec = data['ligand'].pos[tor_bonds[1]] - data['ligand'].pos[tor_bonds[0]]
        tor_bond_attr = lig_node_attr[tor_bonds[0]] + lig_node_attr[tor_bonds[1]]

        tor_bonds_sh = o3.spherical_harmonics("2e", tor_bond_vec, normalize=True, normalization='component')
        tor_edge_sh = self.final_tp_tor(tor_edge_sh, tor_bonds_sh[tor_edge_index[0]])

        tor_edge_attr = torch.cat([tor_edge_attr, lig_node_attr[tor_edge_index[1], :self.ns],
                                   tor_bond_attr[tor_edge_index[0], :self.ns]], -1)
        tor_pred = self.tor_bond_conv(lig_node_attr, tor_edge_index, tor_edge_attr, tor_edge_sh,
                                  out_nodes=data['ligand'].edge_mask.sum(), reduce='mean')
        tor_pred = self.tor_final_layer(tor_pred).squeeze(1)
        edge_sigma = tor_sigma[data['ligand'].batch][data['ligand', 'ligand'].edge_index[0]][data['ligand'].edge_mask]

        if self.scale_by_sigma:
            tor_pred = tor_pred * torch.sqrt(torch.tensor(torus.score_norm(edge_sigma.cpu().numpy())).float()
                                             .to(data['ligand'].x.device))
        
        atom_close_labels, receptor_close_labels = self.get_ligand_labels(data)
        ligand_pred = self.bind_or_not_predictor(lig_node_attr)
        residue_pred =  self.residue_bind_or_not_predictor(rec_node_attr)
        return tr_pred, rot_pred, tor_pred, ligand_pred, atom_close_labels.unsqueeze(1), residue_pred, receptor_close_labels.unsqueeze(1)

    def get_ligand_labels(self, data):
        atom_close_labels = []
        receptor_close_labels = []
        if 'orig_pos' not in data['ligand']:
            device = data['receptor']['atoms_pos'].device
            return torch.tensor([], dtype=torch.int, device=device), torch.tensor([], dtype=torch.int, device=device)
        
        for i in range(len(data['ligand']['orig_pos'])):
            start = data['receptor']['ptr'][i].item()
            end = data['receptor']['ptr'][i+1].item()
            receptor_sub_pos = data['receptor']['atoms_pos'][start:end].float()  # Shape: (num_residues, num_atoms_per_residue, 3)
            
            # Number of residues and atoms per residue from the shape
            num_residues = receptor_sub_pos.shape[0]
            num_atoms_per_residue = receptor_sub_pos.shape[1]
            
            # Identify valid receptor atoms (no NaNs in any coordinate)
            valid_mask = ~torch.isnan(receptor_sub_pos[..., 0])  # Shape: (num_residues, num_atoms_per_residue)
            ligand_pos = torch.from_numpy(data['ligand']['orig_pos'][i]).float().to(receptor_sub_pos.device)  # Shape: (num_ligand_atoms, 3)
            
            if receptor_sub_pos.numel() > 0 and torch.any(valid_mask):
                # Reshape receptor positions for distance calculation: (num_residues * num_atoms_per_residue, 3)
                receptor_atoms = receptor_sub_pos.reshape(-1, 3)  # Shape: (num_residues * num_atoms_per_residue, 3)
                valid_mask_flat = valid_mask.flatten()  # Shape: (num_residues * num_atoms_per_residue,)
                valid_receptor_atoms = receptor_atoms[valid_mask_flat]  # Shape: (num_valid_atoms, 3)
                
                if len(valid_receptor_atoms) > 0:
                    # Compute pairwise distances between all ligand and receptor atoms
                    distances = torch.cdist(ligand_pos, valid_receptor_atoms)  # Shape: (num_ligand_atoms, num_valid_atoms)
                    
                    # For ligand atoms: Find minimum distance to any receptor atom
                    min_distances_ligand = torch.min(distances, dim=1).values  # Shape: (num_ligand_atoms,)
                    ligand_labels = (min_distances_ligand < 4.5).int().tolist()
                    atom_close_labels.extend(ligand_labels)
                    
                    # For receptor residues: Find minimum distance per receptor atom to any ligand atom
                    min_distances_receptor = torch.min(distances, dim=0).values  # Shape: (num_valid_atoms,)
                    
                    # Map back to all receptor atoms (including invalid ones)
                    receptor_atom_labels = torch.zeros(receptor_atoms.shape[0], device=receptor_sub_pos.device)
                    receptor_atom_labels[valid_mask_flat] = (min_distances_receptor < 4.5).float()
                    
                    # Reshape to per-residue level and take max per residue (1 if any atom is close)
                    receptor_residue_labels = receptor_atom_labels.reshape(num_residues, num_atoms_per_residue).max(dim=1).values  # Shape: (num_residues,)
                    receptor_close_labels.extend(receptor_residue_labels.int().tolist())
                else:
                    # No valid receptor atoms
                    atom_close_labels.extend([0] * ligand_pos.shape[0])
                    receptor_close_labels.extend([0] * num_residues)
            else:
                # No receptor atoms or all invalid
                atom_close_labels.extend([0] * ligand_pos.shape[0])
                receptor_close_labels.extend([0] * num_residues)
        
        atom_close_labels = torch.tensor(atom_close_labels, dtype=torch.int, device=receptor_sub_pos.device)
        receptor_close_labels = torch.tensor(receptor_close_labels, dtype=torch.int, device=receptor_sub_pos.device)
        return atom_close_labels, receptor_close_labels

    def build_lig_conv_graph(self, data):
        # builds the ligand graph edges and initial node and edge features
        data['ligand'].node_sigma_emb = self.timestep_emb_func(data['ligand'].node_t['tr'])

        # compute edges
        radius_edges = radius_graph(data['ligand'].pos, self.lig_max_radius, data['ligand'].batch)
        edge_index = torch.cat([data['ligand', 'ligand'].edge_index, radius_edges], 1).long()
        edge_attr = torch.cat([
            data['ligand', 'ligand'].edge_attr,
            torch.zeros(radius_edges.shape[-1], self.in_lig_edge_features, device=data['ligand'].x.device)
        ], 0)

        # compute initial features
        edge_sigma_emb = data['ligand'].node_sigma_emb[edge_index[0].long()]
        edge_attr = torch.cat([edge_attr, edge_sigma_emb], 1)
        node_attr = torch.cat([data['ligand'].x, data['ligand'].node_sigma_emb], 1)

        src, dst = edge_index
        edge_vec = data['ligand'].pos[dst.long()] - data['ligand'].pos[src.long()]
        edge_length_emb = self.lig_distance_expansion(edge_vec.norm(dim=-1))

        edge_attr = torch.cat([edge_attr, edge_length_emb], 1)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')

        return node_attr, edge_index, edge_attr, edge_sh


    def build_surface_conv_graph(self, data):
        # builds the receptor initial node and edge embeddings
        # tr = data['receptor'].node_t['tr']
        tr = data['receptor'].node_t['tr'][0]
        # data
        # data['surface'].node_t['tr'] = tr * torch.ones(data['surface'].num_nodes).to(tr.device)

        data['surface'].node_sigma_emb = self.timestep_emb_func(tr * torch.ones(data['surface'].num_nodes).to(tr.device)) # tr rot and tor noise is all the same
        # surface may have nan in features
        node_attr = torch.cat([torch.nan_to_num(data['surface'].x), data['surface'].node_sigma_emb], 1)

        # this assumes the edges were already created in preprocessing since protein's structure is fixed
        edge_index = data['surface','surface_edge','surface'].edge_index
        src, dst = edge_index
        edge_vec = data['surface'].pos[dst.long()] - data['surface'].pos[src.long()]

        edge_length_emb = self.surface_distance_expansion(edge_vec.norm(dim=-1))

        edge_sigma_emb = data['surface'].node_sigma_emb[edge_index[0].long()]

        edge_attr = torch.cat([data['surface','surface_edge','surface'].edge_attr,edge_sigma_emb, edge_length_emb], 1).float()
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')

        return node_attr, edge_index, edge_attr, edge_sh
    def build_rec_conv_graph(self, data):
        # builds the receptor initial node and edge embeddings
        # data['receptor'].node_sigma_emb = self.timestep_emb_func(data['receptor'].node_t['tr']) # tr rot and tor noise is all the same
        node_attr = data['receptor'].x

        # this assumes the edges were already created in preprocessing since protein's structure is fixed
        edge_index = data['receptor', 'receptor'].edge_index
        src, dst = edge_index
        edge_vec = data['receptor'].pos[dst.long()] - data['receptor'].pos[src.long()]

        edge_length_emb = self.rec_distance_expansion(edge_vec.norm(dim=-1))
        # edge_sigma_emb = data['receptor'].node_sigma_emb[edge_index[0].long()]

        edge_attr = torch.cat([data['receptor', 'rec_contact', 'receptor'].edge_attr, edge_length_emb], 1).float()
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')

        return node_attr, edge_index, edge_attr, edge_sh

    def build_surface_cross_conv_graph(self, data, cross_distance_cutoff):
        # builds the cross edges between ligand and receptor
        if torch.is_tensor(cross_distance_cutoff):
            # different cutoff for every graph (depends on the diffusion time)
            edge_index = radius(data['surface'].pos / cross_distance_cutoff[data['surface'].batch],
                                data['ligand'].pos / cross_distance_cutoff[data['ligand'].batch], 1,
                                data['surface'].batch, data['ligand'].batch, max_num_neighbors=30)
        else:
            edge_index = radius(data['surface'].pos, data['ligand'].pos, cross_distance_cutoff,
                            data['surface'].batch, data['ligand'].batch, max_num_neighbors=30)
        src, dst = edge_index
        edge_vec = data['surface'].pos[dst.long()] - data['ligand'].pos[src.long()]

        edge_length_emb = self.cross_distance_expansion(edge_vec.norm(dim=-1))
        edge_sigma_emb = data['ligand'].node_sigma_emb[src.long()]
        edge_attr = torch.cat([edge_sigma_emb, edge_length_emb], 1)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')

        return edge_index, edge_attr, edge_sh
    def build_surface_rec_cross_conv_graph(self, data, cross_distance_cutoff = 15):
        edge_index = radius(data['surface'].pos, data['receptor'].pos, cross_distance_cutoff,
                        data['surface'].batch, data['receptor'].batch, max_num_neighbors=30)
        src, dst = edge_index
        edge_vec = data['surface'].pos[dst.long()] - data['receptor'].pos[src.long()]

        edge_length_emb = self.cross_distance_expansion(edge_vec.norm(dim=-1))
        # edge_sigma_emb = data['receptor'].node_sigma_emb[src.long()]
        edge_attr = edge_length_emb#torch.cat([edge_sigma_emb, edge_length_emb], 1)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')

        return edge_index, edge_attr, edge_sh


    def build_center_conv_graph(self, data):
        # builds the filter and edges for the convolution generating translational and rotational scores
        edge_index = torch.cat([data['ligand'].batch.unsqueeze(0), torch.arange(len(data['ligand'].batch)).to(data['ligand'].x.device).unsqueeze(0)], dim=0)

        center_pos, count = torch.zeros((data.num_graphs, 3)).to(data['ligand'].x.device), torch.zeros((data.num_graphs, 3)).to(data['ligand'].x.device)
        center_pos.index_add_(0, index=data['ligand'].batch, source=data['ligand'].pos)
        center_pos = center_pos / torch.bincount(data['ligand'].batch).unsqueeze(1)

        edge_vec = data['ligand'].pos[edge_index[1]] - center_pos[edge_index[0]]
        edge_attr = self.center_distance_expansion(edge_vec.norm(dim=-1))
        edge_sigma_emb = data['ligand'].node_sigma_emb[edge_index[1].long()]
        edge_attr = torch.cat([edge_attr, edge_sigma_emb], 1)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')
        return edge_index, edge_attr, edge_sh

    def build_bond_conv_graph(self, data):
        # builds the graph for the convolution between the center of the rotatable bonds and the neighbouring nodes
        bonds = data['ligand', 'ligand'].edge_index[:, data['ligand'].edge_mask].long()
        bond_pos = (data['ligand'].pos[bonds[0]] + data['ligand'].pos[bonds[1]]) / 2
        bond_batch = data['ligand'].batch[bonds[0]]
        edge_index = radius(data['ligand'].pos, bond_pos, self.lig_max_radius, batch_x=data['ligand'].batch, batch_y=bond_batch)

        edge_vec = data['ligand'].pos[edge_index[1]] - bond_pos[edge_index[0]]
        edge_attr = self.lig_distance_expansion(edge_vec.norm(dim=-1))

        edge_attr = self.final_edge_embedding(edge_attr)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')

        return bonds, edge_index, edge_attr, edge_sh


class GaussianSmearing(torch.nn.Module):
    # used to embed the edge distances
    def __init__(self, start=0.0, stop=5.0, num_gaussians=50):
        super().__init__()
        offset = torch.linspace(start, stop, num_gaussians)
        self.coeff = -0.5 / (offset[1] - offset[0]).item() ** 2
        self.register_buffer('offset', offset)

    def forward(self, dist):
        dist = dist.view(-1, 1) - self.offset.view(1, -1)
        return torch.exp(self.coeff * torch.pow(dist, 2))
