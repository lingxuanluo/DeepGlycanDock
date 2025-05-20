import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig
from tqdm import tqdm
import os

def read_fasta(file_path):
    """Read FASTA file and return dictionary of labels and sequences"""
    sequences = {}
    current_label = None
    current_sequence = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_label and current_sequence:
                    sequences[current_label] = ''.join(current_sequence)
                current_label = line[1:]  # Remove '>'
                current_sequence = []
            else:
                current_sequence.append(line)
    
    # Add the last sequence
    if current_label and current_sequence:
        sequences[current_label] = ''.join(current_sequence)
    
    return sequences

def process_sequence(sequence, client):
    """Process a single sequence and return transformed embeddings"""
    # Create ESMProtein object
    protein = ESMProtein(sequence=sequence)
    protein_tensor = client.encode(protein)
    
    # Get logits and embeddings
    logits_output = client.logits(
        protein_tensor, 
        LogitsConfig(sequence=True, return_embeddings=True)
    )
    
    # Transform embeddings
    # Original shape: [1, seq_len+2, 1152]
    embeddings = logits_output.embeddings.squeeze(0)  # [seq_len+2, 1152]
    # Remove begin and end tokens
    embeddings = embeddings[1:-1, :]  # [seq_len, 1152]
    # Pad to 1280 dimensions with zeros, move padding to GPU
    padding = torch.zeros(embeddings.shape[0], 1280 - 1152).to('cuda')  # [seq_len, 128]
    embeddings = torch.cat([embeddings, padding], dim=1)  # [seq_len, 1280]
    
    # Detach from GPU and move to CPU
    embeddings = embeddings.detach().cpu()
    
    return embeddings

def process_fasta_to_individual_pts(fasta_path, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load model
    client = ESMC.from_pretrained("esmc_600m").to("cuda")
    
    # Read FASTA file
    sequences = read_fasta(fasta_path)
    
    # Process sequences and save individual files
    for label, sequence in tqdm(sequences.items(), desc="Processing sequences"):
        # Define output path for the .pt file
        output_path = os.path.join(output_dir, f"{label}.pt")
        
        # Skip if the .pt file already exists
        if os.path.exists(output_path):
            print(f"Skipping {label}: {output_path} already exists")
            continue
        
        # Process the sequence to generate embeddings
        embeddings = process_sequence(sequence, client)
        
        # Create data structure for this sequence
        data = {
            'label': label,
            'representations': {33: embeddings}
        }
        
        # Save to individual .pt file
        torch.save(data, output_path)


def main():
    # Set up argument parser
    import argparse
    parser = argparse.ArgumentParser(description="Process FASTA file to individual .pt files")
    parser.add_argument('--fasta_path', type=str, required=True, 
                        help="Path to the input FASTA file")
    parser.add_argument('--output_dir', type=str, required=True, 
                        help="Directory to save output .pt files")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Call the processing function
    process_fasta_to_individual_pts(args.fasta_path, args.output_dir)
    print(f"Completed processing {args.fasta_path} to {args.output_dir}")

if __name__ == "__main__":
    main()