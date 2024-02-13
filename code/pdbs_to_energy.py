import os
import itertools
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def parse_pdb(file_path):
    c3_atoms = []

    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") and line[12:16].strip() == "C3'":
                atom_info = {
                    'residue_name': line[17:20].strip(),
                    'chain_id': line[21],
                    'residue_seq': int(line[22:26].strip()),
                    'x_coord': float(line[30:38].strip()),
                    'y_coord': float(line[38:46].strip()),
                    'z_coord': float(line[46:54].strip()),
                }
                c3_atoms.append(atom_info)

    return c3_atoms

def compute_intra_chain_distances(c3_atoms):
    distances = []

    for atom1, atom2 in itertools.combinations(c3_atoms, 2):
        if atom1['chain_id'] == atom2['chain_id'] and abs(atom1['residue_seq'] - atom2['residue_seq']) >= 3:
            distance = np.linalg.norm(np.array([atom1['x_coord'], atom1['y_coord'], atom1['z_coord']]) -
                                      np.array([atom2['x_coord'], atom2['y_coord'], atom2['z_coord']]))
            distances.append(distance)

    return distances

def compute_frequencies(distances):
    bins = np.linspace(0, 20, 21)  # 20 intervals from 0 to 20 Å
    counts, _ = np.histogram(distances, bins=bins)
    frequencies = counts / len(distances)
    return frequencies

def compute_log_ratios(observed_frequencies, reference_frequencies):
    log_ratios = -np.log(np.divide(observed_frequencies, reference_frequencies))
    log_ratios = np.where(np.isinf(log_ratios), 10, log_ratios)
    return log_ratios

def calculate_pseudo_energy(distances_dict: dict):
    reference_frequencies = compute_frequencies([value for values in distances_dict.values() for value in values])
    energy_dict = {}

    for key, value in distances_dict.items():
        observed_frequencies = compute_frequencies(value)
        log_ratios = compute_log_ratios(observed_frequencies, reference_frequencies)
        energy_dict[key] = log_ratios

    return energy_dict

def write_energy_dict(energy_dict: dict, file_type='txt', output_dir='pseudo_energies'):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if file_type == 'txt':
        for key, value in energy_dict.items():
            output_path = os.path.join(output_dir, f'{key}_scores.txt')
            np.savetxt(output_path, value, delimiter='\n', fmt='%s')

    elif file_type == 'csv':
        energy_df = pd.DataFrame.from_dict(energy_dict)
        energy_df.to_csv(os.path.join(output_dir, 'pseudo_energies.csv'))

    else:
        print(f'File type {file_type} is not supported')

def plot_distances_boxplots(distances_for_base_pair: dict, output_dir: str):
  # Convert the dictionary to a DataFrame
  df = pd.DataFrame.from_dict(distances_for_base_pair, orient='index').transpose()

  # Plot boxplots using Seaborn
  plt.figure(figsize=(10, 6))
  sns.set(style="whitegrid")
  sns.boxplot(data=df)
  plt.title('Boxplots of Distances Distributions')
  plt.savefig(f'{output_dir}/distances_boxplots.png')
  plt.close()

def plot_distances_histograms(distances_for_base_pair: dict, output_dir: str):
  fig, axs = plt.subplots(2, 5, figsize=(15, 6))
  for idx, (key, value) in enumerate(distances_for_base_pair.items()):
      row = idx // 5
      col = idx % 5
      sns.histplot(value, ax=axs[row, col])
      axs[row, col].set_title(key)

  plt.tight_layout()
  plt.savefig(os.path.join(output_dir, 'distances_histograms.png'))
  plt.close()

def main(pdb_dir, output_dir):
    base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']
    distances_for_base_pair = {base_pair: [] for base_pair in base_pairs}  # Initialize distances dict

    for pdb_file in os.listdir(pdb_dir):
        if pdb_file.endswith(".pdb"):
            pdb_file_path = os.path.join(pdb_dir, pdb_file)
            c3_atoms = parse_pdb(pdb_file_path)

            for base_pair in base_pairs:
                chain1_atoms = [atom for atom in c3_atoms if atom['residue_name'] == base_pair[0] and atom['chain_id'] == 'A']
                chain2_atoms = [atom for atom in c3_atoms if atom['residue_name'] == base_pair[1] and atom['chain_id'] == 'B']

                distances = compute_intra_chain_distances(chain1_atoms + chain2_atoms)
                distances_for_base_pair[base_pair].extend(distances)

    for base_pair in base_pairs:
        mean_distance = np.mean(distances_for_base_pair[base_pair])
        std_dev = np.std(distances_for_base_pair[base_pair])

        print(f'{base_pair} - Mean distance: {mean_distance:.2f}, Standard Deviation: {std_dev:.2f}')

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Save distances_for_base_pair as a JSON file
    with open(os.path.join(output_dir, 'distances_for_base_pair.json'), 'w') as json_file:
        json.dump(distances_for_base_pair, json_file)
    print(f"Distances saved to {output_dir}")

    plot_distances_boxplots(distances_for_base_pair, output_dir)
    plot_distances_histograms(distances_for_base_pair, output_dir)

    # Calculate pseudo energies and save them
    energy_dict = calculate_pseudo_energy(distances_for_base_pair)
    write_energy_dict(energy_dict, output_dir=output_dir)
    write_energy_dict(energy_dict, 'csv', output_dir=output_dir)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute distances and pseudo-energies for RNA structures.")
    parser.add_argument("pdb_dir", type=str, help="Path to the directory containing PDB files.")
    parser.add_argument("output_dir", type=str, help="Path to the output directory to save results.")
    args = parser.parse_args()

    main(args.pdb_dir, args.output_dir)