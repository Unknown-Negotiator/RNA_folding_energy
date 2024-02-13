import argparse
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt

def parse_pdb(file_path):
    atoms = []

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
                atoms.append(atom_info)

    return atoms

def compute_intra_chain_distances(atoms):
    distances = []

    for atom1, atom2 in itertools.combinations(atoms, 2):
        if atom1['chain_id'] == atom2['chain_id'] and abs(atom1['residue_seq'] - atom2['residue_seq']) >= 3:
            distance = np.linalg.norm(np.array([atom1['x_coord'], atom1['y_coord'], atom1['z_coord']]) -
                                      np.array([atom2['x_coord'], atom2['y_coord'], atom2['z_coord']]))
            distances.append(distance)

    return distances

def pdb_to_distances(pdb_file_path: str):
    if pdb_file_path.endswith(".pdb"):
        atoms = parse_pdb(pdb_file_path)
        base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']
        distances_for_base_pair = {}

        for base_pair in base_pairs:
            chain1_atoms = [atom for atom in atoms if atom['residue_name'] == base_pair[0] and atom['chain_id'] == 'A']
            chain2_atoms = [atom for atom in atoms if atom['residue_name'] == base_pair[1] and atom['chain_id'] == 'B']

            distances = compute_intra_chain_distances(chain1_atoms + chain2_atoms)
            distances_for_base_pair[base_pair] = distances

            mean_distance = np.mean(distances_for_base_pair[base_pair])
            std_dev = np.std(distances_for_base_pair[base_pair])

            print(f'{base_pair} - Mean distance: {mean_distance:.2f}, Standard Deviation: {std_dev:.2f}')

        return distances_for_base_pair

def interpolate_energy(energy_dict, residue_pair, distance):
    energies = energy_dict[residue_pair]
    lower_index = int(distance) - 1
    upper_index = lower_index + 1

    # Check if the lower and upper index are within the valid range of energy values
    if lower_index >= 0 and upper_index < len(energies):
        lower_energy = energies[lower_index]
        upper_energy = energies[upper_index]

        # Check if the reference score is not NaN
        if not np.isnan(lower_energy) and not np.isnan(upper_energy):
            lower_distance = lower_index + 0.5
            interpolated_energy = lower_energy + (upper_energy - lower_energy) * (distance - lower_distance)
            return interpolated_energy
    return np.nan

def estimate_structure_energy(energy_dict, test_distances):
    total_energy = 0.0
    dist_energy_dict = {}
    for residue_pair, distances in test_distances.items():
        valid_distances, valid_energies = [], []
        for distance in distances:
            if distance <= 20:
                energy = interpolate_energy(energy_dict, residue_pair, distance)
                valid_distances.append(distance)
                valid_energies.append(energy)
                if not np.isnan(energy):
                    total_energy += energy

        # Sort both lists based on distance
        sorted_indices = sorted(range(len(valid_distances)), key=lambda k: valid_distances[k])
        sorted_distances = [valid_distances[i] for i in sorted_indices]
        sorted_energies = [valid_energies[i] for i in sorted_indices]

        dist_energy_dict[residue_pair] = {'distances':sorted_distances, 'energies':sorted_energies}

    return total_energy, dist_energy_dict

def csv_to_dict(csv_file_path, orient='list'):
    """
    Convert a CSV file to a dictionary.

    Args:
        csv_file_path (str): Path to the CSV file.
        orient (str): Orientation of the output dictionary.
            Valid values are 'list', 'dict', 'series', 'split', and 'records'.
            Default is 'list'.

    Returns:
        dict: Dictionary representation of the CSV file.
    """
    df = pd.read_csv(csv_file_path, index_col=0)
    dict_from_csv = df.to_dict(orient=orient)
    return dict_from_csv

def plot_estimated_energy_all(residue_dict:dict, output_dir='.'):
    """
    Plot line graphs for all pairs of residues separately on one plot.

    Args:
        residue_dict (dict): Dictionary of dictionaries where the primary key is a residue pair.
            Each sub-dictionary contains keys "distances" and "energies" with corresponding lists.
    """
    # Create a new figure for the plot
    plt.figure(figsize=(10, 6))

    # Iterate over each residue pair
    for pair, data in residue_dict.items():
        distances = data['distances']
        energies = data['energies']

        # Plot the line graph for the current pair
        plt.plot(distances, energies, label=pair)

    # Add labels and legend
    plt.xlabel('Distance')
    plt.ylabel('Energy')
    plt.title('Energy vs Distance for Residue Pairs')
    plt.legend()
    plt.grid(True)

    plt.savefig(f'{output_dir}/predicted_pseudo_energies_all.png')

def plot_estimated_energy_each(residue_dict:dict, output_dir='.'):
    """
    Plot multiple subplots in one picture, with each subplot containing the line graph for a separate residue pair.

    Args:
        residue_dict (dict): Dictionary of dictionaries where the primary key is a residue pair.
            Each sub-dictionary contains keys "distances" and "energies" with corresponding lists.
    """
    # Define the number of rows and columns for subplots
    num_rows = 2
    num_cols = 5

    # Create a new figure
    fig, ax = plt.subplots(num_rows, num_cols, figsize=(20, 8))

    # Flatten the axes to iterate over them
    ax = ax.flatten()

    # Iterate over each subplot
    for i, (pair, data) in enumerate(residue_dict.items()):
        distances = data['distances']
        energies = data['energies']

        # Plot the line for the corresponding pair of residues
        ax[i].plot(distances, energies)

        # Set labels and title
        ax[i].set_ylabel("Pseudo-energy")
        ax[i].set_xlabel("Distance")
        ax[i].set_title(pair)

    # Adjust layout and save the plot
    plt.tight_layout()
    plt.savefig(f'{output_dir}/predicted_pseudo-energies_each.png')

def main(args):
    # Call functions based on arguments
    test_distances = pdb_to_distances(args.pdb_file)
    energy_dict = csv_to_dict(args.energy_csv)

    estimated_energy, dist_energy_dict = estimate_structure_energy(energy_dict, test_distances)
    print("Estimated Gibbs free energy:", estimated_energy)

    if args.plot_all:
        plot_estimated_energy_all(dist_energy_dict, args.output_dir)
    if args.plot_each:
        plot_estimated_energy_each(dist_energy_dict, args.output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files and plot estimated energies.")
    parser.add_argument("pdb_file", help="Path to the PDB file")
    parser.add_argument("energy_csv", help="Path to the CSV file containing energy data")
    parser.add_argument("output_dir", default=".", help="Output directory for plots (default: current directory)")
    parser.add_argument("--plot-all", action="store_true", help="Plot all residue pairs on one plot")
    parser.add_argument("--plot-each", action="store_true", help="Plot each residue pair on separate subplots")
    args = parser.parse_args()

    main(args)