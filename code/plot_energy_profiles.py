import pandas as pd
import matplotlib.pyplot as plt

def plot_energy_profiles(input_csv, output_dir):
    # Read the CSV file containing energy scores
    energies = pd.read_csv(input_csv, index_col=0)
    energies.index = [i for i in range(1, 21)]
    energies = energies.T

    # Get the indexes and values from the DataFrame
    indexes = energies.index.tolist()
    values = energies.values

    # Define the number of rows and columns for subplots
    num_rows = 2
    num_cols = 5

    # Create a new figure
    fig, ax = plt.subplots(num_rows, num_cols, figsize=(20, 8))

    # Flatten the axes to iterate over them
    ax = ax.flatten()

    # Iterate over each subplot
    for i, ax_i in enumerate(ax):
        # Plot the line for the corresponding pair of residues
        ax_i.plot(range(1, 21), values[i])

        # Set labels and title
        ax_i.set_ylabel("Pseudo-energy")
        ax_i.set_xlabel("Distance")
        ax_i.set_title(indexes[i])

    # Adjust layout and save the plot
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pseudo-energies.png')
    plt.close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot energy scores profiles from a CSV file.")
    parser.add_argument("input_csv", type=str, help="Path to the input CSV file containing energy scores.")
    parser.add_argument("output_dir", type=str, help="Path to the output directory to save the plot.")
    args = parser.parse_args()

    plot_energy_profiles(args.input_csv, args.output_dir)