# RNA_folding_energy
## [Full Project Code in Interactive Mode](https://colab.research.google.com/drive/1P4vSfwAUmnN3j4nZ8-AhFFI1x49XBaTW?usp=sharing)

# RNA Structure Analysis Pipeline

This Python script allows you to analyze RNA structure data from PDB files, compute distances between C3' atoms, calculate pseudo-energies, and generate visualizations of the data.

## Requirements

- Python 3.x
- NumPy
- Pandas
- Matplotlib
- Seaborn

## Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/your_username/rna-structure-analysis.git
    ```

2. Navigate to the cloned directory:

    ```bash
    cd rna-structure-analysis
    ```

3. Install the required dependencies:

    ```bash
    pip install -r requirements.txt
    ```

## Usage

1. Run the script to calculate distances and pseudo-energies:

    ```bash
    python pdbs_to_energy.py /path/to/pdb/files /path/to/output/directory
    ```

    Replace `/path/to/pdb/files` with the directory containing your PDB files and `/path/to/output/directory` with the directory where you want to save the results.

2. Optional arguments:
   
    - `--write_distances_json`: Specify whether to save distances data as a JSON file. Default is True.

## Plotting Pseudo-energy Profiles

The script `plot_energy_profiles.py` plots pseudo-energy profiles from a CSV file containing energy scores.

### Requirements

- Python 3.x
- Pandas
- Matplotlib

### Usage

1. Run the script to plot energy profiles:

    ```bash
    python plot_energy_profiles.py /path/to/input/csv /path/to/output/directory
    ```

    Replace `/path/to/input/csv` with the path to the input CSV file containing energy scores and `/path/to/output/directory` with the directory where you want to save the plot.

### Output

The script generates a plot named `pseudo-energies.png`, which displays pseudo-energy profiles for each pair of residues.

## Estimating Structure Energy

The script `estimate_struct_energy.py` estimates the energy of RNA structures.

### Requirements

- Python 3.x
- NumPy
- Pandas

### Usage

1. Run the script to estimate structure energy:

    ```bash
    python estimate_struct_energy.py /path/to/input/pdb /path/to/input/energy_csv /path/to/output/directory --plot-all --plot-each
    ```

    Replace `/path/to/input/pdb` with the path to the input PDB file, `/path/to/input/energy_csv` with the path to the input CSV file containing energy data, and `/path/to/output/directory` with the directory where you want to save the results.

2. Optional arguments:
   
    - `--plot-all`: Plot all residue pairs on one plot.
    - `--plot-each`: Plot each residue pair on separate subplots.

### Output

The script generates plots displaying the estimated pseudo-energy profiles for RNA structure pairs.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)
