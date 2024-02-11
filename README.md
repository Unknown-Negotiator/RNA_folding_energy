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

1. Run the script with the following command:

    ```bash
    python pdbs_to_energy.py /path/to/pdb/files /path/to/output/directory
    ```

    Replace `/path/to/pdb/files` with the directory containing your PDB files and `/path/to/output/directory` with the directory where you want to save the results.

2. Optional arguments:
   
    - `--write_distances_json`: Specify whether to save distances data as a JSON file. Default is True.

## Output

The script generates the following output:

- `distances_for_base_pair.json`: JSON file containing computed distances for each base pair.
- `distances_boxplots.png`: Boxplot visualizations of distances distributions.
- `distances_histograms.png`: Histogram visualizations of distances distributions.
- Pseudo-energy scores saved as text and CSV files for each base pair.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)
