## Setup Instructions

### Creating the 'HMM Analysis' Environment

Install the required tools in a new Conda environment named "HMM" using the following command:

```bash
conda create -y -n HMM -c bioconda -c conda-forge python=3.11 seqkit hmmer muscle=3.8.1551 weblogo transdecoder easel diamond trimal pymsaviz alive-progress
```

To activate the newly created environment, execute:

```bash
conda activate HMM
```

### Installing Conda
For Conda installation, visit the following [link](https://docs.anaconda.com/free/miniconda/).

### Setting Up 'HMMer'

Select a desired working directory and clone the HMMerMe repository:

```bash
git clone git@github.com:plantgenomicslab/HMMerMe.git
```

If an SSH key is not yet set up, please generate one. Navigate to the 'HMMerMe' folder:

```bash
cd HMMerMe
```

## Executing 'HMMerMe'

Usage:

```bash
python run.py --input [INPUT DIRECTORY] --db [DATABASE DIRECTORY] --output [OPTIONAL: NAME YOUR OUTPUT DIRECTORY] --E [OPTIONAL: E-value FOR SEQUENCE. DEFAULT SET TO 1e-5] --domE [OPTIONAL: E-value FOR DOMAIN. DEFAULT SET TO 1e-10] --CPU [OPTIONAL: NUMBER OF CPU CORES TO UTILIZE. DEFAULT SET TO 2] --logging [OPTIONAL: LOG ALL COMMANDS, INCLUDING SUCCESS AND FAILURE] --visualization [OPTIONAL: CALL WEBLOGO AND PYMSAVIZ]
```

Sample usage:

```bash
python run.py --input Input/ --db Database/ --output two_species --E 1e-10 --logging --visualization
```

The main script is `run.py`, which requires two arguments: `--input` for the input directory which includes {species_name}.fasta. **It is crucial that the filenames for these FASTA files do not contain spaces or special characters like "_".** .

`--db` for the database directory which includes {domain_name}.hmm. Optional arguments include `--CPU` to specify the number of cores (default is 2), `--visualization` to generate data visualizations.
**It's essential to ensure that the filenames for both domain and FASTA files are free from spaces or special characters such as "_". Additionally, the names of the domain and FASTA files must match exactly.**

`--output` to generate a custom output directory based on users choice, `--E` to specifcy the E-value for sequence search, and `--domE` to specify the E-value for domain searches.

Execute the script as follows:

```bash
python run.py --input Input/ --db Database/ [OPTIONS --CPU, --output, --visualization]
```

For data visualization, include the `--visualization` flag:

```bash
python run.py --input Input/ --db Database/ --visualization
```
## Additional features to enhance the analysis include:
- An `--output` option for specifying the directory where output files will be stored, defaulting to `output`.
- A `--CPU` option to set the number of cores used, with a default value of 2.
- An  `--visualization` option to enable or disable WebLogo generation, which is disabled by default for streamlined analysis.
- `--E` to specifcy the E-value for sequence search, and `--domE` to specify the E-value for domain searches.
## Input file
The required format for input sequences is the FASTA format. The system is designed to accommodate multiple FASTA files in same directory simultaneously.
For the database, the expected file format is HMM, and similarly, the system supports processing multiple HMM files in same directory concurrently.

## Output Files Overview

The analysis generates several types of output files, outlined as follows:
- `{Species}.hmm_results`: Contains the results of homology searches with significant E-value scores.
- `{Species}_{Homology_domain}.list`: Lists identified domain genes within a species.
- `{Species}_{Homology_domain}.fasta`: Sequences of domain genes.
- `{Species}_{Homology_domain}_domain.bed`: BED file indicating the location of domain genes.
- `{Species}_{Homology_domain}_domain.fasta`: Fasta sequences derived from BED locations.
- `{Species}_conflict.list`: Lists genes with conflicting domains.
- `{Species}_conflict.fasta`: Sequences of genes with conflicts domains.
- `{Species}_{Homology_domain}_muslced_domain.fa`: Aligned sequences via Muscle.
- `{Species}_{Homology_domain}_muscled_trimal_domain.fa`: Trimmed alignments to remove excessive gaps.
- `{Species}_{Homology_domain}_muscled_combined_domain.afa`: Combined alignment files.
- `{Species}_{Homology_domain}_muscled_combined_trimal_domain.afa`: Combined and trimmed alignment files.
- `{Species}_counts.tab`: Counts of domain genes identified.
- the `counts.tab`. This file will give you the total amount of Domain genes that were distinguished from your `.hmm_search` file. This tab delimited file, you could open in Excel.
  
| AaegyptiLVPWY | RYamideLuqin | Prothoracicotropichormone | SIfamide | CCHamide1 |
|---------------|--------------|---------------------------|----------|-----------|
| AaegyptiLVPWY | 1            | 3                         | 1        | 1         |

- Visualization files: `.pdf` and `.png` formats for visual representation.
- `.pdf` and `.png` files are made. If you called `--visualziation` it will call Weblogo v3 (Crooks et al 2004) and Pymsaviz (https://pypi.org/project/pyMSAviz/) to create the `.pdf` and `.png` files respectively. Here are some visualization:
[AaegyptiLVPWY_Adipokinetic_Corazonin_muscled_combined_trimal_domain.pdf](https://github.com/plantgenomicslab/HMMerMe/files/14514243/AaegyptiLVPWY_Adipokinetic_Corazonin_muscled_combined_trimal_domain.pdf)
![image](https://github.com/plantgenomicslab/HMMerMe/assets/907041/e8dc4ad5-a2a2-4a74-8bc5-b99301bc080f)
![AaegyptiLVPWY_Adipokinetic_Corazonin_muscled_combined_trimal_domain](https://github.com/plantgenomicslab/HMMerMe/assets/137996393/75e025cf-6219-40db-9465-78a11d80f7c4)


### Detailed Explanation
The process begins within the `Database` directory, hosting `.hmm` profile files crucial for identifying homology between your query sequences and extensive sequence databases. The `Input` directory should contain species-specific fasta or fas files. This setup enables `run.py` to discover homologous sequences between your fasta files and the predefined HMM profiles.

Consider the example of the species AaegyptiLVPWY, with its fasta file containing 28,392 genes. The database directory houses 49 HMM profiles, utilized by `run.py` to find potential homologous matches. The output files provide a comprehensive set of data ranging from homology results to detailed lists and sequences of identified domain genes, including handling conflicts where domain genes share sequences but have different names.

Furthermore, the analysis details the alignment of sequences via the Muscle program and the subsequent trimming of alignments to remove excess gaps, culminating in a set of combined alignment files. The `--visualization` option, if utilized, employs Weblogo and Pymsaviz tools to create graphical representations of the aligned sequences, enhancing the interpretability of the results.

## Custom Database Creation
For custom database creation, follow the steps below using `muscle` for alignment, `trimal` for trimming alignments based on conservation and gap thresholds, `esl-reformat` for format conversion, and `hmmbuild` to construct the HMM profile:
```bash
muscle -in {Homology_domain}.fasta -out {Homology_domain}.aln -clw
trimal -in {Homology_domain}.aln -out {Homology_domain}_trimmed.aln -gt 0.50 -cons 60
esl-reformat stockholm {Homology_domain}_trimmed.aln > {Homology_domain}.sto
hmmbuild {Homology_domain}.hmm {Homology_domain}.sto
cp {Homology_domain}.hmm ./Database
cp {Homology_domain}_trimmed.aln ./Database_fasta/{Homology_domain}.fasta
```
**It is essential to ensure that the filenames for both domain and FASTA files are free from spaces or special characters such as "_". Additionally, the names of the domain and FASTA files must match exactly.**
**{Homology_domain}.hmm need to be in database folder and {Homology_domain}.fasta need to be in `./Database_fasta` folder.**

## Transcripts input
At present, our system exclusively supports protein sequences. For those interested in analyzing transcript sequences, we recommend utilizing TransDecoder, available at https://github.com/TransDecoder/TransDecoder/wiki.
