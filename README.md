# <i>Staphylococcus aureus</i> Lineage Typer (SaLTy)
SaLTy assigns a lineage to <i>Staphylococcus aureus</i> WGS data and is suitable for describing large-scale <i>S. aureus</i> genomic epidemiology.

SALTy typing is highly accurate and can quickly analyse large volumes of <i>S. aureus </i> illumina reads (fastq.gz) or genome assemblies (fasta).

Now [published]([https://www.medrxiv.org/content/10.1101/2024.01.21.24301506v1](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001250)) in Microbial Genomics!


---
# Dependencies
1. python (v3.6 or greater)
2. kma (v1.4.9 or greater)
3. pandas  (v1.5.0 or greater)
4. mlst (v2.23.0 or greater)

---

# Installation
```commandline
conda install -c conda-forge -c bioconda salty
```
**See below for Mac M1 installation.

[![Anaconda-Server Badge](https://anaconda.org/bioconda/salty/badges/downloads.svg)](https://anaconda.org/bioconda/salty)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/salty/badges/version.svg)](https://anaconda.org/bioconda/salty)

# Usage
```commandline
salty.py -i <input_folder> -o <output_folder>

GENERAL:
  -t THREADS, --threads THREADS
                        Number of threads (speeds up parsing raw reads).
  -f, --force           Overwite existing output folder.
  --report              Only generate summary report from previous SALTy
                        outputs.
  --check               Checks that required dependencies are available.
  -v, --version         Returns current SaLTy version.


INPUT:
  -i INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Folder of genomes (*.fasta or *.fna) and/or pair end
                        reads (each accession must have *_1.fastq.qz and
                        *_2.fastq.

OUTPUT:
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output Folder to save result.
  -c, --csv_format    Output file in csv format.
  -s, --summary         Concatenate all output assignments into single file.

DB & Program Paths:
  -l LINEAGES, --lineages LINEAGES
                        Path to specific alleles for each lineage.
  -k KMA_INDEX, --kma_index KMA_INDEX
                        Path to indexed KMA database.
  -m, --mlstPrediction  Explained in ReadMe. Used as backup when lineage is unable to be called
                        through SaLTy screening. Marked with *.
```


# Example Usage
Run SaLTy on a folder containing pairs of fastq files and/or assemblies.

```commandline
salty --input_folder /input/reads/folder/ --output_folder /output/folder/name/ --threads INT
```
Specify SaLTy output assigned lineages in comma-separated (.csv) format. 
```commandline
salty --input_folder /input/reads/folder/ --output_folder /output/folder/name/ --csv_format
```

# Example Output
````
Genome		Lineage		SACOL0451	SACOL1908	SACOL2725
SRR992071338	15		24		20		24
SRR992044718	15		24		20		24
SRR10591328	15		24		20		69
````

## Column Descriptions
| Column | Description |
| ----------- | ----------- |
|Sample|Input strain ID (extracted from input files).|
|Lineage|The lineage assigned by the SaLTy algorithim.|
|SACOL0451|The identified allele for the SACOL0451 locus.|
|SACOL1908|The identified allele for the SACOL1908 locus.|
|SACOL2725|The identified allele for the SACOL2725 locus.|

# Installing SaLTy with Conda on Mac M1 Installation (Aarch64)
As of March 2023 not all required SaLTy dependencies have been developed for both Intel and Mac M1 chipsets. 
The conda environments created on Mac M1 can be configured to collect dependencies developed only for Intel (osx-64).
Additionally, python v3.9.0 should be used to avoid compatability issues with Numpy (required for Pandas which is required for SaLTy).

1. Create a conda environment and configure the environment for Intel (osx-64).
```commandline
conda create -n salty
conda activate salty
conda config --env --set subdir osx-64
```

2. Install SaLTy with Python v3.9.0.
```commandline
conda install -c bioconda salty python=3.9.0
```

# SaLTy Prediction with Multi-Locus Sequence Types (Marked with an *)
It is possible to infer a SaLTy lineage through the Multi-Locus Sequencing Type (MLST). During the development of SaLTy a select number of MLST types were associated with SaLTy lineages. This list of MLST types is not exhaustive and will not be continually updated.


MLST type can be used in some instances to infer the SaLTy lineage. Referred to as mlstPrediction (in the SaLTy usage), when SaLTy has attempted to use MLST to predict a lineage an asterisk (*) is marked.


Below are three cases of SaLTy analysis and the use of mlstPrediction is explained.

#1 SaLTy predicted lineage only using three-gene markers (no asterisk).

#2 SaLTy predicted lineage based on MLST type. Three-gene markers unable to infer lineage. Instead associated MLST type used.

#3 SaLTy unable to predict lineage. Tried both three-gene markers and mlst prediction (marker by asterisk).

````
# Genome        Lineage   SACOL0451 SACOL1908 SACOL2725
1 SRR9920718    15        20        24        24
2 ERR109478     *4        13	    -	      16
3 ERR1213758    *No lin.  -         -         -
````
