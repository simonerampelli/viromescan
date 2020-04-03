# viromescan
ViromeScan allows the user to explore and taxonomically characterize the virome from metagenomic reads, efficiently denoising samples from reads of other microorganisms.

A new module is now included for detection of SARS-CoV-2 virus.

Below the easier way for installing Viromescan.

Please refer to [this publication](https://doi.org/10.1186/s12864-016-2446-3) for further information.

## README
SPACE REQUIREMENTS : ViromeScan requires about 45 GB of free space to work on your device

### FIRST STEPS
Install conda
Download the viromescan2.tar.gz folder and the viromescan2.yml file 

### 1) CREATE CONDA ENVIRONMENT FROM SOURCE
```
conda env create -f viromescan2.yml -p /your/conda/path/envs/viromescan
```

### 2) UNTAR VIROMESCAN
```
tar -zxvf viromescan2.tar.gz
```

### 2) MOVE THE VIROMESCAN FOLDER IN THE CONDA DIRECTORY
```
mv viromescan /your/conda/path/envs/viromescan
```

### 5) CREATE A VIRTUAL LINK FOR THE viromescan.sh AND viromescan_covid19.sh SCRIPT
```
ln -s /your/conda/path/envs/viromescan/viromescan/viromescan.sh  /your/conda/path/envs/viromescan/bin/viromescan
ln -s /your/conda/path/envs/viromescan/viromescan/viromescan_covid19.sh  /your/conda/path/envs/viromescan/bin/viromescan_covid19
```

### 6) ACTIVATE AND DEACTIVATE THE VIROMESCAN ENVIRONMENT
```
conda activate viromescan
conda deactivate viromescan
```

### 7) USAGE AND HELP 

For usage and help information please digit "viromescan" or "viromescan_covid19" without any option on your command line


### 8) EXAMPLES OF USAGE

A) Computing the analysis on a single-end .fastq file for the human DNA viruses
```
viromescan -p 3 -m /your/conda/path/envs/viromescan -d human_DNA -1 sample_A.fastq -o viromescan_sample_A 
```
B) Computing the analysis on paired-end .fastq files for the global viruses
```
viromescan -p 3 -m /your/conda/path/envs/viromescan -d virus_ALL -1 sample_A_r1.fastq -2 sample_A_r2.fastq -o viromescan_sample_A
```
C) Computing the analysis using the COVID19 module for SARS-CoV-2 virus detection
```
viromescan -p 3 -m /your/conda/path/envs/viromescan -1 sample_A_r1.fastq -2 sample_A_r2.fastq -o viromescan_sample_A_covid19
```
