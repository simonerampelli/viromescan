# viromescan
ViromeScan allows the user to explore and taxonomically characterize the virome from metagenomic reads, efficiently denoising samples from reads of other microorganisms.

A new module is now included for detection of SARS-CoV-2 virus.

Below the easier way for installing Viromescan.

Please refer to [this publication](https://doi.org/10.1186/s12864-016-2446-3) for further information.

## README
SPACE REQUIREMENTS : ViromeScan requires about 45 GB of free space to work on your device

### FIRST STEPS
Install conda and clone the git folder
```
git clone https://github.com/simonerampelli/viromescan.git
mv viromescan viromescan2
```

### 1) CREATE CONDA ENVIRONMENT FROM SOURCE
```
conda env create -f viromescan2.yml -p /your/conda/path/envs/viromescan
```

### 2) DOWNLOAD AND UNTAR VIROMESCAN
Download the native versione of ViromeScan and its database from [this site](https://sourceforge.net/projects/viromescan/files/)
Then unzip and untar the folder.
```
wget https://sourceforge.net/projects/viromescan/files/viromescan.tar.gz 
tar -zxvf viromescan.tar.gz
```

### 3) SUBSTITUTE THE NEW VERSION OF bmtagger.sh PROVIDED HERE IN GITHUB TO THE OLD VERSION WITHIN THE VIROMESCAN FOLDER
```
rm -fr viromescan/tools/bmtagger.sh
mv bmtagger.sh viromescan/tools/
```

### 4) THE DATABASES
```
cd $viromescan_path/viromescan/database

gzip -d Bacteria_custom/*
gzip -d  bowtie2/*
gzip -d hg19/*

cd hg19/
```
Now you need to create the indexes database for hg19.
*N.B. Bmtagger scripts (bmfilter, srprism) require about 8.5Gb RAM memory and three times as much hard-disk space for index data.*

-  Make indexes for bmfilter. 
```
bmtool -d hg19reference.fa -o hg19reference.bitmask -A 0 -w 18
```
- Make index for srprism
```
srprism mkindex -i hg19reference.fa -o hg19reference.srprism -M 7168
```
- Make blastdb for blast
```
makeblastdb -in hg19reference.fa -dbtype nucl
```

### 5) MOVE THE VIROMESCAN FOLDER IN THE CONDA DIRECTORY
```
mv viromescan /your/conda/path/envs/viromescan
```

### 6) OPTIONAL: INSTALL THE COVID19 MODULE
Follow the instructions contained within the viromescan_covid19 folder.


### 7) CREATE A VIRTUAL LINK FOR THE viromescan.sh AND viromescan_covid19.sh SCRIPT
```
ln -s /your/conda/path/envs/viromescan/viromescan/viromescan.sh  /your/conda/path/envs/viromescan/bin/viromescan
ln -s /your/conda/path/envs/viromescan/viromescan/viromescan_covid19.sh  /your/conda/path/envs/viromescan/bin/viromescan_covid19
```

### 8) ACTIVATE AND DEACTIVATE THE VIROMESCAN ENVIRONMENT
```
conda activate viromescan
conda deactivate viromescan
```

### 9) USAGE AND HELP 

For usage and help information please digit "viromescan" or "viromescan_covid19" without any option on your command line


### 10) EXAMPLES OF USAGE

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
