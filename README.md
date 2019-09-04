# viromescan
ViromeScan allows the user to explore and taxonomically characterize the virome from metagenomic reads, efficiently denoising samples from reads of other microorganisms.

Please see also the [SourceForge web page](https://sourceforge.net/projects/viromescan).

## README
SPACE REQUIREMENTS : ViromeScan requires about 45 GB of free space to work on your device

### FIRST STEPS
After downloading the tool, in order to perform the analysis, you need to unzip the downloaded file and compile the database. 

### 1) UNTAR VIROMESCAN
```
tar -zxvf viromescan.tar.gz
```
### 2) MOVE THE DATABASE DIRECTORY IN THE VIROMESCAN FOLDER
```
cd $viromescan_path/viromescan/database
```
### 3) UNZIP THE ZIPPED FILES
```
gzip -d Bacteria_custom/*
gzip -d  bowtie2/*
gzip -d hg19/*
```
### 4) BUILD THE HUMAN DATABASE FOR THE FILTERING PROCEDURE
```
cd hg19/
```
Now you need to create the indexes database for running bmtagger. bmtagger and other bmtools necessary to run viromescan are already present in the $PWD/viromescan/tools folder. You can eventually download a compatible version of bmtools for your operating system at ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger.

*N.B. Bmtagger scripts (bmfilter, srprism) require about 8.5Gb RAM memory and three times as much hard-disk space for index data.*
  
If you downloaded Bmtools, you will have to move bmtagger.sh, bmfilter, bmtool, extract_fullseq and srprism scripts in the directory $PWD/viromescan/tools inside the viromescan folder.*

Now you need to:
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

### 5) CREATE A VIRTUAL LINK OF THE viromescan.sh SCRIPT

On your command line type:
```
ln -s $viromescan_path/viromescan/viromescan.sh  $your_binary_folder/viromescan
```

### 6) USAGE AND HELP 

For usage and help information please digit "viromescan" without any option on your command line


### 7) EXAMPLES OF USAGE

A) Computing the analysis on a single-end .fastq file for the human DNA viruses
```
viromescan -p 3 -m /user/matteo/ -d human_DNA -1 /user/matteo/seqs/sample_A.fastq -o viromescan_sample_A 
```
B) Computing the analysis on paired-end .fastq files for the global viruses
```
viromescan -p 3 -m /user/matteo/ -d virus_ALL -1 /user/matteo/seqs/sample_A_r1.fastq -2 /user/matteo/seqs/sample_A_r2.fastq -o viromescan_sample_A
```
