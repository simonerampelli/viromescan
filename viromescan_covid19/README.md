# viromescan_covid19 
This is the module implemented by the ViromeScan Researchers for identified viral sequences of SARS-CoV-2 in metagenomic data.
Simply copy and paste the contents of this folder in your viromescan folder (pay attention to put the provided materials in the correspondent location of your viromescan folder without substituting the original contents)

## Instruction after download:
```
unzip database/bowtie2.zip
cp database/bowtie2/* path_to_your_viromescan_folder/database/bowtie2/
cp var/* path_to_your_viromescan_folder/var/
cp viromescan_covid19.sh path_to_your_viromescan_folder/
```

Create the new link typing:
```
ls -s path_to_your_viromescan_folder/viromescan_covid19.sh your_bin_folder/viromescan_covid19
```

Visualize help:
```
viromescan_covid19
```
