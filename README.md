# Timo Variant Calling Pipeline

Identifying minority variants in viral sequence data can highlight regions of the genome 
that are under selection or indicate regions with increased mutational tolerance. 
Additonally, investigation of the minority variants may be able to inform transmission 
chains or detect subtle shifts in the viral population before consensus changes occur. 

Confident prediction of minority variants in complicated by errors introduced during the 
amplification pre-processing steps required to sequence many viral genomes, 
and also can be introduced during the sequencing process itself. Many software packages exist to 
identify minority variants and each package differs slightly in its bioinformatic and statistical approach.

The timo script takes de-duplicated bam files and identifies minority variants above user given 
coverage and allele frequency cutoffs, outputting user friendly csv files.

We also provide two additional scripts to generate consensus sequences and to add amino acid sequences 
to the output files from timo.

#### Requirements

python3  
numpy  
pysam  
SciPy  

## Running timo:

Clone the git repository and edit the following items

#### **run_IndexVariants.sh**

The run_IndexVariants.sh script provides run information and calls the index_variants.sh script described below. Once the pipeline is set up on your system, the run_IndexVariants.sh script is the only script that should need editing

*JOB_NAME:* Give your run an identifiable name for output purposes

*RUNDIR:* Path to the directory where your scripts are located and where you will be running the pipeline from

*BAM_DIR:* Path to the directory where your bam files are located

*STRAIN:* The strain of the virus you are using. We recommend 'cov19' for SARS-CoV-2 variant detection

*REFSEQ:* Path to and file name of your reference sequence. Ex. /data/$USER/SARS-CoV-2.fa

*USER:* email address for running information

*array:* 0-#, zero based number of samples to run the pipeline on. Used for submitting array jobs for faster output.

#### **index_variants.sh**

index_variants.sh is the run script that sorts, indexes, marks Duplicates and de-duplicates the bam files. Once this is complete, the script calls the timo script to identify variants.

While nothing in this file needs to be changed, if you would like to run timo with different paraments than the defaults, you can change this here. 

The *timo* parameters are as follows:



--ref (-r): Full path to the reference file. This is passed in as a variable provided in run_IndexVariants.sh  

--infile (-i): Input single bamfile. Needs full path if not local. This is passed in as a variable provided in run_IndexVariants.sh using the bam directory and the array number. 

--seqment (-s): Input single segment. Only needed in multi-segment genomes 

--qual (-q): Phred quality cutoff, default is 25.

--cutoff (-c): Minor variant frequency cutoff for called variants, default is 0.01 (1%). We recommend passing in a value of 0.001 and filtering for higher frequency variants from the output.

--strain (-T): Strain of the virus being analyzed. This is passed in as a variable provided in run_IndexVariants.sh  

--covercutoff (-C): Coverage cutoff for consensus, default is 200. We recommend passing in a value of 1 and filtering for high coverage regions from the output  

### *timo* output

a csv file with the variant information for each nucleotide position will be stored in a directory called FILES/fullvarlist that *timo* creates in the run directory

### Additional options:  
##### Creating consensus and coverage files and addition of amino acid information to *timo* output

We provide two additional scripts for
1. Creating consensus files and a coverage file
2. Adding amino acid information to the fullvarlist files created from timo

### Running ConsensusFasta.Coverage.v4.py

This code does three different things:
1. Checks the snplist file to make sure all positions are present, and not duplicated
    - if there is missing data (no reads aligning to position), then the position won't be present in snplist file
    - the updated csv files will be saved where the original snplist csv were

2. Generates a consensus sequence using checked/updated snplist file
    - Generates by 'segment'
    - If coverage is below provided cutoff, the nt will be replaced with an "N"

3. Pulls coverage information- that can then be used for coverage figures, etc

4. Will make separate fasta files for seqs

**Input parameters:**

--ref: path to reference, and name of reference fasta file. Should be exact same as used for *timo*

--var: path to the variant files (snplist)

--cov: coverage cutoff for calling a major nt. Default is 10x.

--minfreq: minorfreq cutoff used to call minor variants using the readreport script.
            This will be used to pull out files, and name files. So it must be
            the same as what was used to generate snplist files.

--strain: the strain named that was used to generate the snplist files. For example:
            COV19, H1N1, H3N2. These should match the strain in the name of the snplist files. These are case-sensitive

--savecov: The directory where you want to save your coverage files. Default is
            the location where the python script is running

--savecon: The directory where you want to save your consensus files. Default is
            the location where the python script is running.
            
**Example command for ConsensusFasta.Coverage.v4.py:**

> python3 ./ConsensusFasta.Coverage.v4.py  
    --ref /data/$USER/SARS-CoV-2.fa  
    --var ./FILES/fullvarlist  
    --cov 5  
    --minfreq 0.001  
    --strain COV19  
    --savecov ./FILES/coverage  
    --savecon ./FILES/consensus  

### Running AddAminoGene.5.py

This program adds amino acid and gene information to the snplist files using a 'features' file

We recommend creating a subdirectory in the FILES directory called 'aasnplist'.

Requires the UPDATED snplist files (ie. after running ConsensusFasta.Coverage.v4.y) and csv file with features interested in

You must have a features file for this script: 
The csv file with features that are of interest should have the following columns:
SEGMENT,START,END,NAME

A feature file for SARS-CoV-2 is provided. 

**Input parameters**

--ref: path to reference, and name of reference fasta file. Should be exact same as used for *timo*

--var: path to the variant files (snplist)

--freqcut: minorfreq cutoff used to call minor variants using the readreport script. This will be used to pull out files, and name files. So it must be the same as what was used to generate snplist files.

--strain: the strain named that was used to generate the snplist files. For example: COV19, H1N1, H3N2. These should match the strain in the name of the snplist files. These are case-sensitive.

--features: path to the features file

--savedir: The directory where you want to save your aasnplist files. Default is
            the location where the python script is running

**Example command for ConsensusFasta.Coverage.v4.py:**

> python3 AddAminoGene.5.py   
-- freqcut 0.001  
--ref /data/$USER/SARS-CoV-2.fa  
--var ./FILES/fullvarlist  
--strain COV19  
--features ./sars-cov-2-features.4.csv  
--save_dir ./FILES/aasnplist/. 


#### Contributors:

Tim Song  
Kate Johnson   
David Chen  
Chang Wang  
Allison Roder  


```

```
