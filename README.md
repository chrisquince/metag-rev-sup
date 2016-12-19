# metag-rev-sup

#Supplementary materials for Quince et al. Metagenome analysis Review

<a name="assembly"/>
##Assembly based metagenomics analysis

This analysis regenerates the results in the left panel of Figure 3. The results will not 
be identical due to some differences in software principally the switch from idba_ud to megahit
for the assembly step. This was changed for reasons of efficiency.

We took twenty samples from healthy control children and those with Crohn's disease from a 
larger study looking at the effect of a treatment (EEN) on the Crohn's microbiota 
[Quince et al. 2015] (https://www.ncbi.nlm.nih.gov/pubmed/26526081).
The sequences comprised a combination of MiSeq 2X250bp and HiSeq 2X150bp paired end reads. The 
sequences used in this tutorial can be downloaded [here](https://metagexample.s3.climb.ac.uk/Reads.tar.gz).

Reads were trimmed and filtered with sickle to remove sequencing adaptors and 
regions with an average quality of < 20. Human DNA was then removed using DeconSeq. 
The minimum and maximum read numbers following trimming and filtering 
were 2,456,000 and 14,720,000 respectively with a median of 7,797,000.

We are now going to perform a basic assembly based metagenomics analysis of these same samples. This will involve 
a collection of different software programs:

1. [megahit](https://github.com/voutcn/megahit): A highly efficient metagenomics assembler currently our default for most studies

2. [bwa](https://github.com/lh3/bwa): Necessary for mapping reads onto contigs

3. [samtools] (http://www.htslib.org/download/): Utilities for processing mapped files

4. [CONCOCT](https://github.com/BinPro/CONCOCT): Our own contig binning algorithm

5. [prodigal] (https://github.com/hyattpd/prodigal/releases/): Used for calling genes on contigs

6. [gnu parallel] (http://www.gnu.org/software/parallel/): Used for parallelising rps-blast

7. [standalone blast] (http://www.ncbi.nlm.nih.gov/books/NBK52640/): Needs rps-blast

8. [COG RPS database] (ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/): Cog databases

9. [GFF python parser] (https://github.com/chapmanb/bcbb/tree/master/gff)


###Co-assembly

Create an example directory to work in:
```
mkdir ~/Example
cd ~/Example
```

Begin by downloading the reads:
```
wget https://metagexample.s3.climb.ac.uk/Reads.tar.gz
tar -xvzf Reads.tar.gz
```
We will then perform a co-assembly of these samples using megahit. First we 
need to get a comma separated list of forward and reverse read files:
```
ls Reads/*R1*fasta | tr "\n" "," | sed 's/,$//' > R1.csv
ls Reads/*R2*fasta | tr "\n" "," | sed 's/,$//' > R2.csv
```
Then set assembler running may take several hours...
```
nohup megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 96 -o Assembly --presets meta > megahit.out&
```

We can have a look at how good the assembly was:
```
contig-stats.pl < Assembly/final.contigs.fa
```
Output should be similar too:
```
sequence #: 1316028	total length: 964326614	max length: 318273	N50: 1088	N90: 299
```
Not bad, a N50 of 1088bp. All scripts are available in this repo in the scripts dir. They will 
need to be added to your path though.

We will now perform CONCOCT binning of these contigs. As explained in 
[Alneberg et al.](http://www.nature.com/nmeth/journal/v11/n11/full/nmeth.3103.html) 
there are good reasons to cut up contigs prior to binning. We will use a script from 
CONCOCT to do this. For convenience we 
will create an environmental variables that points to the CONCOCT install directory 
change as appropriate:
```
export CONCOCT=~/Installed/CONCOCT
```

###Mapping

First we cut up contigs and place in new dir:

```bash
mkdir contigs
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m Assembly/final.contigs.fa > contigs/final_contigs_c10K.fa
```

Having cut-up the contigs the next step is to map all the reads from each sample back onto them. 
First index the contigs with bwa:

```bash
cd contigs
bwa index final_contigs_c10K.fa
cd ..
```

Then perform the actual mapping you may want to put this in a shell script and adjust the number of threads 
according to your resources:

```bash
mkdir Map

for file in Reads/*R1.fasta
do
    stub=${file%_R1.fasta}
    stub=${stub#Reads\/}
    
    bwa mem -t 64 contigs/final_contigs_c10K.fa $file Reads/${stub}_R2.fasta > Map/${stub}.sam
    echo $stub
done
```

Then we need to calculate our contig lengths.

```bash
python ~/bin/Lengths.py -i contigs/final_contigs_c10K.fa > contigs/final_contigs_c10K.len
```

Then we calculate coverages for each contig in each sample:

```bash
for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub	
    samtools view -h -b -S $file > ${stub}.bam 
    samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam
    samtools sort ${stub}.mapped.bam -o ${stub}.mapped.sorted
    bedtools genomecov -ibam ${stub}.mapped.sorted -g contigs/final_contigs_c10K.len > ${stub}_cov.txt
done
```

and use awk to aggregate the output of bedtools:

```bash
for i in Map/*_cov.txt 
do 
   echo $i
   stub=${i%_cov.txt}
   stub=${stub#Map\/}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > Map/${stub}_cov.csv
done
```

and finally run the following perl script to collate the coverages across samples, where we have simply adjusted the format 
from csv to tsv to be compatible with CONCOCT:

```bash
Collate.pl Map | tr "," "\t" > Coverage.tsv
```


We will only run concoct with default settings here. The only one we 
vary is the cluster number which should be at least twice the number of genomes in your 
co-assembly. 

```
mkdir Concoct
cd Concoct
mv ../Coverage.tsv .
concoct --coverage_file Coverage.tsv --composition_file ../contigs/final_contigs_c10K.fa 
cd ..
```

The output of concoct is principally an assignment of each contig to a bin as a comma separated file, 
clustering_gt1000.csv, with a simple format of contig,assignment:
```
tail clustering_gt1000.csv
```
We can count the total number of bins as follows:
```
cut -d"," -f2 clustering_gt1000.csv | sort | uniq -c | wc
```
We should give something like 264.

###Annotate genes on contigs

First we call genes on the contigs using prodigal.
```
mkdir Annotate
cd Annotate/
LengthFilter.pl ../contigs/final_contigs_c10K.fa 1000 > final_contigs_gt1000_c10K.fa
prodigal -i final_contigs_gt1000_c10K.fa -a final_contigs_gt1000_c10K.faa -d final_contigs_gt1000_c10K.fna  -f gff -p meta -o final_contigs_gt1000_c10K.gff
```

Then we assign COGs with RPSBlast. The script in the CONCOCT distro for doing this requires 
an environment variable pointing to the COG database to be set. We also set another variable 
pointing to the CONCOCT scripts for convenience:
```
export COGSDB_DIR=/home/opt/rpsblast_db
export CONCOCT=/home/chris/Installed/CONCOCT
```

Then we run this on 32 cores:
```
$CONCOCT/scripts/RPSBLAST.sh -f final_contigs_gt1000_c10K.faa -p -c 32 -r 1
```

This script performs an RPS-BLAST of translated 
sequences against the NCBI COG database using an e-value cut-off of 1.0e-3. 
Each query was assigned to the top RPS-BLAST and only if it covered at least 50% of the target sequence. 

###Evaluate clustering based on SCGs

The script COG_table.py was then used to generate a table of counts for 36 COGs that we 
previously identified as being found in all bacterial genome with a single copy. We refer to these 
COGs as single-copy core genes (SCGs): 
```
$CONCOCT/scripts/COG_table.py -b final_contigs_gt1000_c10K.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c ../Concoct/clustering_gt1000.csv --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_gt1000_scg.tsv
```

and a pdf:
```
$CONCOCT/scripts/COGPlot.R -s clustering_gt1000_scg.tsv -o clustering_gt1000_scg.pdf
```
This enables us to estimate the completeness of each CONCOCT cluster, in total we found 88 clusters which 
were at least 75% complete and pure in this data set. We will consider these clusters to be 
partially complete metagenome assembled genomes or MAGs. In total we obtained 54 clusters 
that were 75% pure and complete.

![Single copy core gene plot](./Figures/clustering_gt1000_scg.pdf)

### Taxonomic classification of contigs

There are many ways to taxonomically classify assembled sequence. We suggest a gene based approach. 
The first step is 
to call genes on all contigs that are greater than 1,000 bp. Shorter sequences are unlikely to contain complete 
coding sequences. The following requires that you have a Diamond formatted version of the NCBI NR on your system. 
To ensure compatibility with the files below this can be downloaded by:

```
wget http://nrdatabase.s3.climb.ac.uk/nr.dmnd
```

Set the environment variable NR_DMD to point to the location of this file:
```
export NR_DMD=$HOME/native/Databases/nr/FASTA/nr.dmnd
```

Return to the example directory and make a new directory...
```
cd ..
mkdir AssignTaxa
cd AssignTaxa
cp ../Annotate/final_contigs_gt1000_c10K.faa .
diamond blastp -p 32 -d $NR_DMD -q final_contigs_gt1000_c10K.faa -a final_contigs_gt1000_c10K > d.out
diamond view -a final_contigs_gt1000_c10K.daa -o final_contigs_gt1000_c10K_nr.m8
```
To classify the contigs we need two files a gid to taxid mapping file and a mapping of taxaid to full lineage:

1. gi_taxid_prot.dmp

2. all_taxa_lineage_notnone.tsv

These can also be downloaded from the Dropbox:
``` 
wget https://www.dropbox.com/s/x4s50f813ok4tqt/gi_taxid_prot.dmp.gz
wget https://www.dropbox.com/s/honc1j5g7wli3zv/all_taxa_lineage_notnone.tsv.gz
```

To perform the classification we use the ClassifyContigNR.py script which is included in this repo.
The path to these files are default in the ClassifyContigNR.py script as the variables:
```
DEF_DMP_FILE = "/home/chris/native/Databases/nr/FASTA/gi_taxid_prot.dmp"

DEF_LINE_FILE = "/home/chris/native/Databases/nr/FASTA/all_taxa_lineage_notnone.tsv"
```

We calculate the gene length in amino acids before running this.
Then we can assign the contigs and genes called on them:
```
python Lengths.py -i final_contigs_gt1000_c10K.faa > final_contigs_gt1000_c10K.len
python ClassifyContigNR.py final_contigs_gt1000_c10K_nr.m8 final_contigs_gt1000_c10K.len -o final_contigs_gt1000_c10K_nr -l /mypath/all_taxa_lineage_notnone.tsv -g /mypath/gi_taxid_prot.dmp
```


