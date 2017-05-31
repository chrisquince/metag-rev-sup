# metag-rev-sup

# Supplementary materials for Quince et al. Metagenome analysis Review

<a name="assembly"/>
## Assembly based metagenomics analysis

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


### Co-assembly

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
nohup megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 96 -o Assembly > megahit.out&
```

We can have a look at how good the assembly was:
```
contig-stats.pl < Assembly/final.contigs.fa
```
Output should be similar too:
```
sequence #: 863051	total length: 942845899	max length: 490055	N50: 1959	N90: 397
```
Not bad, a N50 of 1959bp. All scripts are available in this repo in the scripts dir. They will 
need to be added to your path though.

We will now perform CONCOCT binning of these contigs. As explained in 
[Alneberg et al.](http://www.nature.com/nmeth/journal/v11/n11/full/nmeth.3103.html) 
there are good reasons to cut up contigs prior to binning. We will use a script from 
CONCOCT to do this. For convenience we 
will create an environmental variables that points to the CONCOCT install directory 
and the scripts directory of the repo change as appropriate:
```
export CONCOCT=~/Installed/CONCOCT
export METAG=/mnt/data-chris/chris/Projects/metag-rev-sup/scripts
```

### Mapping

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
python $METAG/Lengths.py -i contigs/final_contigs_c10K.fa > contigs/final_contigs_c10K.len
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
$METAG/Collate.pl Map | tr "," "\t" > Coverage.tsv
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
In our clustering we obtained 272.

### Annotate genes on contigs

First we call genes on the contigs using prodigal.
```
cd ..
mkdir Annotate
cd Annotate/
$METAG/LengthFilter.pl ../contigs/final_contigs_c10K.fa 1000 > final_contigs_gt1000_c10K.fa
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

### Evaluate clustering based on SCGs

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
This enables us to estimate the completeness of each CONCOCT cluster, in total we found 80 clusters which 
were at least 75% complete and pure in this data set. We will consider these clusters to be 
partially complete metagenome assembled genomes or MAGs. 

![Single copy core gene plot](./Figures/clustering_gt1000_scg.pdf)

### Clusters with differential abundance between Crohn's and healthy children

To explore differences in community composition between the CD and control children, we calculated per sample 
normalised abundances for each cluster. Move into the Concoct directory and using a script from the CONCOCT repo:

```
cd ../Concoct
tr "\t" "," < Coverage.tsv > Coverage.csv
python $METAG/ClusterMeanCov.py Coverage.csv clustering_gt1000.csv ../contigs/final_contigs_c10K.fa > cluster_freq.csv
```

An NMDS plot of these abundances are shown in the figure below:

![NMDS](./Figures/NMDS.pdf)

This can be generated by first parsing the cluster frequency file a little bit:
```
sed 's/Map\///g' cluster_freq.csv | sed 's/^/D/' > cluster_freqR.csv
```
A copy of this same file is in the Results directory and then source the R commands in the file scripts/NMDS.R.
From this it is apparent that there are significant differences in community composition between the two types. In fact 11.2% of the 
variance in community composition was explained by type (perm. ANOVA, p-value = 0.002). 
There was also a higher variance in community composition in the CD children (p-value = 0.0002525). 

To determine those clusters responsible for this difference we performed Kruskal-Wallis non-parametric ANOVA on the 
log-transformed abundances in each group. 
Benjamini-hochberg was used to correct for multiple comparisons. 
This revealed 70 clusters with a q-value < 0.05 of which 26 were at least 75% 
complete MAGs. The majority of these were negatively associated with Crohn's disease. 


