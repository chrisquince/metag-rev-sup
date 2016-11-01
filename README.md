# metag-rev-sup

#Supplementary materials for Quince et al. Metagenome analysis Review

<a name="assembly"/>
##Assembly based metagenomics analysis

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

Begin by downloading the reads:
```
wget ****
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
contig-stats.pl < Coassembly/final.contigs.fa
```
Not great, but expected given the low read number.

We will now perform CONCOCT binning of these contigs. As explained in [Alneberg et al.](http://www.nature.com/nmeth/journal/v11/n11/full/nmeth.3103.html) 
there are good reasons to cut up contigs prior to binning. We will use a script from CONCOCT to do this. For convenience we 
have created an environmental variables that points to the CONCOCT install directory:
```
echo $CONCOCT
```

###Mapping

First we cut up contigs and place in new dir:

```bash
mkdir contigs
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m Coassembly/final.contigs.fa > contigs/final_contigs_c10K.fa
```

Having cut-up the contigs the next step is to map all the reads from each sample back onto them. First index the contigs with bwa:

```bash
cd contigs
bwa index final_contigs_c10K.fa
cd ..
```

Then perform the actual mapping you may want to put this in a shell script:

```bash
mkdir Map

for file in MetaTutorial/{C,H}*R12.fasta
do 
   
   stub=${file%_R12.fasta}
   stub=${stub#MetaTutorial\/}
   echo $stub

   bwa mem -t 8 contigs/final_contigs_c10K.fa $file > Map/${stub}.sam
done
```
