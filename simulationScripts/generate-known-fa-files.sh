# This script builds the transcriptome annotations from the reference annotations for all methods.
# You should have already downloaded in Bowtie2Index directory (read Bowtie2Index/README file for this). 
# Required threads: 4

# First unzip gtf files:
unzip genes.gtf.zip
unzip genesNew.gtf.zip
# create transcriptome_data directory with tophat
tophat -o A1.tophat.out --no-novel-juncs --no-discordant --no-mixed -x 100 -T -p 4 -G genes.gtf --transcriptome-index=transcriptome_data/known Bowtie2Index/genome testFile_1.fastq testFile_2.fastq
rm -r A1.tophat.out
cd transcriptome_data
rsem-prepare-reference --bowtie2 known.fa ref #this is for rsem
sed -e 's/\,.*//' known.fa > output.fasta #this is for making tigar names work
sailfish index -t known.fa -o indDir -p 4 #this is for sailfish



