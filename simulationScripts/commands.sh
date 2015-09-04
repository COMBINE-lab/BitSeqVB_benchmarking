# Note: For reasonable computing time this script should be split into parallel ones according to the following steps.
# At least 4 cores should be available (tigar needs 12). 

#################################################################################
# STEP 0. Choose dataset							#
# STEP 1. Generate RPK values							#
# STEP 2. Simulate fastq files with spanki.					#
# STEP 3. Read mapping procedure:						#
#	3.a Align reads with bowtie						#
#	3.b Align reads with tophat						#
#	Note: Steps 3.a and 3.b can run in parallel if possible (recommended).	#
# STEP 4. Run methods and post-process necessary output:			#
#	4.a BitSeqMCMC								#
#	4.b BitSeqVB								#
#	4.c Casper								#
#	4.d Cufflinks								#
#	4.e RSEM								#
#	4.f Sailfish								#
#	4.g Tigar2								#
#	4.h eXpress								#
#	Note: Steps 4.a,...,4.h can run in parallel if possible (recommended).	#
# STEP 5. Produce graphs							#
#################################################################################


#################################################################################
# STEP 0. Choose dataset by uncommenting one of the following four lines:       #
# cd random-scenario								#
# cd bitseq-scenario								#
# cd poisson-glm-mix-scenario 							#
# cd rsem-scenario								#
#################################################################################

# Now your working directory should be one of the previous lines. 

# STEP 1: Generate rpk values
R CMD BATCH trFiles.R


# STEP 2. Simulate fastq files with spanki.					
#################################################################################################################################
# NOTE: Each of the following commands can also run in parallel for each replicate (highly recommended).			#
# replicate 1															#
spankisim_transcripts -o A1 -g ../genesNew.gtf -f ../Bowtie2Index/genome.fa -bp 76 -t tr_File_A1.tr -ends 2			#
cd A1 && R CMD BATCH ../../getTranscriptNames.R											#
./../../shuffleReadsTMP.sh 		#shuffle read pairs (essential if your are  using sailfish and/or eXpress)		#
cd ..																#
# replicate 2															#
spankisim_transcripts -o A2 -g ../genesNew.gtf -f ../Bowtie2Index/genome.fa -bp 76 -t tr_File_A2.tr -ends 2			#
cd A2 && R CMD BATCH ../../getTranscriptNames.R											#
./../../shuffleReadsTMP.sh 		#shuffle read pairs (essential  if your are using sailfish and/or eXpress)		#
cd ..																#
# replicate 3															#
spankisim_transcripts -o A3 -g ../genesNew.gtf -f ../Bowtie2Index/genome.fa -bp 76 -t tr_File_A3.tr -ends 2			#
cd A3 && R CMD BATCH ../../getTranscriptNames.R											#
./../../shuffleReadsTMP.sh 		#shuffle read pairs (essential  if your are using sailfish and/or eXpress)		#
cd ..																#
# replicate 4															#
spankisim_transcripts -o A4 -g ../genesNew.gtf -f ../Bowtie2Index/genome.fa -bp 76 -t tr_File_A4.tr -ends 2			#
cd A4 && R CMD BATCH ../../getTranscriptNames.R											#
./../../shuffleReadsTMP.sh 		#shuffle read pairs (essential  if your are using sailfish and/or eXpress)		#
cd ..																#
# replicate 5															#
spankisim_transcripts -o A5 -g ../genesNew.gtf -f ../Bowtie2Index/genome.fa -bp 76 -t tr_File_A5.tr -ends 2			#
cd A5 && R CMD BATCH ../../getTranscriptNames.R											#
./../../shuffleReadsTMP.sh 		#shuffle read pairs (essential  if your are using sailfish and/or eXpress)		#
cd ..																#
#################################################################################################################################



# STEP 3.a Bowtie						
#################################################################################################################################################
mkdir bowtie																	#
cd bowtie																	#
# Run Bowtie and compute alignment probabilities for BitSeq. 											#
# NOTE: Each of the following commands can also run in parallel for each replicate (highly recommended).					# 
# replicate 1																	#
bowtie2 -q -k 100 --threads 4 --no-mixed --no-discordant -x ../../transcriptome_data/known -1 ../A1/out1.fq -2 ../A1/out2.fq -S A1.sam		#
parseAlignment A1.sam -o A1.prob --trSeqFile ../../transcriptome_data/known.fa --trInfoFile data.tr --uniform --verbose				#
# replicate 2																	#
bowtie2 -q -k 100 --threads 4 --no-mixed --no-discordant -x ../../transcriptome_data/known -1 ../A2/out1.fq -2 ../A2/out2.fq -S A2.sam		#
parseAlignment A2.sam -o A2.prob --trSeqFile ../../transcriptome_data/known.fa --uniform --verbose						#
# replicate 3																	#
bowtie2 -q -k 100 --threads 4 --no-mixed --no-discordant -x ../../transcriptome_data/known -1 ../A3/out1.fq -2 ../A3/out2.fq -S A3.sam		#
parseAlignment A3.sam -o A3.prob --trSeqFile ../../transcriptome_data/known.fa --uniform --verbose						#
# replicate 4																	#
bowtie2 -q -k 100 --threads 4 --no-mixed --no-discordant -x ../../transcriptome_data/known -1 ../A4/out1.fq -2 ../A4/out2.fq -S A4.sam		#
parseAlignment A4.sam -o A4.prob --trSeqFile ../../transcriptome_data/known.fa --uniform --verbose						#
# replicate 5																	#
bowtie2 -q -k 100 --threads 4 --no-mixed --no-discordant -x ../../transcriptome_data/known -1 ../A5/out1.fq -2 ../A5/out2.fq -S A5.sam		#
parseAlignment A5.sam -o A5.prob --trSeqFile ../../transcriptome_data/known.fa --uniform --verbose						#
cd ..																		#
#################################################################################################################################################


# STEP 3.b Tophat
#########################################################################################################################################################
# Run Tophat and sort bam files for casper.														#
mkdir tophat																		#
cd tophat																		#
# NOTE: Each of the following commands can also run in parallel for each replicate (highly recommended).						# 
#replicate 1																		#
tophat -o A1.tophat.out -p 4 -g 100 --read-mismatches 10 --read-gap-length 10 --read-edit-dist 10 ../../Bowtie2Index/genome ../A1/sim_1.fastq ../A1/sim_2.fastq																		#
cd A1.tophat.out																	#
samtools sort accepted_hits.bam accepted_hits.sorted													#
samtools index accepted_hits.sorted.bam															#
cd ..																			#
#replicate 2																		#
tophat -o A2.tophat.out -p 4 -g 100 --read-mismatches 10 --read-gap-length 10 --read-edit-dist 10 ../../Bowtie2Index/genome ../A2/sim_1.fastq ../A2/sim_2.fastq																		#
cd A2.tophat.out																	#
samtools sort accepted_hits.bam accepted_hits.sorted													#
samtools index accepted_hits.sorted.bam															#
cd ..																			#
#replicate 3																		#
tophat -o A3.tophat.out -p 4 -g 100 --read-mismatches 10 --read-gap-length 10 --read-edit-dist 10 ../../Bowtie2Index/genome ../A3/sim_1.fastq ../A3/sim_2.fastq																		#
cd A3.tophat.out																	#
samtools sort accepted_hits.bam accepted_hits.sorted													#
samtools index accepted_hits.sorted.bam															#
cd ..																			#
#replicate 4																		#
tophat -o A4.tophat.out -p 4 -g 100 --read-mismatches 10 --read-gap-length 10 --read-edit-dist 10 ../../Bowtie2Index/genome ../A4/sim_1.fastq ../A4/sim_2.fastq																		#
cd A4.tophat.out																	#
samtools sort accepted_hits.bam accepted_hits.sorted													#
samtools index accepted_hits.sorted.bam															#
cd ..																			#
#replicate 5																		#
tophat -o A5.tophat.out -p 4 -g 100 --read-mismatches 10 --read-gap-length 10 --read-edit-dist 10 ../../Bowtie2Index/genome ../A5/sim_1.fastq ../A5/sim_2.fastq																		#
cd A5.tophat.out																	#
samtools sort accepted_hits.bam accepted_hits.sorted													#
samtools index accepted_hits.sorted.bam															#
cd ../..																		#
#########################################################################################################################################################


# STEP 4.a BitSeqMCMC
#################################################################################################################		
# Run BitSeqMCMC												#
mkdir bitseq													#
cd bitseq													#
# NOTE: Each of the following commands can also run in parallel for each replicate (highly recommended).	#
estimateExpression ../bowtie/A1.prob -o A1 -p A1.txt -P 4  #replicate 1						#
estimateExpression ../bowtie/A2.prob -o A2 -p A2.txt -P 4  #replicate 2						#
estimateExpression ../bowtie/A3.prob -o A3 -p A3.txt -P 4  #replicate 3						#
estimateExpression ../bowtie/A4.prob -o A4 -p A4.txt -P 4  #replicate 4						#
estimateExpression ../bowtie/A5.prob -o A5 -p A5.txt -P 4  #replicate 5						#
# Process output												#
R CMD BATCH ../../withinGeneBitSeq.R										#
cd ..														#
#################################################################################################################


# STEP 4.b BitSeqVB
#########################################################################
# Run BitSeqVB								#
mkdir vb								#
cd vb									#
estimateVBExpression -o A1 ../bowtie/A1.prob  #replicate 1		#
estimateVBExpression -o A2 ../bowtie/A2.prob  #replicate 2		#
estimateVBExpression -o A3 ../bowtie/A3.prob  #replicate 3		#
estimateVBExpression -o A4 ../bowtie/A4.prob  #replicate 4		#
estimateVBExpression -o A5 ../bowtie/A5.prob  #replicate 5		#
# Process output							#
R CMD BATCH ../../withinGeneVB.R					#
cd ..									#
#########################################################################



# STEP 4.c Casper
#########################################################################
# Run casper								#
mkdir casper								#
cd casper								#
R CMD BATCH ../../casper.R						#
# Process output							#
R CMD BATCH ../../casperNames.R						#
R CMD BATCH ../../withinGeneCasper.R					#
cd ..									#
#########################################################################


# STEP 4.d Cufflinks
#################################################################################################################################
# Run Cufflinks															#
mkdir cufflinks															#
cd cufflinks															#
cufflinks -q -o A1.cufflinks.out -u -p 4 -G ../../genes.gtf ../tophat/A1.tophat.out/accepted_hits.bam >> A1_cuff.log		#
cufflinks -q -o A2.cufflinks.out -u -p 4 -G ../../genes.gtf ../tophat/A2.tophat.out/accepted_hits.bam >> A2_cuff.log		#
cufflinks -q -o A3.cufflinks.out -u -p 4 -G ../../genes.gtf ../tophat/A3.tophat.out/accepted_hits.bam >> A3_cuff.log		#
cufflinks -q -o A4.cufflinks.out -u -p 4 -G ../../genes.gtf ../tophat/A4.tophat.out/accepted_hits.bam >> A4_cuff.log		#
cufflinks -q -o A5.cufflinks.out -u -p 4 -G ../../genes.gtf ../tophat/A5.tophat.out/accepted_hits.bam >> A5_cuff.log		#
# Process output														#
R CMD BATCH ../../cufflinksNames.R												#
R CMD BATCH ../../withinGeneCufflinks.R												#
cd ..																#
#################################################################################################################################


# STEP 4.e RSEM
#################################################################################################################################################################
# Run RSEM																			#
mkdir rsem																			#
cd rsem																				#
# NOTE: Each of the following commands can also run in parallel for each replicate (highly recommended).							#
rsem-calculate-expression -p 4 --paired-end --bowtie2 --no-bam-output ../A1/out1.fq ../A1/out2.fq ../../transcriptome_data/ref expressionA1			#
rsem-calculate-expression -p 4 --paired-end --bowtie2 --no-bam-output ../A2/out1.fq ../A2/out2.fq ../../transcriptome_data/ref expressionA2			#
rsem-calculate-expression -p 4 --paired-end --bowtie2 --no-bam-output ../A3/out1.fq ../A3/out2.fq ../../transcriptome_data/ref expressionA3			#
rsem-calculate-expression -p 4 --paired-end --bowtie2 --no-bam-output ../A4/out1.fq ../A4/out2.fq ../../transcriptome_data/ref expressionA4			#
rsem-calculate-expression -p 4 --paired-end --bowtie2 --no-bam-output ../A5/out1.fq ../A5/out2.fq ../../transcriptome_data/ref expressionA5			#
# uncomment the following lines if you want to use rsem-PME													#
#rsem-calculate-expression -p 4 --calc-pme --paired-end --bowtie2 --no-bam-output ../A1/out1.fq ../A1/out2.fq ../../transcriptome_data/ref expressionA1		#
#rsem-calculate-expression -p 4 --calc-pme --paired-end --bowtie2 --no-bam-output ../A2/out1.fq ../A2/out2.fq ../../transcriptome_data/ref expressionA2		#
#rsem-calculate-expression -p 4 --calc-pme --paired-end --bowtie2 --no-bam-output ../A3/out1.fq ../A3/out2.fq ../../transcriptome_data/ref expressionA3		#
#rsem-calculate-expression -p 4 --calc-pme --paired-end --bowtie2 --no-bam-output ../A4/out1.fq ../A4/out2.fq ../../transcriptome_data/ref expressionA4		#
#rsem-calculate-expression -p 4 --calc-pme --paired-end --bowtie2 --no-bam-output ../A5/out1.fq ../A5/out2.fq ../../transcriptome_data/ref expressionA5		#
# Process output																		#
R CMD BATCH ../../rsemNames.R																	#
R CMD BATCH ../../withinGeneRSEM.R																#
cd ..																				#
#################################################################################################################################################################


# STEP 4.f Sailfish
#################################################################################################################
# Run Sailfish													#
mkdir sailfish													#
cd sailfish													#
sailfish quant -i ../../transcriptome_data/indDir -p 4 -l IU -1 ../A1/out1.fq -2 ../A1/out2.fq -o A1	#
sailfish quant -i ../../transcriptome_data/indDir -p 4 -l IU -1 ../A2/out1.fq -2 ../A2/out2.fq -o A2	#
sailfish quant -i ../../transcriptome_data/indDir -p 4 -l IU -1 ../A3/out1.fq -2 ../A3/out2.fq -o A3	#
sailfish quant -i ../../transcriptome_data/indDir -p 4 -l IU -1 ../A4/out1.fq -2 ../A4/out2.fq -o A4	#
sailfish quant -i ../../transcriptome_data/indDir -p 4 -l IU -1 ../A5/out1.fq -2 ../A5/out2.fq -o A5	#
# Process output												#
R CMD BATCH ../../sailfishNames.R										#
R CMD BATCH ../../withinGeneSailfish.R										#
cd ..														#
#################################################################################################################


# STEP 4.g Tigar2
# IMPORTANT: The user should manually provide the full path to Tigar2 source code.
# IMPORTANT: 12 cores should be available here
#################################################################################################################################################
# Run Tigar2																	#
mkdir tigar																	#
cd tigar																	#
# NOTE: Each of the following commands can also run in parallel for each replicate (highly recommended).					#
java -jar YOUR_FULL_PATH_TO_TIGAR2/Tigar2.jar ../../transcriptome_data/output.fasta ../bowtie/A1.sam --is_paired --alpha_zero 0.1 A1.txt	#
java -jar YOUR_FULL_PATH_TO_TIGAR2/Tigar2.jar ../../transcriptome_data/output.fasta ../bowtie/A2.sam --is_paired --alpha_zero 0.1 A2.txt	#
java -jar YOUR_FULL_PATH_TO_TIGAR2/Tigar2.jar ../../transcriptome_data/output.fasta ../bowtie/A3.sam --is_paired --alpha_zero 0.1 A3.txt	#
java -jar YOUR_FULL_PATH_TO_TIGAR2/Tigar2.jar ../../transcriptome_data/output.fasta ../bowtie/A4.sam --is_paired --alpha_zero 0.1 A4.txt	#
java -jar YOUR_FULL_PATH_TO_TIGAR2/Tigar2.jar ../../transcriptome_data/output.fasta ../bowtie/A5.sam --is_paired --alpha_zero 0.1 A5.txt	#
# Process output																#
R CMD BATCH ../../tigarNames.R															#
R CMD BATCH ../../withinGeneTigar.R														#
cd ..																		#
#################################################################################################################################################


# STEP 4.h eXpress
#########################################################################################
# Run eXpress										#
mkdir xpress										#
cd xpress										#
express -o A1 --no-update-check ../../transcriptome_data/known.fa ../bowtie/A1.sam 	#
express -o A2 --no-update-check ../../transcriptome_data/known.fa ../bowtie/A2.sam 	#
express -o A3 --no-update-check ../../transcriptome_data/known.fa ../bowtie/A3.sam 	#
express -o A4 --no-update-check ../../transcriptome_data/known.fa ../bowtie/A4.sam 	#
express -o A5 --no-update-check ../../transcriptome_data/known.fa ../bowtie/A5.sam 	#
# Process output									#
R CMD BATCH ../../xpressNames.R								#
R CMD BATCH ../../withinGeneXpress.R							#
cd ..											#
#########################################################################################



# STEP 5. Produce graphs
#################################################################################################
# Summarize results										#
R CMD BATCH ../getResults.R									#
R CMD BATCH ../withinGenePlots.R  #this produces the scatterplots of within Gene Estimates	#
R CMD BATCH ../criteria.R	  #this produces the graph with the ranking of methods		#
#################################################################################################



