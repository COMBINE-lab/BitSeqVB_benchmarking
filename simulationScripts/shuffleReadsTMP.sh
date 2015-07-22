#!/bin/bash
awk '{OFS="\t"; getline seq; \
                getline sep; \
                getline qual; \
                print $0,seq,sep,qual}' sim_1.fastq > sim_1.fastq.temp
awk '{OFS="\t"; getline seq; \
                getline sep; \
                getline qual; \
                print $0,seq,sep,qual}' sim_2.fastq > sim_2.fastq.temp
paste sim_1.fastq.temp sim_2.fastq.temp > pasteFile.txt
./../../../../../../mr01-home01/mqbssppe/coreutils-8.23/src/sort --parallel=4 -R pasteFile.txt > randomMerged.txt
awk '{OFS="\n"; print $1,$2,$3,$4}' randomMerged.txt > out1.fq
awk '{OFS="\n"; print $5,$6,$7,$8}' randomMerged.txt > out2.fq
rm sim_1.fastq.temp sim_2.fastq.temp randomMerged.txt pasteFile.txt

