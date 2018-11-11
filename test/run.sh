../miccr.py -d ./Refseq_viruses.fa     -i Refseq_viruses.fa   -p TEST_1 -v
../miccr.py -d ./Refseq_viruses.fa.mmi -i Refseq_viruses.fa   -p TEST_2 -v
../miccr.py -d ./Refseq_viruses.fa     -i Refseq_viruses.fa   -p TEST_3 -v  -x asm5    
../miccr.py -d ./Refseq_viruses.fa     -i ERR2538127.fastq.gz -c            -x map-ont 
../miccr.py -d ./Refseq_viruses.fa     -i ERR2538127.fastq.gz -p TEST_5 -v  -x map-ont 
../miccr.py -d ./Refseq_viruses.fa     -i ERR2538127.fastq.gz -p TEST_6 -v  -x map-ont -mp 0.05
../miccr.py -d ./Refseq_viruses.fa     -i ERR2538127.fastq.gz -p TEST_7 -v  -x map-ont -if 0.5
