#!/bin/bash

ACC_FILE=$(find /home -name "rna_sample_acc.txt" | grep "home")
echo "ACC file is $ACC_FILE"
ACC_PATH=$(dirname $ACC_FILE)
GTF_DIR="/home/tom/Documents/Spring2021/Honors_Project/Tools/STAR/human_GRCh38_13/GCF_000001405.39_GRCh38.p13_genomic.gtf" #change this so it finds the genome directory
FC_OPTS=" -T 22 -p -a $GTF_DIR "
cd $ACC_PATH
echo "Moving working directory to $(pwd)"

LIST_PRESENT=$(find . -name "rna_acc_completed_list.txt"| grep -c "")
if (( $LIST_PRESENT != 0 ))
then
	echo "list of completed accessions found"
else
	echo "list of completed accessions not found,"
	echo "starting from the beginning...."
	touch rna_acc_completed_list.txt
fi
COMPLETED_LIST="rna_acc_completed_list.txt"
ACCESSIONS_NUM=$(wc -l $ACC_FILE | cut -f1 -d ' ')
ACCESSIONS_DONE=$(wc -l $COMPLETED_LIST | cut -f1 -d ' ')
ACCESSIONS_LEFT=$(($ACCESSIONS_NUM - $ACCESSIONS_DONE))
echo "$ACCESSIONS_LEFT accessions left to process in $ACC_FILE"
while(($ACCESSIONS_LEFT != 0))
do
	i=$(($ACCESSIONS_NUM - $ACCESSIONS_LEFT))
	((i++))
	CURRENT_ACC=$(head -n $i $ACC_FILE | tail -n 1)
	echo "Downloading $CURRENT_ACC reads"
	READ_1=$(grep $CURRENT_ACC List_of_accessions_and_links.txt | cut -f 7 | cut -d ";" -f 1)
	READ_2=$(grep $CURRENT_ACC List_of_accessions_and_links.txt | cut -f 7 | cut -d ";" -f 2)
	READ_1_NAME=$(echo $READ_1 | cut -d "/" -f 7)
	READ_2_NAME=$(echo $READ_2 | cut -d "/" -f 7)
	curl --max-time 300 $READ_1 -o $READ_1_NAME
	curl --max-time 300 $READ_2 -o $READ_2_NAME
	NUM_DOWNLOADED=$(($(ls | grep -c $READ_1_NAME) + $(ls | grep -c $READ_2_NAME)))
	if(($NUM_DOWNLOADED != 2))
	then
		COUNT_DOWN=5
		while(($NUM_DOWNLOADED != 2))
		do
			if(($COUNT_DOWN == 0))
			then
				echo "Could not download acc. from $READ_1 or $READ_2"
				exit 1
			fi
			curl --max-time 150 $READ_1 -o $READ_1_NAME
			curl --max-time 150 $READ_2 -o $READ_2_NAME
			NUM_DOWNLOADED=$(($(ls | grep -c $READ_1_NAME) + $(ls | grep -c $READ_2_NAME)))
			((COUNT_DOWN--))
		done
	fi
	trim_galore --cores 8 --paired --basename "trimmed_read" $READ_1_NAME $READ_2_NAME
	rm $READ_1_NAME $READ_2_NAME
	bowtie2 -q --threads 22 -x ref_gen -1 trimmed_read_val_1.fq.gz -2 trimmed_read_val_2.fq.gz -S temp.sam
	featureCounts $FC_OPTS -o temp_out.txt temp.sam
	cut -f1,7,8,9,10,11,12 temp_out.txt > "${CURRENT_ACC}_matrix.txt"
	rm trimmed_read_val_1.fq.gz trimmed_read_val_2.fq.gz temp.sam temp_out.txt.summary temp_out.txt *trimming_report.txt 
	
	OUTPUT_SIZE=$(stat --format=%s "${CURRENT_ACC}_matrix.txt")
	if(($OUTPUT_SIZE != 0))
	then
		echo $CURRENT_ACC >> rna_acc_completed_list.txt
		((ACCESSIONS_LEFT--))
		echo "$ACCESSIONS_LEFT left"
	fi
done
