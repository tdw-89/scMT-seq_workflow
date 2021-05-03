#!/bin/bash

ACC_FILE=$(find /home -name "met_sample_acc.txt" | grep "home")
echo "ACC file is $ACC_FILE"
ACC_PATH=$(dirname $ACC_FILE)
GENOME_DIR="/home/tom/Documents/Spring2021/ref_genome/human" #change this so it finds the genome directory
BIS_OPTS="--parallel 2 --non_directional --genome $GENOME_DIR"
cd $ACC_PATH
echo "Moving working directory to $(pwd)"

LIST_PRESENT=$(find . -name "met_acc_completed_list.txt"| grep -c "")
if (( $LIST_PRESENT != 0 ))
then
	echo "list of completed accessions found"
else
	echo "list of completed accessions not found,"
	echo "starting from the beginning...."
	touch met_acc_completed_list.txt
fi
COMPLETED_LIST="met_acc_completed_list.txt"
ACCESSIONS_NUM=$(wc -l $ACC_FILE | cut -f1 -d ' ')
ACCESSIONS_DONE=$(wc -l $COMPLETED_LIST | cut -f1 -d ' ')
ACCESSIONS_LEFT=$(($ACCESSIONS_NUM - $ACCESSIONS_DONE))
echo "$ACCESSIONS_LEFT accessions left to process in $ACC_FILE"
while(($ACCESSIONS_LEFT != 0))
do
	i=$(($ACCESSIONS_NUM - $ACCESSIONS_LEFT))
	((i++))
	CURRENT_ACC=$(head -n $i $ACC_FILE | tail -n 1)
	READ_MET=$(grep $CURRENT_ACC List_of_accessions_and_links.txt | cut -f 7 | cut -d ";" -f 1)
	READ_MET_NAME=$(echo $READ_MET | cut -d "/" -f 7)
	
	echo "Downloading $READ_MET_NAME"
	NUM_DOWNLOADED=0
	DOWN_LENGTH=0
	TEST_LENGTH=$(curl -sI $READ_MET | grep -i Content-Length | cut -d " " -f 2)
	TEST_LENGTH=$(printf %d $TEST_LENGTH)
	echo "Test length: $TEST_LENGTH"
	curl $READ_MET -o $READ_MET_NAME #REPLACE WITH: curl --max-time 150 $READ_MET -o $READ_MET_NAME
	DOWN_LENGTH=$(wc -c $READ_MET_NAME | cut -d " " -f 1)
	DOWN_LENGTH=$(printf %d $DOWN_LENGTH)
	echo "Download length: $DOWN_LENGTH"
	if [ $DOWN_LENGTH -eq $TEST_LENGTH ]
	then
		echo "Download Complete"
		NUM_DOWNLOADED=1
	fi
	
	if(($NUM_DOWNLOADED != 1))
	then
		COUNT_DOWN=5
		while(($NUM_DOWNLOADED != 1))
		do
			if(($COUNT_DOWN == 0))
			then
				echo "Could not download acc. from $READ_MET"
				exit 1
			fi
			curl --max-time 150 $READ_MET -o $READ_MET_NAME
			if(($(ls | grep -c $READ_MET_NAME) == 1))
			then
				DOWN_LENGTH=$(wc -c $READ_MET_NAME | cut -d " " -f 1)
				DOWN_LENGTH=$(printf %d $DOWN_LENGTH)
			fi
			if [ $DOWN_LENGTH -eq $TEST_LENGTH ]
			then
				NUM_DOWNLOADED=1
			fi
			((COUNT_DOWN--))
		done
	fi
	
	trim_galore --cores 8 --clip_R1 6 --basename "temp_met" $READ_MET_NAME
	bismark $BIS_OPTS temp_met_trimmed.fq.gz
	rm temp_met_trimmed.fq.gz $READ_MET_NAME *trimming_report.txt temp_met_trimmed_bismark_bt2_SE_report.txt
	deduplicate_bismark --single "temp_met_trimmed_bismark_bt2.bam"
	rm temp_met_trimmed_bismark_bt2.bam
	mkdir "${READ_MET_NAME}_met"
	bismark_methylation_extractor --parallel 6 --single-end --comprehensive --no_header --output "./${READ_MET_NAME}_met" temp_met_trimmed_bismark_bt2.deduplicated.bam
	echo $CURRENT_ACC >> met_acc_completed_list.txt
	((ACCESSIONS_LEFT--))
	echo "$ACCESSIONS_LEFT left"
done
