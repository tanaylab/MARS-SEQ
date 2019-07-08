#!/bin/sh
#$ -S /bin/sh

if [ "$#" -ne 2 ]
then
	echo usage : run_mapping.sh scdb_path njobs
	exit
fi

scdb_path=`readlink -f $1`

FASTQ_LIST=$scdb_path/_temp/fastq_list.txt
AMP_BATCHES_TO_PROCESS=$scdb_path/config/amp_batches_to_process.txt
R_HOME=`grep R_HOME $scdb_path/config/config.txt | cut -f2 -d"="`
scRNA_scripts=`grep scRNA_scripts $scdb_path/config/config.txt | cut -f2 -d"="`

mkdir -p $scdb_path/output/
mkdir -p $scdb_path/output/QC/
mkdir -p $scdb_path/_labeled_raw_reads/
mkdir -p $scdb_path/_trimmed_mapped_reads/
mkdir -p $scdb_path/_temp
mkdir -p $scdb_path/_logs
\rm -r $scdb_path/_logs/*

annot_problems_flag=$($R_HOME/Rscript $scRNA_scripts/check_annots.r $scdb_path/)
if [ "$annot_problems_flag" != 0 ]
then
	echo ANNOTATION ERROR!!!
	cat $scdb_path/_logs/check_annots.log
	exit
fi

fastq_problems_flag=$($R_HOME/Rscript $scRNA_scripts/create_fastq_list.r $scdb_path/ $AMP_BATCHES_TO_PROCESS $FASTQ_LIST)
if [ "$fastq_problems_flag" != 0 ]
then
	echo ANNOTATION ERROR!!!
	exit
fi


NJOBS=`wc -l $FASTQ_LIST | cut -f 1 -d ' '`
NJOBS=`echo $NJOBS-1 | bc`
NTASKS=$NJOBS;

MAX_JOBS=$2
echo Mapping $NTASKS fastq files
touch qsub_log
i=1
while [ "$i" -le $NTASKS ]
do
  running_jobs=`qstat | tail -n +3 | wc -l `
  while [ $running_jobs -ge $MAX_JOBS ]
	do
	sleep 10
	running_jobs=`qstat | tail -n +3 | wc -l `
  done
  to=`echo $i+$MAX_JOBS-$running_jobs-1 | bc`

  if [ "$to" -gt "$NTASKS" ];then
	to=$NTASKS
  fi
  echo submitting $i-$to
  qsub -q all.q@@dell6420-384g -pe threads 20 -wd $scdb_path -t $i-$to $scRNA_scripts/distrib_mapping.sh >>qsub_log


  i=`echo $to+1| bc`
  sleep 30

done


running_jobs=`qstat | tail -n +3 | wc -l `
while [ $running_jobs -gt 0 ]
do
  sleep 10
  running_jobs=`qstat | tail -n +3 | wc -l `
done


$R_HOME/Rscript $scRNA_scripts/check_status.r $scdb_path/ mapping $NTASKS

$R_HOME/Rscript $scRNA_scripts/merge_seq_batches_tables.r $scdb_path/
