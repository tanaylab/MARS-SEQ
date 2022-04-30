#!/bin/sh
#$ -N bt
#$ -S /bin/sh

if [ "$#" -eq 3 ]
then
    suffix=$3
elif [ "$#" -eq 2 ]
then
    suffix=""
else
    echo usage : run_demultiplexing.sh scdb_path njobs suffix
    exit
fi

if [ "$SGE_ROOT" = "" ]
then
    QCMD="bjobs"
else
    suffix=""
    QCMD="qstat"
fi

scdb_path=`readlink -f $1`
AMP_BATCHES_TO_PROCESS=$scdb_path/config/amp_batches_to_process.txt
scRNA_scripts=`grep scRNA_scripts $scdb_path/config/config.txt | cut -f2 -d"="`
R_HOME=`grep R_HOME $scdb_path/config/config.txt | cut -f2 -d"="`

mkdir -p $scdb_path/_debug/
mkdir -p $scdb_path/output/
mkdir -p $scdb_path/output/QC/
mkdir -p $scdb_path/output/QC/rd/
mkdir -p $scdb_path/_temp
mkdir -p $scdb_path/_logs
rm -f $scdb_path/_logs/demultiplex*


N_AMP_BATCHES=`wc -l $AMP_BATCHES_TO_PROCESS | cut -f 1 -d ' '`
echo "Demultiplexing $N_AMP_BATCHES amplificaiton batches" ;


MAX_JOBS=$2
i=1
cd $scdb_path
while [ "$i" -le $N_AMP_BATCHES ]
do
  running_jobs=`$QCMD | grep -c sc_dm$suffix `
  while [ $running_jobs -ge $MAX_JOBS ]
    do
    sleep 10
  running_jobs=`$QCMD | grep -c sc_dm$suffix `
  done
  to=`echo $i+$MAX_JOBS-$running_jobs-1 | bc`
 
  if [ "$to" -gt "$N_AMP_BATCHES" ];then
    to=$N_AMP_BATCHES
  fi
  echo submitting $i-$to

  if [ "$QCMD" = "qstat" ]
  then
    qsub -wd $scdb_path -t $i-$to  $scRNA_scripts/distrib_demultiplex.sh >>qsub_log
  else
     if [ "$LSB_QUEUE" = "" ]
     then
        $scRNA_scripts/cmds_runner.py \
        "bsub -oo _logs/demultiplex_{c1}.out -eo _logs/demultiplex_{c1}.err -J {c1}_sc_dm$suffix -R rusage[mem=32000] $scRNA_scripts/distrib_demultiplex.sh {c1} " $i:$to >> qsub_log
     else
        $scRNA_scripts/cmds_runner.py \
        "bsub -q $LSB_QUEUE -oo _logs/demultiplex_{c1}.out -eo _logs/demultiplex_{c1}.err -J {c1}_sc_dm$suffix -R rusage[mem=32000] $scRNA_scripts/distrib_demultiplex.sh {c1} " $i:$to >> qsub_log
     fi
  fi
  i=`echo $to+1| bc`
  sleep 30  

done


running_jobs=`$QCMD | grep -c sc_dm$suffix `

while [ $running_jobs -gt 0 ]
do
  sleep 10
  running_jobs=`$QCMD | grep -c sc_dm$suffix `
done


$R_HOME/Rscript $scRNA_scripts/check_status.r $scdb_path demultiplex $N_AMP_BATCHES
