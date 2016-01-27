#!/bin/bash
#$ -N scdb_demult
#$ -S /bin/bash
#$ -e _logs/demultiplex_err.$TASK_ID
#$ -o _logs/demultiplex_out.$TASK_ID

AMP_BATCH_ID=`head -n $SGE_TASK_ID config/amp_batches_to_process.txt | tail -n 1`
R_HOME=`grep R_HOME config/config.txt | cut -f2 -d"="`
DATA_DIR=`grep data_dir config/qc_config.txt | cut -f2 -d"="`
scRNA_scripts=`grep scRNA_scripts config/config.txt | cut -f2 -d"="`

echo demultiplexing amplification batch $AMP_BATCH_ID
mkdir -p _debug/$AMP_BATCH_ID
echo perl $scRNA_scripts/demultiplex.pl . $AMP_BATCH_ID

perl $scRNA_scripts/demultiplex.pl . $AMP_BATCH_ID
echo $R_HOME/Rscript $scRNA_scripts/make_batch_qc.r $AMP_BATCH_ID $DATA_DIR
$R_HOME/Rscript $scRNA_scripts/make_batch_qc.r $AMP_BATCH_ID $DATA_DIR

status_file=_logs/demultiplex_status.$SGE_TASK_ID
touch $status_file
echo OK > $status_file
