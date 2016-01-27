#!/bin/bash
#$ -N scdb_map
#$ -S /bin/bash
#$ -e _logs/mapping_err.$TASK_ID
#$ -o _logs/mapping_out.$TASK_ID

echo TASK_ID=$SGE_TASK_ID 
echo HOSTNAME=$HOSTNAME
R_HOME=`grep R_HOME config/config.txt | cut -f2 -d"="`
scRNA_scripts=`grep scRNA_scripts config/config.txt | cut -f2 -d"="`

$R_HOME/Rscript $scRNA_scripts/distrib_mapping.r $SGE_TASK_ID $scRNA_scripts _logs/mapping_status.$SGE_TASK_ID

