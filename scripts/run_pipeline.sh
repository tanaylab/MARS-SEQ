#!/bin/sh
#$ -S /bin/sh

if [ "$#" -ne 1 ]
then
	echo usage : run_pipeline.sh scdb_path
	exit
fi

scdb_path=$1

rm -f $scdb_path/_logs/mapping_status
rm -f $scdb_path/_logs/demultiplex_status

R_HOME=`grep R_HOME $scdb_path/config/config.txt | cut -f2 -d"="`
DATA_DIR=`grep data_dir $scdb_path/config/qc_config.txt | cut -f2 -d"="`
scRNA_scripts=`grep scRNA_scripts $scdb_path/config/config.txt | cut -f2 -d"="`

dos2unix -q annotations/*.txt
dos2unix -q config/amp_batces_to_process.txt

echo "--------------------------------------------------"
echo "Running mapping"
echo
$scRNA_scripts/run_mapping.sh $scdb_path 30
mapping_status=`cat $scdb_path/_logs/mapping_status`
if [[ $mapping_status != "OK" ]]
then
	exit
fi
echo
echo "--------------------------------------------------"
echo "Running demultiplex"
echo
$scRNA_scripts/run_demultiplexing.sh $scdb_path 30
demultiplex_status=`cat $scdb_path/_logs/demultiplex_status`
if [[ $demultiplex_status != "OK" ]]
then
	exit
fi
echo
echo "--------------------------------------------------"
echo "Merging QC reports"
$R_HOME/Rscript $scRNA_scripts/make_qc_report.r $scdb_path $DATA_DIR > $scdb_path/_logs/make_qc_report.log  2> $scdb_path/_logs/make_qc_report.log

echo
echo "--------------------------------------------------"
echo "Done"

