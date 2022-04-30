argv <- commandArgs(TRUE)
sleep_time=sample(1:180,1)
sleep_time=1
message("Sleeping for ",sleep_time," Seconds")
Sys.sleep(sleep_time)

stop1=function(msg){
  cmd_rm=paste("\rm -f /tmp/scRNA_",seq_batch_id,"_",prefix,".*",sep="")
  system(cmd_rm)
  stop(msg)
}

task_id=argv[1]
scripts_path=argv[2]
status_fn=argv[3]
#task_id=1
config=as.data.frame(t(read.table("config/config.txt",sep="=",header=F,stringsAsFactors=F,row.names=1)), stringsAsFactors=F)
fastq_list=read.table("_temp/fastq_list.txt",header=T,stringsAsFactors=F)
fastq_details=fastq_list[task_id,]
seq_batches=read.delim("annotations/seq_batches.txt",header=T,stringsAsFactors=F)
seq_batch_id=fastq_details$Seq_batch_ID
seq_batch_details=seq_batches[match(seq_batch_id,seq_batches$Seq_batch_ID),]
seq_batch_details[is.na(seq_batch_details)]=""
R1_file_path=paste("raw_reads",seq_batch_id,fastq_details$R1_file,sep="/")
R2_file_path=paste("raw_reads",seq_batch_id,fastq_details$R2_file,sep="/")
prefix=gsub(".fastq.gz","",fastq_details$R1_file)
fn=paste(prefix,".fastq.gz",sep="")
labeled_fastq_path=paste("_labeled_raw_reads",seq_batch_id,fn,sep="/")
extract_labels_stats_path=paste("_temp/QC",seq_batch_id,paste(prefix,".txt",sep=""),sep="/")
mapped_sam_path=paste("/tmp/scRNA_",seq_batch_id,"_",prefix,".sam",sep="")
trimmed_mapped_sam_path=paste("_trimmed_mapped_reads/",seq_batch_id,"/",prefix,".sam",sep="")

if (!file.exists(paste("_labeled_raw_reads/",seq_batch_id,sep=""))){
  dir.create(paste("_labeled_raw_reads/",seq_batch_id,sep=""))
}

if (!file.exists("_temp/QC/")){
  dir.create("_temp/QC/")
}
if (!file.exists(paste("_temp/QC/",seq_batch_id,sep=""))){
  dir.create(paste("_temp/QC/",seq_batch_id,sep=""))
}

gzip_r1_flag=F
#if (substr(R1_file_path,nchar(R1_file_path)-2,nchar(R1_file_path))==".gz"){
#  new_R1_file_path=paste("/tmp/scRNA_",seq_batch_id,"_",substr(fastq_details$R1_file,1,nchar(fastq_details$R1_file)-3),sep="")
#  cmd01=paste("gzip -dc",R1_file_path,">",new_R1_file_path)
#  message(cmd01)
#  system(cmd01)
#  R1_file_path=new_R1_file_path
#  gzip_r1_flag=T
#}

gzip_r2_flag=F
#if (substr(R2_file_path,nchar(R2_file_path)-2,nchar(R2_file_path))==".gz"){
#  new_R2_file_path=paste("/tmp/scRNA_",seq_batch_id,"_",substr(fastq_details$R2_file,1,nchar(fastq_details$R2_file)-3),sep="")
#  cmd02=paste("gzip -dc",R2_file_path,">",new_R2_file_path)
#  message(cmd02)
#  system(cmd02)
#  R2_file_path=new_R2_file_path
#  gzip_r2_flag=T
#}



cmd1=paste("perl",paste(scripts_path,"/extract_labels.pl",sep=""),R1_file_path,R2_file_path,seq_batch_details$R1_design,seq_batch_id,config$oligos_file,labeled_fastq_path,extract_labels_stats_path,".",seq_batch_details$I5_design)
message(cmd1)
ok=system(cmd1)
if (ok!=0){
  stop1("ERROR: extract_labels.pl crashed!")
}
additional_bowtie_params=""
if (!is.null(config$bowtie_params)){
  additional_bowtie_params=paste(" ",config$bowtie_params,sep="")
}

cmd2=paste(config$bowtie_bin, "-x",config$bowtie_index,"-U",labeled_fastq_path,additional_bowtie_params,sep=" ")

# use star mapper if defined in config and overwrite bowtie
if (!is.null(config$star_bin))
{
  additional_star_params=""
  if (!is.null(config$star_params))
  {
  additional_star_params=paste(" ",config$star_params,sep="")
  }
  cmd2 = paste(config$star_bin,
               "--genomeDir", config$star_index,
               "--readFilesIn", labeled_fastq_path,
               additional_star_params,
               "--outFileNamePrefix", paste0("_logs/star.",prefix,"."),
               "--outStd SAM --outFilterMultimapScoreRange 2 --readFilesCommand zcat",
               sep=" ")
  cmd2 = paste('rm -f', mapped_sam_path, ';', cmd2)
}

if (!is.null(config$maternal_bcf) || !is.null(config$paternal_bcf)) {
    bcfs = character()
    if (!is.null(config$maternal_bcf)) {
        bcfs = c(bcfs, config$maternal_bcf)
    }
    if (!is.null(config$paternal_bcf)) {
        bcfs = c(bcfs, config$paternal_bcf)
    }
    hyppy = paste0("--bcf", 1:length(bcfs), " ",  bcfs, collapse=" ")
    hyppy = paste0("~/src/tlsrc/tools/hyppy/hyp.py ", hyppy)
    cmd2 = paste0(cmd2, " | ", hyppy)
}
cmd2 = paste0("/bin/bash -c 'source ~/.bashrc;", cmd2, " > ", mapped_sam_path, "'")
message(cmd2)
ok=system(cmd2)
if (ok!=0){
  stop1("ERROR: STAR / bowtie crashed!")
}

#cmd2.5=paste("gzip -f",labeled_fastq_path)
#message(cmd2.5)
#ok=system(cmd2.5)
#if (ok!=0){
#  stop1("ERROR!")
#}

flag_counts=table(read.table(pipe(paste("grep -v ^@ ",mapped_sam_path," | cut -f2",sep="")))[,1])

stats=read.table(extract_labels_stats_path,header=T)
stats$mapped=flag_counts["0"]+flag_counts["16"]

write.table(file=extract_labels_stats_path,stats,sep="\t",col.names=T,row.names=F,quote=F)


if (!file.exists(paste("_trimmed_mapped_reads/",seq_batch_id,sep=""))){
  dir.create(paste("_trimmed_mapped_reads/",seq_batch_id,sep=""))
}
cmd3=paste("cut -f1-9,12- ",mapped_sam_path," > " ,trimmed_mapped_sam_path,sep="")
message(cmd3)
ok=system(cmd3)
if (ok!=0){
  stop1("ERROR!")
}


#if (gzip_r1_flag){
#  cmd4=paste("rm -f",R1_file_path)
#  message(cmd4)
#  ok=system(cmd4)
#} else{
#  cmd4.1=paste("gzip",R1_file_path)
#   message(cmd4.1)
#  ok=system(cmd4.1)
#}
#if (ok!=0){
#  stop1("ERROR!")
#}
#
#if (gzip_r2_flag){
#  cmd5=paste("rm -f",R2_file_path)
#  message(cmd5)
#  ok=system(cmd5)
#} else{
#  cmd5.1=paste("gzip",R2_file_path)
#  message(cmd5.1)
#  ok=system(cmd5.1)
#}
#if (ok!=0){
#  stop1("ERROR!")
#}

cmd6=paste("rm -f",mapped_sam_path)
message(cmd6)
ok=system(cmd6)
if (ok!=0){
  stop1("ERROR!")
}
write.table(file=status_fn,"OK",row.names=F,col.names=F,quote=F)
