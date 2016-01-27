args=commandArgs(trailingOnly = TRUE)
scdb_path=args[1]
amp_batches_fn=args[2]
output_fn<- args[3]


seq_batches=read.table(paste(scdb_path,"/annotations/seq_batches.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
amp_batches=read.delim(paste(scdb_path,"/annotations/amp_batches.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
do_one_batch=function(seq_batch_id){
  files_R1=list.files(path=paste(scdb_path,"/raw_reads/",seq_batch_id,sep=""),pattern="_R1_")
  files_R1=files_R1[substr(files_R1,nchar(files_R1)-7,nchar(files_R1))=="fastq.gz"|substr(files_R1,nchar(files_R1)-4,nchar(files_R1))=="fastq"]
  files_R2=list.files(path=paste(scdb_path,"/raw_reads/",seq_batch_id,sep=""),pattern="_R2_")
  files_R2=files_R2[substr(files_R2,nchar(files_R2)-7,nchar(files_R2))=="fastq.gz"|substr(files_R2,nchar(files_R2)-4,nchar(files_R2))=="fastq"]
  
  if (any(gsub("_R1_","_R2_",files_R1)!=files_R2)){
    cat(1)
    stop("There is a discrepancy betweeen R1 and R2 files of sequencing batch ",seq_batch_id)
  }
  if (length(files_R1)==0){
    cat(1)
    stop(paste(scdb_path,"/raw_reads/",seq_batch_id,sep="")," doesn't contain any fastq files")
  }
  return(rbind(files_R1=files_R1,files_R2=files_R2,Seq_batch_ID=seq_batch_id))
}


amp_batches_to_process=read.delim(amp_batches_fn,stringsAsFactors=F,header=F)[,1]
#print (wells_cells$Amp_batch_ID)

seq_batches_to_process=unique(amp_batches$Seq_batch_ID[amp_batches$Amp_batch_ID%in%amp_batches_to_process])

m=t(as.data.frame(lapply(seq_batches_to_process,do_one_batch)))
colnames(m)=c("R1_file","R2_file","Seq_batch_ID")
write.table(m,file=output_fn,col.names=T,sep="\t",quote=F,row.names=F)

cat("0")
