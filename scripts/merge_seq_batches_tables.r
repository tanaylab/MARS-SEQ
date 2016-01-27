args=commandArgs(trailingOnly = TRUE)
scdb_path=args[1]

fastqs=read.table(paste(scdb_path,"/_temp/fastq_list.txt",sep=""),stringsAsFactors=F,header=T)
l=list()
for (i in 1:dim(fastqs)[1]){
 # print(i)
  a=t(read.table(paste(scdb_path,"/_temp/QC/",fastqs[i,3],"/",gsub("fastq","txt",gsub("fastq.gz","txt",fastqs[i,1])),sep=""),header=T)[1,])
  if (is.null(l[[fastqs[i,3]]])){
    l[[fastqs[i,3]]]=a
  }
  l[[fastqs[i,3]]]=l[[fastqs[i,3]]]+a
}
df=t(as.data.frame(l))
rownames(df)=names(l)
colnames(df)=rownames(l[[1]])
write.table(file=paste(scdb_path,"/output/QC/seq_batches_stats.txt",sep=""),df,quote=F,sep="\t",row.names=T,col.names=T)
