
library(gplots)
if (! exists("output_dir")){
  args=commandArgs(trailingOnly = TRUE)
  if (length(args)==2){
    scdb_path=args[1]
    output_dir=args[2]
  }else{
    stop("usage: Rscript make_qc_report.r scdb_path data_dir")
    
  }
}
get_stats_per_seq_batch=function(seq_batch){
  n_umis_per_cell=colSums(umitab)
  cell_mask=intersect(rownames(sample_list)[sample_list$sequencing_batch%in%seq_batch],single_cell_mask)
  n_cells=length(cell_mask)
  avg_mapping_stats=colMeans(t(mapping_stats)[cell_mask,],na.rm=T)
  avg_n_RMT=mean( n_umis_per_cell[cell_mask])
  total_n_RMT=sum(n_umis_per_cell[cell_mask])
  names(avg_mapping_stats)=c("avg_n_mouse_reads","avg_n_ecoli_reads","avg_n_ercc_reads","avg_n_raw_reads")
  total_n_raw_reads=as.numeric(avg_mapping_stats["avg_n_raw_reads"]*n_cells)
  avg_n_reads_per_RMT=total_n_raw_reads/total_n_RMT
  return(round(c(n_cells=n_cells,avg_mapping_stats,avg_n_RMTs=avg_n_RMT,total_n_RMTs=total_n_RMT,total_n_raw_reads=total_n_raw_reads,avg_n_reads_per_RMT=avg_n_reads_per_RMT)))
}

try(system(paste("mkdir -p ",scdb_path,"/",output_dir,"/QC_reports",sep=""),ignore.stderr=T))

sample_list=read.table(paste(scdb_path,"/annotations/wells_cells.txt",sep=""),stringsAsFactors=F,header=T,sep="\t")
amp_batches=read.delim(paste(scdb_path,"/annotations/amp_batches.txt",sep=""),stringsAsFactors=F,header=T,sep="\t")

rownames(sample_list)=sample_list$Well_ID

main_stats_list=list()

for (batch in amp_batches$Amp_batch_ID){
   try({
    load(paste(scdb_path,"/",output_dir,"/QC/rd/",batch,".rd",sep=""))
    main_stats_list[[batch]]=main_stats
  },silent=T)
}
m=t(sapply(main_stats_list,function(x){return(x[names(main_stats_list[[1]])])}))

pdf(paste(scdb_path,"/_temp/additional_qc.pdf",sep=""))
m2=t(sapply(main_stats_list,function(x){return(x[c("gene_umis","gene_umis_neg_control","avg_noffsets_per_rmt","avg_noffsets_per_rmt_neg_ctrl")])}))


# RMTs vs. QC1
try({
plot(amp_batches$QC1..pooled.cDNA.,m2[match(amp_batches[,1],rownames(m2)),1],xlab="QC1",ylab="average #RMTs",pch="",log="y")
text(amp_batches$QC1..pooled.cDNA.,m2[match(amp_batches[,1],rownames(m2)),1],labels=amp_batches[,1],cex=.5)
})

# RMTs vs. QC2
try({
plot(amp_batches$QC2..RT.2.,m2[match(amp_batches[,1],rownames(m2)),1],xlab="QC2",ylab="average #RMTs",pch="",log="y")
text(amp_batches$QC2..RT.2.,m2[match(amp_batches[,1],rownames(m2)),1],labels=amp_batches[,1],cex=.5)
})

# RMTs vs. QC3
try({
plot(amp_batches$QC3..library.Ct.,m2[match(amp_batches[,1],rownames(m2)),1],xlab="QC3",ylab="average #RMTs",pch="",log="y")
text(amp_batches$QC3..library.Ct.,m2[match(amp_batches[,1],rownames(m2)),1],labels=amp_batches[,1],cex=.5)
})
dev.off()

for (seq_batch in unique(amp_batches$Seq_batch_ID)){
  cur_amp_batches=amp_batches$Amp_batch_ID[amp_batches$Seq_batch_ID==seq_batch]
  pdfs=paste(scdb_path,"/",output_dir,"/QC/report_per_amp_batch/",cur_amp_batches,".pdf",sep="")
  pdfs=pdfs[file.exists(pdfs)]
  cmd=paste(paste("gs -dBATCH -dNOPAUSE -dAutoRotatePages=/None -q -sDEVICE=pdfwrite -sOutputFile=",scdb_path,"/",output_dir,"/QC_reports/qc_",seq_batch,".pdf",sep=""),paste(pdfs,collapse=" "))
  
  system(cmd)
}

df=data.frame(amp_batch=rownames(m),m)
write.table(file=paste(scdb_path,"/",output_dir,"/amp_batches_stats.txt",sep=""),df,row.names=F,col.names=T,sep="\t",quote=F)

rownames(amp_batches)=amp_batches[,"Amp_batch_ID"]
m4=cbind(amp_batch_ID=rownames(m),amp_batches[rownames(m),c("Experiment_ID","Description","QC2..RT.2.")],m[,match(c("reads","spike_yield","gene_umis","gene_umis_neg_control","noise_estimation","avg_noffsets_per_rmt","avg_reads_per_rmt"),colnames(m))])
 write.table(file=paste(scdb_path,"/",output_dir,"/amp_batches_summary.txt",sep=""),m4,row.names=F,col.names=T,sep="\t",quote=F)


