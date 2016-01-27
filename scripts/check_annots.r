args=commandArgs(trailingOnly = TRUE)
scdb_path=args[1]


wells=read.table(paste(scdb_path,"/annotations/wells_cells.txt",sep=""),stringsAsFactors=F,header=T,sep="\t")
amp_batches=read.delim(paste(scdb_path,"/annotations/amp_batches.txt",sep=""),stringsAsFactors=F,header=T,sep="\t")
seq_batches=read.delim(paste(scdb_path,"/annotations/seq_batches.txt",sep=""),stringsAsFactors=F,header=T,sep="\t")
flag=F
msgl=list()
if (length(unique(wells[,1]))!=dim(wells)[1]){msgl[[length(msgl)+1]]="Well_ID is non-unique";flag=T}
if (length(unique(amp_batches[,1]))!=dim(amp_batches)[1]){msgl[[length(msgl)+1]]=("Amp_batch_ID is non-unique");flag=T}
if (length(unique(seq_batches[,1]))!=dim(seq_batches)[1]){msgl[[length(msgl)+1]]=("Seq_batch_ID is non-unique");flag=T}


uniq_well_barcodes=sapply(split(wells$Cell_barcode,wells$Amp_batch_ID),function(x){length(unique(x))==length(x)})[amp_batches[,1]]
if (sum(!uniq_well_barcodes)){msgl[[length(msgl)+1]]=paste("Well barcodes are non-unique in batches",paste(names(which(!uniq_well_barcodes)),collapse=", "));flag=T}

if (sum(!wells$Amp_batch_ID%in%amp_batches$Amp_batch_ID)>0){
  msgl[[length(msgl)+1]]=("Some wells are linked to undefined Amp_batch_ID")
  flag=T
}

if (sum(!amp_batches$Seq_batch_ID%in%seq_batches$Seq_batch_ID)>0){
    msgl[[length(msgl)+1]]=("Some amplification batches are linked to undefined Seq_batch_ID")
  flag=T
}
if (!flag){
  msgl[[length(msgl)+1]]="Annotation files are ok"
}

write.table(file=paste(scdb_path,"/_logs/check_annots.log",sep=""),unlist(msgl),quote=F,row.names=F,col.names=F,sep="\t")
cat(as.numeric(flag))
#if (flag){
#  sapply(msgl,message)
#  stop("Correct your annotation tables!")
#}else {
#  message("Annotation tables are OK")
#}
