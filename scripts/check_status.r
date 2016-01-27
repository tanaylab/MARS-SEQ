args=commandArgs(trailingOnly = TRUE)
scdb_path=args[1]
component=args[2]
ntasks=args[3]

error_status=file.exists(paste(scdb_path,"/_logs/",component,"_status.",1:ntasks,sep=""))
if (sum(!error_status)==0){
  status="OK"
  message(component, ": finished successfully")
} else {
  message(component, ": ERROR! ",sum(!error_status),"/",ntasks," tasks failed")
  status=paste("ERROR: ",sum(!error_status),"/",ntasks," (err/total)", sep="") 
 }

write.table(file=paste(scdb_path,"/_logs/",component,"_status",sep=""),status,quote=F,row.names=F,col.names=F)
