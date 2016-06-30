#' Loads and merges data sets of umis. 
#' Merging occurs only if loading more than one data set.
#'
#' @param scr_file_path (character) Vector or one column file containing paths to umis data sets. 
#' The name of the umis table will be recognized as the corresponding cell batch.
#' @param scr_scdb_list (character) Vector or one column file containing the names of db annotations used
#' for every umis data set. That is, one data base per umis table. If no db annotation file is provided,
#' then this function will look for a 'db' pattern in the path and the corresponding directory will be 
#' an error will be prompted. Default is NULL.
#' @param info_as_files (logical) Logical value indicating if file or vector is provided. Files are expected by
#' default.
#' @param min_umis_n (Integer) Minimum number of umis to consider a cell into the umis data set. Default is 0,
#' then no cells are filtered out.
#' 
#' @return scram_obj with the following fields: 
#' \itemize{
#'  \item{"umis"}{UMI matrix (numeric matrix of umis with genes as rownames and cells as columns)}
#'  \item{"umis_path"}{}
#'  \item{"scr_batch"}{}
#'  \item{"scr_scdb"}{}
#'  \item{"unmatched_rownames_n"}{}
#'  \item{"filtered_cells_n"}{}
#' }
#'
#' @examples
#' 
#' scram.load_umis(scr_file_path, scr_scdb_list=NULL, Info_as_files=TRUE, min_umis_n=0)
#'
#' @export
scram.load_umis <- function(scr_file_path, scr_scdb_list=NULL, info_as_file=TRUE, min_umis_n=0, umitab_db="/home/eladch/proj/scRNA/analysis/e8.5_2015-06-25/umitab/umitab.txt"){

	rnames_miss_total <- 0
	
	if(info_as_file){
		files_list <- readLines(con=as.character(scr_file_path))
		if (!is.null(scr_scdb_list)){
			scdb_list <- readLines(con=as.character(scr_scdb_list))
		}else{
			scdb_list <- NULL
		}
	} else {
		files_list <- as.character(scr_file_path)
		if (!is.null(scr_scdb_list)){
			scdb_list <- as.character(scr_scdb_list)
		}else{
			scdb_list <- NULL
		}
	}
    
	if(!is.null(scdb_list) & length(files_list)!=length(scdb_list)){
		stop("Error: List of files and list of db annotations have different length. Each file must include the corresponding db annotation")
	}
	
	for (i in 1:length(files_list)){
		if (is.null(scdb_list) & length(grep("db",files_list[i]))<1){
			stop(paste("Error: No db annotation list was provided as well as no 'db' pattern was found in ",files_list[i],sep=""))
		}
		message("\tReading umis table number ",i)
		cell_names <- scan(umitab_db, what=character(), nlines=1)
		umis_tmp   <- data.table::fread(umitab_db, header=FALSE,skip=1)
		umis_tmp   <- as.data.frame(umis_tmp)
		if (length(cell_names)==length(colnames(umis_tmp))-1){
			colnames(umis_tmp) <- c("Row.names",cell_names)
		}else{
			colnames(umis_tmp)    <- cell_names
			colnames(umis_tmp)[1] <- "Row.names"
		}
		rownames(umis_tmp)    <- umis_tmp$Row.names
		if (length(unique(umis_tmp$row.names))!=length(umis_tmp$row.names)){
			stop(paste("Error: Duplicated row names in ",files_list[i],sep="")				
		}
		if (i==1){
			umis <- umis_tmp
			rm(umis_tmp)
			message("Umis table: ",files_list[i])
			message("dim(umis) ", paste(nrow(umis),ncol(umis)-1, collapse=" "))
			finfo_tmp <- strsplit(files_list[i],split="/")
			finfo_tmp <- unlist(finfo_tmp)
			scr_batch <- gsub(".txt","",finfo_tmp[length(finfo_tmp)])
			scr_scdb  <- ifelse(is.null(scdb_list),finfo_tmp[grep("db",finfo_tmp)[1]],scdb_list[i]) 
		}
		if (i!=1){
			umis <- merge(umis[,-1], umis_tmp[,-1], by.x="row.names", by.y=0, all=T)
			rnames_union   <- length(union(umis$row.names,umis_tmp$row.names))
			rnames_inter   <- length(intersect(umis$row.names,umis_tmp$row.names))
			rnames_missing <- rnames_union-rnames_inter
			rm(umis_tmp)
			message("After merging umis table: ",files_list[i])
			message("Number of unmatched row names = ",rnames_missing)
			message("New dim(umis) ", paste(nrow(umis),ncol(umis)-1, collapse=" "))
			rnames_miss_total <- rnames_miss_total+rnames_missing
			finfo_tmp     <- strsplit(files_list[i],split="/")
			finfo_tmp     <- unlist(finfo_tmp)
			scr_batch_tmp <- gsub(".txt","",finfo_tmp[length(finfo_tmp)]) 
			scr_scdb_tmp  <- ifelse(is.null(scdb_list),finfo_tmp[grep("db",finfo_tmp)[1]],scdb_list[i]) 
			scr_batch <- c(scr_batch,scr_batch_tmp)
			scr_scdb  <- c(scr_scdb,scr_scdb_tmp)			
		}
	}
	no_ercc <- grep("ERCC", rownames(umis), invert=T)
	umis    <- umis[no_ercc, -1]
	f       <- colSums(umis)>min_umis_n
	umis    <- umis[,f]
	f_cells <- sum(!f)
	return(list(umis=umis,
        umis_path=files_list,
        scr_batch=scr_batch,
        scr_scdb=scr_scdb,
        unmatched_rownames_n=rnames_miss_total,
        filtered_cells_n=f_cells))
}

#' Performs downsampling in a given cell
#'
#' Brief explanation of the implementation:
#' In .downsamp_one, rep(1:length(v),times=v) repeat each gene index times the number of umis found in that gene
#' in order to sample n elements from it (usign sample(rep(1:length(v),times=v),replace=F,size=n)). 
#' Then, hist(sample(rep(1:length(v),times=v),replace=F,size=n),0.5+0:length(v),plot=F)$counts outputs how many
#' times a gene index was sampled (Breaks are used to fit gene indexes, so it counts umis sampled per gene). 
#' That is, the resulting downsampled cell.
#'
#' @param v (numeric vector) cell column in a umis table
#' @param n (integer) target number of molecules to sample from v without replacement
#' @param normalize normalize using a division by the sum of each cell
#' 
#' @return downsampled vector
#'
#' @examples
#'
.downsamp_one <- function(v, n){
  hist(sample(rep(1:length(v),times=v),replace=F,size=n),0.5+0:length(v),plot=F)$counts
}

#' Downsample and normalize umis
#'
#' @param scram_obj scram object with umis field that contains umis matrix (numeric matrix of umis with genes as rownames and cells as columns)
#' @param dsamp_n target number of umis per each cell
#' @param normalize normalize using a division by the sum of each cell
#' 
#' @return scram_obj with the following additional fields: 
#' \itemize{
#'  \item{"scr_umis_ds"}{downsampled UMI matrix}
#'  \item{"scr_umis_norm"}{downsampled and normalized UMI matrix}
#'  \item{"target_dsamp_n"}{downsampling N}
#' }
#'
#' @examples
#'
#' @export
scram.downsamp <- function(scram_obj, dsamp_n, normalize=TRUE){
    scr_umis <- scram_obj$umis
	scr_umis_ds   <- apply(scr_umis[, colSums(scr_umis)>=dsamp_n], 2, .downsamp_one, dsamp_n)
	rownames(scr_umis_ds)   <- rownames(scr_umis)
	if(normalize){
		scr_umis_norm <- round(1000*t(t(scr_umis)/colSums(scr_umis)),2)
		rownames(scr_umis_norm) <- rownames(scr_umis)
	} else { 
		scr_umis_norm=NULL
	}
    
    scram_obj <- c(scram_obj, list(
            scr_umis_ds = scr_umis_ds,
            scr_umis_norm = scr_umis_norm,
            target_dsamp_n = dsamp_n))
            
	return(scram_obj)
}
