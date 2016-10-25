#!/usr/bin/env Rscript
suppressPackageStartupMessages(library('dplyr'))
library('optparse')

fread <- function(...) as.data.frame(data.table::fread(...), stringsAsFactors=FALSE)


########################################################################
main <- function(argv)
{
    command <- parse_command(argv)
    
    intervals <- fread(command$intervals, sep='\t', header=TRUE)
    intervals <- intervals %>% filter(substr(chrom, 1, 4) != 'ERCC')
    
    gene_names <- intervals %>% 
                  group_by(gene_name, chrom) %>%
                  do(split_name(.$gene_name[1])) %>%
                  select(gene_name, geneSymbol, chrom) %>%
                  ungroup()
          
    write.table(gene_names, command$output, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
} 


########################################################################
parse_command <- function(argv)
{
    options_list <- list(
    )
    
    parser  <- OptionParser(usage='%prog  [OPTIONS]  <gene_intervals.txt>  <output.txt>', option_list=options_list)
    parsed  <- parse_args(parser, argv, positional_arguments=TRUE)
    
    command <- list();
    
    command$intervals <- parsed$args[1] 
    command$output    <- parsed$args[2]
    
    
    if (is.na(command$intervals) || is.na(command$output)) {
        stop('The names of gene_intervals and gene_names files must be specified.')
    }
    
    return(command);
}


########################################################################
split_name <- function(gene_name)
{
    names <- strsplit(gene_name, ';', fixed=TRUE)[[1]]
    return(data.frame(gene_name=gene_name, geneSymbol=names, stringsAsFactors=FALSE))
}


########################################################################
main(commandArgs(TRUE))