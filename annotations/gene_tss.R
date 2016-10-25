#!/usr/bin/env Rscript
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
source('../../../analysis/common/ucsc.r', chdir=TRUE)

fread <- function(...) as.data.frame(data.table::fread(...), stringsAsFactors=FALSE)


########################################################################
main <- function(argv)
{
    command <- parse_command(argv)
    
    known_gene <- read_known_gene(file.path(command$ucsc, 'knownGene.txt.gz'))
    kg_xref <- read_kg_xref(file.path(command$ucsc, 'kgXref.txt.gz'))
    
    gene_names <- fread(command$names, sep='\t', header=TRUE)
    
    gene_tss <- gene_names %>% 
                inner_join(kg_xref, by='geneSymbol') %>%
                inner_join(known_gene, by=c(kgID='name', 'chrom')) %>%
                mutate(end=txStart+1) %>%
                select(chrom, start=txStart, end, strand, gene_name) %>%
                distinct() %>%
                mutate(strand=recode(strand, '+'=+1, '-'=-1))
        
    write.table(gene_tss, command$output, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
} 



########################################################################
parse_command <- function(argv)
{
    options_list <- list(
            make_option( c('--ucsc', '-u'), type='character', dest='ucsc',
                    help="Directory holding UCSC tables (knownGene.txt.gz and kgXref.txt.gz)." )
    )
    
    parser  <- OptionParser(usage='%prog  [OPTIONS]  <gene_names.txt>  <gene_names.txt>', option_list=options_list)
    parsed  <- parse_args(parser, argv, positional_arguments=TRUE)
    
    command <- list();
    
    command$ucsc   <- parsed$options$ucsc
    command$names  <- parsed$args[1] 
    command$output <- parsed$args[2]
    
    
    if (is.na(command$names) || is.na(command$output)) {
        stop('The names of the gene_names and output files must be specified.')
    }
    
    return(command);
}


########################################################################
main(commandArgs(TRUE))