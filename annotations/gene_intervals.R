#!/usr/bin/env Rscript
suppressPackageStartupMessages(library('dplyr'))
library('optparse')
library('misha')


########################################################################
KNOWN_GENE_COLS <- c('name',
                     'chrom',
                     'strand',
                     'txStart',
                     'txEnd',
                     'cdsStart',
                     'cdsEnd',
                     'exonCount',
                     'exonStarts',
                     'exonEnds',
                     'proteinID',
                     'alignID')

KG_XREF_COLS <- c('kgID',
                  'mRNA',
                  'spID',
                  'spDisplayID',
                  'geneSymbol',
                  'refseq',
                  'protAcc',
                  'description',
                  'rfamAcc',
                  'tRnaName')


########################################################################
main <- function(argv)
{
    command <- parse_command(argv)
    
    gsetroot(command$root)
    
    # Read USCS tables
    knownGene <- fread(qq('gzip -d -c @{command$ucsc}/knownGene.txt.gz'),
                       sep='\t', header=FALSE, col.names=KNOWN_GENE_COLS)
       
    kgXref <- fread(qq('gzip -d -c @{command$ucsc}/kgXref.txt.gz'),
                    sep='\t', select=1:8, header=FALSE, col.names=KG_XREF_COLS[1:8])      
    

    # Remove non-chromosome contigs and correct strand
    knownGene <- knownGene %>%
                 filter(!grepl('_', chrom, fixed=TRUE)) %>%
                 mutate(strand=ifelse(strand=='+', +1, -1))        

    # Add gene symbol
    kgXref <- kgXref %>% select(kgID, geneSymbol)
    knownGene <- knownGene %>%
                 select(name, chrom, strand, exonStarts, exonEnds) %>%
                 left_join(kgXref, by=c('name'='kgID'))     

    # Expand to one-line per exon
    exon_start <- strsplit(knownGene$exonStarts, ',')
    exon_end <- strsplit(knownGene$exonEnds, ',')
    exon_tab <- data.frame(name=rep(knownGene$name, sapply(exon_start, length)), 
                           start=as.numeric(unlist(exon_start)), 
                           end=as.numeric(unlist(exon_end)),
                           stringsAsFactors=FALSE)
    
    knownGene <- knownGene %>% 
                 left_join(exon_tab, by='name') %>%
                 select(chrom, start, end, strand, gene_name=geneSymbol)
         
    # Merge overlapping exons
    knownGene_pos <- knownGene %>% 
                     filter(strand==+1) %>%
                     merge_genes() %>% 
                     mutate(strand=+1) %>%
                     select(chrom, start, end, strand, gene_name)
    knownGene_neg <- knownGene %>% 
                     filter(strand==-1) %>%
                     merge_genes() %>% 
                     mutate(strand=-1) %>%
                     select(chrom, start, end, strand, gene_name)
    
    knownGene <- rbind(knownGene_pos, knownGene_neg)
    
    # Write table to disk
    write.table(knownGene, file=command$output, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
}

########################################################################
parse_command <- function(argv)
{
    options_list <- list(
            make_option( c('--root', '-r'), type='character', dest='root',
                    help="Genomic database's root directory." ),
            
            make_option( c('--ucsc', '-u'), type='character', dest='ucsc',
                    help="Directory holding UCSC tables (knownGene.txt.gz and kgXref.txt.gz)." )
    )
    
    parser  <- OptionParser(usage='%prog  [OPTIONS]  <intervals_file.txt>', option_list=options_list)
    parsed  <- parse_args(parser, argv, positional_arguments=TRUE)
    
    command <- list();
    
    command$root   <- parsed$options$root
    command$ucsc   <- parsed$options$ucsc
    command$output <- parsed$args[1]
    
    
    if (is.null(command$root) || is.null(command$ucsc) || is.null(command$output)) {
        print_help(parser)
        q()
    }    
    
    return(command);
}


########################################################################
fread <- function(...)
{
    return(as.data.frame(data.table::fread(...)))
}


########################################################################
qq <- GetoptLong::qq 


########################################################################
merge_genes <- function(gene_tab)
{
    # Crete canonic intervals
    merged_tab <- gintervals.canonic(gene_tab)
    canonic <- attr(merged_tab, 'mapping')
    
    gene_name <- gene_tab %>%
                 select(gene_name) %>%
                 mutate(canonic=attr(merged_tab, 'mapping'))

    # Map gene names to all other names that overlap in a canonic interval
    name_map <- gene_name %>%
                rename(target_name=gene_name) %>%
                right_join(gene_name, by='canonic') %>%
                group_by(gene_name, target_name) %>%
                summarize() %>%
                ungroup()
    
    while(TRUE) {
        before <- nrow(name_map)
        name_map <- name_map %>%
                    left_join(name_map, by=c('target_name'='gene_name')) %>%
                    select(gene_name, target_name=target_name.y) %>%
                    group_by(gene_name, target_name) %>%
                    summarize() %>%
                    ungroup()
        after <- nrow(name_map)
        if (after == before) {
            break                    
        }
    }
        
    # Map each gene name to a semicolon sepearated list of overlapping genes
    name_map <- name_map %>%
                group_by(gene_name) %>% 
                summarize(target_name=paste0(sort(unique(target_name)), collapse=';')) %>%
                ungroup()
    
    # Add list of gene names tp canonic intervals
    merged_tab <- merged_tab %>%
                  mutate(canonic=row_number()) %>%
                  left_join(gene_name, by='canonic') %>%
                  left_join(name_map, by='gene_name') %>%
                  group_by(chrom, start, end) %>%
                  summarize(gene_name=paste0(unique(target_name), collapse='--')) %>%
                  ungroup()
    
    # Order genes based on start position, maintaining exons together
    merged_tab$group <- group_indices(merged_tab, chrom, gene_name)
    order_tab <- merged_tab %>%
                 group_by(chrom, gene_name, group) %>%
                 summarize(start=min(start)) %>%
                 ungroup() %>%
                 arrange(chrom, start) %>%
                 mutate(order=row_number()) %>%
                 select(group, order) %>%
                 right_join(merged_tab, by='group') %>%
                 arrange(order, chrom, start) %>%
                 select(chrom, start, end, gene_name)

    return(merged_tab)
}


########################################################################
main(commandArgs(TRUE))