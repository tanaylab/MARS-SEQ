#!/usr/bin/env Rscript
library('optparse')
library('misha')
source('../../../analysis/common/tgverse.r', chdir=TRUE)
source('../../../analysis/common/ucsc.r', chdir=TRUE)

########################################################################
main <- function(argv)
{
    command <- parse_command(argv)
    
    gsetroot(command$root)
    
    # Read USCS tables
    knownGene <- read_known_gene(file.path(command$ucsc, 'knownGene.txt.gz'))
    kgXref <- read_kg_xref(file.path(command$ucsc, 'kgXref.txt.gz'))

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
    exon_tab <- tibble(name=rep(knownGene$name, sapply(exon_start, length)), 
                       start=as.numeric(unlist(exon_start)), 
                       end=as.numeric(unlist(exon_end)))
    
    knownGene <- knownGene %>% 
                 left_join(exon_tab, by='name') %>%
                 select(chrom, start, end, strand, gene_name=geneSymbol) %>%
                 distinct()
         
    # Merge overlapping exons
    knownGene_pos <- knownGene %>% 
                     filter(strand==+1) %>%
                     merge_genes() %>% 
                     select(chrom, start, end, strand, gene_name)
    knownGene_neg <- knownGene %>% 
                     filter(strand==-1) %>%
                     merge_genes() %>% 
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
    
    parser  <- OptionParser(usage='%prog  [OPTIONS]  <output.txt>', option_list=options_list)
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
merge_genes <- function(gene_tab)
{
    # Match each interval into its canonic interval
    merged_tab <- gintervals.canonic(gene_tab)
    gene_tab <- gene_tab %>% 
                mutate(canonic=attr(merged_tab, 'mapping'))
    merged_tab <- merged_tab %>%
                  mutate(canonic=row_number())
          
    # Keep one representative of each (canonic, gene) pairs
    pairs <- gene_tab %>%
             select(gene_name, canonic) %>%
             distinct()

    # Generates all pairs of genes that share the same canonic interval
    pairs <- pairs %>%
             group_by(canonic) %>%
             filter(n() > 1) %>%
             do(extract_pairs(.$gene_name)) %>%
             ungroup() %>%
             select(gene1, gene2) %>%
             distinct()
    
    # Convert to an adjacency matrix
    names <- unique(c(pairs$gene1, pairs$gene2)) %>% 
             sort()
    adjacency <- matrix(data=0, nrow=length(names), ncol=length(names))
    colnames(adjacency) <- names
    rownames(adjacency) <- names
    adjacency[cbind(pairs$gene1, pairs$gene2)] <- 1
    adjacency <- adjacency + t(adjacency)
    diag(adjacency) <- 1
    
    # Calculate exponentially increasing powers of the adjacency table until all
    # Higher order adjacencies are found
    last_count <- sum(adjacency != 0)
    while(TRUE) {
        adjacency <- adjacency %*% adjacency
        adjacency[adjacency != 0] <- 1
        count <- sum(adjacency != 0)
        if (count == last_count) {
            break
        }
        last_count <- count
    }
    
    # Convert gene names into names of overlapping groups
    adjacency <- adjacency > 0
    name_map <- apply(adjacency, 1, function(r) paste0(names(r)[r], collapse=';'))
    name_map <- tibble(gene_name=names(name_map), target_name=name_map)
    gene_tab <- gene_tab %>%
                left_join(name_map, by='gene_name') %>%
                mutate(gene_name=ifelse(is.na(target_name), gene_name, target_name)) %>%
                select(-target_name)
        
    # Make sure that each canonic interval now has a single name
    gene_tab <- gene_tab %>%
                group_by(canonic, strand) %>%
                summarize(count=length(unique(gene_name)), gene_name=first(gene_name)) %>%
                ungroup()
    stopifnot(nrow(gene_tab %>% filter(count > 1))  == 0)        
    
    # Canonicize the intervals of each group
    gene_tab <- gene_tab %>%
                select(-count)
    gene_tab <- merged_tab %>%
                left_join(gene_tab, by='canonic')
    
    # Order genes based on start position, maintaining exons together
    gene_tab <- gene_tab %>%
                select(chrom, start, end, strand, gene_name) %>%
                group_by(chrom, gene_name) %>%
                mutate(gene_start=min(start)) %>%
                ungroup() %>%
                arrange(chrom, gene_start, gene_name) %>%
                select(-gene_start)

    return(gene_tab)
}


########################################################################
extract_pairs <- function(names)
{
    pairs <- combn(names, 2)
    return(tibble(gene1=pairs[1,], gene2=pairs[2,]))
}


########################################################################
main(commandArgs(TRUE))