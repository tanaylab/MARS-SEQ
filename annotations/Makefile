ANNOTATIONS := gene_intervals_mm9.txt \
               gene_intervals_mm10.txt \
               gene_names_mm9.txt \
               gene_names_mm10.txt \
               gene_tss_mm9.txt \
               gene_tss_mm10.txt

MM9_ROOT := /net/mraid14/export/data/db/tgdb/mm9/trackdb
MM9_UCSC := /net/mraid14/export/data/db/tgdb/mm9/rawdata/UCSC

MM10_ROOT := /net/mraid14/export/data/db/tgdb/mm10/trackdb
MM10_UCSC := /net/mraid14/export/data/db/tgdb/mm10/rawdata/UCSC

ERCC_INTERVALS := gene_intervals_ercc.txt

.PHONY: all clean

all: $(ANNOTATIONS)

clean:
	rm -f $(ANNOTATIONS)
	
# TSS table uses mm10 names but mm9 positions, as energy/ChIP tracks are currently mostly at mm9
gene_tss_mm10.txt: gene_names_mm10.txt
	./gene_tss.R --ucsc $(MM9_UCSC) $< $@

gene_tss_mm9.txt: gene_names_mm9.txt
	./gene_tss.R --ucsc $(MM9_UCSC) $< $@

gene_names_mm9.txt : gene_intervals_mm9.txt
	./gene_names.R $< $@
	
gene_names_mm10.txt : gene_intervals_mm10.txt
	./gene_names.R $< $@
	
gene_intervals_mm9.txt: $(MM9_UCSC)/knownGene.txt.gz $(MM9_UCSC)/kgXref.txt.gz $(ERCC_INTERVALS)
	./gene_intervals.R --root $(MM9_ROOT) --ucsc $(MM9_UCSC) $@
	cat $(ERCC_INTERVALS) >> $@	

gene_intervals_mm10.txt: $(MM10_UCSC)/knownGene.txt.gz $(MM10_UCSC)/kgXref.txt.gz $(ERCC_INTERVALS)
	./gene_intervals.R --root $(MM10_ROOT) --ucsc $(MM10_UCSC) $@
	cat $(ERCC_INTERVALS) >> $@