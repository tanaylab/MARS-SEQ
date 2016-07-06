ANNOTATIONS := gene_intervals_mm10.txt

MM10_ROOT := /net/mraid14/export/data/db/tgdb/mm10/trackdb
MM10_UCSC := /net/mraid14/export/data/db/tgdb/mm10/rawdata/UCSC

ERCC_INTERVALS := gene_intervals_ercc.txt

.PHONY: all clean

all: $(ANNOTATIONS)

clean:
	rm -f $(ANNOTATIONS)

gene_intervals_mm10.txt: $(MM10_UCSC)/knownGene.txt.gz $(MM10_UCSC)/kgXref.txt.gz $(ERCC_INTERVALS)
	./gene_intervals.R --root $(MM10_ROOT) --ucsc $(MM10_UCSC) $@
	cat $(ERCC_INTERVALS) >> $@