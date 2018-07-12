SHELL=/bin/bash -x

SAMPLES=mono
OUTPUTDIR=$(addprefix processeddata/, $(SAMPLES))

trims: $(addsuffix /trim.fq, $(OUTPUTDIR))
norrnas: $(addsuffix /norrna.fq, $(OUTPUTDIR))
sort: $(addsuffix /accepted_hits.bam.bai, $(OUTPUTDIR))

CUTADAPT=/fh/fast/subramaniam_a/user/rasi/virtualenv/default2/bin/cutadapt
RIBO_CUTADAPT_ARGS=--cut=8 --adapter=TCGTATGCCGTCTTCTGCTTG --minimum-length=15
BOWTIE=/app/bowtie/1.1.1/bowtie
RRNA_BOWTIE_ARGS=--seedlen=23 --threads=8
RRNA_BOWTIE_INDEX=/fh/fast/subramaniam_a/db/rasi/bowtie/sacCer3.contaminants
TOPHAT=/app/tophat/2.0.14/bin/tophat
GENOME_TOPHAT_ARGS=--bowtie1 --num-threads 8 --no-novel-juncs --library-type fr-unstranded --keep-tmp --read-mismatches 2
GENOME_BOWTIE_INDEX=/fh/fast/subramaniam_a/db/rasi/bowtie/sacCer3.genome
TRANSCRIPTOME_BOWTIE_INDEX=/fh/fast/subramaniam_a/db/rasi/bowtie/sacCer3.sgd.exons
SAMTOOLS=/app/samtools/1.0/bin/samtools

# add bowtie to path so that tophat can find it
export PATH

## all: run analysis to generate all final files
all: trims
all: norrnas
all: sort

## clean: removes all processed data files, leaves figures intact
clean:	
	rm --force annotations/*
	rm --force rawdata/*.sra
	rm --force processeddata/*/trim.fq
	rm --force processeddata/*/norrna.fq

## help: displays this help
help: Makefile
	@sed -n 's/^##//p' $<

# prevents error in case files are named with these keywords
.PHONY: clean all help

# 5' and 3' trimming
processeddata/%/trim.fq: rawdata/%.fastq.gz
	mkdir -p processeddata/$*
	$(CUTADAPT) $(RIBO_CUTADAPT_ARGS) --output $@ $< 1> $(addsuffix .log,$@) 2> $(addsuffix .log,$@)

# align against rrna contaminants
processeddata/%/norrna.fq: processeddata/%/trim.fq
	$(BOWTIE) \
	$(RRNA_BOWTIE_ARGS) --un $@ --sam $(RRNA_BOWTIE_INDEX) $< \
	2>$(addsuffix .log,$@) \
	> /dev/null

# align against genome with pre-built transcriptome index
processeddata/%/tophat/accepted_hits.bam: processeddata/%/norrna.fq
	mkdir -p processeddata/$*/tophat
	$(TOPHAT) \
	$(GENOME_TOPHAT_ARGS) \
	--transcriptome-index=$(TRANSCRIPTOME_BOWTIE_INDEX) \
	--output-dir=processeddata/$*/tophat \
	$(GENOME_BOWTIE_INDEX) $< \
	2>$(addsuffix .log,$@)
	cat processeddata/$*/tophat/align_summary.txt > $(addsuffix .log, $@)
	rm processeddata/$*/tophat/align_summary.txt

# index alignments and move alignment and index to parent folder
processeddata/%/accepted_hits.bam.bai: processeddata/%/tophat/accepted_hits.bam
	$(SAMTOOLS) index $< 2>$(addsuffix .log,$@)
	mv $< $(subst tophat/,,$<)
	mv $(addsuffix .bai,$<) $@
