#!/bin/bash
#unzip files from gencode
gunzip -f gencode.v45.basic.annotation.gtf.gz 2>/dev/null
gunzip -f gencode.vM25.basic.annotation.gtf.gz 2>/dev/null

# Extract human TSS with 2kb threshold
awk '$3 == "gene"' gencode.v45.basic.annotation.gtf | awk -v OFS="\t" '{if ($7 == "+") tss = $4; else tss = $5; start = tss - 2000; if (start < 0) start = 0; end = tss + 2000; print $1, start, end;}' | sort -k1,1 -k2,2n > human_tss_2kb.bed

# Extract mouse TSS with 2kb threshold
awk '$3 == "gene"' gencode.vM25.basic.annotation.gtf | awk -v OFS="\t" '{if ($7 == "+") tss = $4; else tss = $5; start = tss - 2000; if (start < 0) start = 0; end = tss + 2000; print $1, start, end;}' | sort -k1,1 -k2,2n > mouse_tss_2kb.bed

# Check  if we are right
wc -l human_tss_2kb.bed mouse_tss_2kb.bed
head -3 human_tss_2kb.bed
head -3 mouse_tss_2kb.bed
