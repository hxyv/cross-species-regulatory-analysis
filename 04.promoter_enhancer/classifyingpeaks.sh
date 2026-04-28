#!/bin/bash

# Load bedtools so we can compare genomic regions
module load bedtools
#we want to classify peaks  we have about 206765 human ATAC peaks.
# if this peak is near  a tss we will classify it as prmoter if not then it is an enhancer
echo "Classifying human peaks..."
bedtools intersect -a human_adrenal_idr_optimal.narrowPeak -b human_tss_2kb.bed -u > humanpromoters.bed
bedtools intersect -a human_adrenal_idr_optimal.narrowPeak -b human_tss_2kb.bed -v > humanenhancers.bed
#now we will classify mouse peaks as promoters or enhancers
echo "Classifying mouse peaks..."
bedtools intersect -a mouse_adrenal_idr_optimal.narrowPeak -b mouse_tss_2kb.bed -u > mouse_promoters.bed
bedtools intersect -a mouse_adrenal_idr_optimal.narrowPeak -b mouse_tss_2kb.bed -v > mouse_enhancers.bed


echo "Human promoters:"; wc -l humanpromoters.bed
echo "Human enhancers:"; wc -l humanenhancers.bed
echo "Mouse promoters:"; wc -l mouse_promoters.bed
echo "Mouse enhancers:"; wc -l mouse_enhancers.bed
# from xingyu's task 2 we know which peaks are shared and specific
#now we want to know among conserved how many are  promoters vs enhancers
#we want to know how  many from the species specific are promoters vs enhacers
echo "Checking conservation..."
# Conserved peaks which are human peaks that have a mouse counterpart
bedtools intersect -a human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed -b mouse_tss_2kb.bed -u > conservedpromoters.bed
bedtools intersect -a human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed -b mouse_tss_2kb.bed -v > conservedenhancers.bed

# Human-specific peaks which are human peaks with NO mouse counterpart
bedtools intersect -a human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed -b mouse_tss_2kb.bed -u > humanspecific_promoters.bed
bedtools intersect -a human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed -b mouse_tss_2kb.bed -v > humanspecific_enhancers.bed

# Mouse-specific peaks which are mouse peaks with NO human counterpart
bedtools intersect -a mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed -b mouse_tss_2kb.bed -u > mousespecific_promoters.bed
bedtools intersect -a mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed -b mouse_tss_2kb.bed -v > mousespecific_enhancers.bed

echo ""
echo "Conserved promoters:"; wc -l conservedpromoters.bed
echo "Conserved enhancers:"; wc -l conservedenhancers.bed
echo "Human-specific promoters:"; wc -l humanspecific_promoters.bed
echo "Human-specific enhancers:"; wc -l humanspecific_enhancers.bed

# we want HUMAN-SPECIFIC IN HUMAN COORDINATES (hg38)
# Using Xingyu's recovered files with original human coordinates
# Classified against human TSS

bedtools intersect -a human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed -b human_tss_2kb.bed -u > human_specific_promoters_hg38.bed
bedtools intersect -a human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed -b human_tss_2kb.bed -v > human_specific_enhancers_hg38.bed

# CONSERVED IN HUMAN COORDINATES (hg38)
bedtools intersect -a human_adrenal_idr_optimal.shared.original_human_coordinates.bed -b human_tss_2kb.bed -u > conserved_promoters_hg38.bed
bedtools intersect -a human_adrenal_idr_optimal.shared.original_human_coordinates.bed -b human_tss_2kb.bed -v > conserved_enhancers_hg38.bed

echo ""
echo "=== Human coordinates (hg38) ==="
echo "Human-specific promoters (hg38):"; wc -l human_specific_promoters_hg38.bed
echo "Human-specific enhancers (hg38):"; wc -l human_specific_enhancers_hg38.bed
echo "Conserved promoters (hg38):"; wc -l conserved_promoters_hg38.bed
echo "Conserved enhancers (hg38):"; wc -l conserved_enhancers_hg38.bed


# we want MOUSE CONSERVED IN MOUSE COORDINATES (mm10)
# Using native mouse conserved peaks
# Classified against mouse TSS

bedtools intersect -a mouse_adrenal_idr_optimal.shared_with_human_mapped.bed -b mouse_tss_2kb.bed -u > mouse_conserved_promoters.bed
bedtools intersect -a mouse_adrenal_idr_optimal.shared_with_human_mapped.bed -b mouse_tss_2kb.bed -v > mouse_conserved_enhancers.bed

echo ""
echo "=== Mouse conserved (mm10) ==="
echo "Mouse conserved promoters:"; wc -l mouse_conserved_promoters.bed
echo "Mouse conserved enhancers:"; wc -l mouse_conserved_enhancers.bed


echo "Mouse-specific promoters:"; wc -l mousespecific_promoters.bed
echo "Mouse-specific enhancers:"; wc -l mousespecific_enhancers.bed



