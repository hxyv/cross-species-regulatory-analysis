** Chester did the following work, some of the git commit and changes were made on Xinyu's computer and account **
# We are choosing Adrenal Gland for both human and mouse for the following reasons:

No data’s NFR passed 0.9

All of data passed(\>94)% properly paired reads test
but Mouse and Human ovary data demonstrate poor periodicity when
comparing to Mouse and Human Adrenal Gland data:

## Mouse ovary

**% Properly Paired Reads:** rep1-94.4, rep2-94.6

![](assets/media/image1.png)

![](assets/media/image3.png)

![](assets/media/image2.png)

Mouse ovary Fragment length statistics (filtered/deduped BAM) plot only shows one peak in both rep, 
lack distinct oscillations and are dominated by short fragments; Even though Both TSS(Transcription start site) enrichment passed 13, 
the filter reads both passed 25 millions(27.41M - rep1, 72.27M - rep2)

## Human ovary

**% Properly Paired Reads:** rep1-96.2, rep2-96.6

![](assets/media/image6.png)

![](assets/media/image4.png)

![](assets/media/image5.png)

Human ovary Fragment length statistics (filtered/deduped BAM) plot lacks distinct oscillations, 
only shows weak periodicity comparing to Adrenal Gland results; rep1’s TSS enrichment did not pass 13 (6.6), 
though the filter reads both passed 25 millions(97.04M - rep1, 275.91M - rep2)

## Mouse AD

**% Properly Paired Reads:** rep1-94.6, rep2-94.7

![](assets/media/image7.png)

![](assets/media/image8.png)

![](assets/media/image9.png)

Compared to the mouse ovary, the fragment length statistics (filtered/deduped BAM) exhibit significantly better periodicity. 
Instead of being dominated by short fragments, the fragment length distributions in both replicates show a more distinct oscillations. 
Additionally, both replicates exhibit high TSS enrichment (18.37 in rep1 and 18.998 in rep2). 
Both duplicates exhibit noticeably higher signal quality than the ovary sample, 
even while the filtered read depth in rep1 is below the optimal 25 million criterion (17.67M), 
while rep2 passes this barrier (27.10M). The mouse adrenal gland dataset is generally of higher quality.

## Human AD

**% Properly Paired Reads:** rep1-96.0, rep2-94.3

![](assets/media/image10.png)

![](assets/media/image11.png)

![](assets/media/image12.png)

Both replicates exhibit significant periodicity in the fragment length statistics (filtered/deduped BAM). 
Instead of being dominated by short fragments, the fragment length distributions in both replicates show a more distinct oscillations. 
Additionally, both replicates exhibit high TSS enrichment (19.19 in rep1 and 20.39 in rep2). 
In addition, both replicates have more than sufficient filtered read depth (364.04M for rep1 and 101.58M for rep2), 
far exceeding the 25 million guideline. When taken as a whole, these findings show that the human adrenal gland dataset is significantly better than the human ovary dataset.

## Poor Periodicity Discussion

Ovary libraries from both humans and mice exhibited weak periodicity, with short fragments dominating fragment length distributions and lacking the distinct mono-/di-nucleosome oscillation typical of high-quality ATAC-seq libraries[1]. Weaker chromatin or nuclei preservation, increased background from damaged cells or free DNA, over-tagmentation, and potentially lower library complexity are all the possible reasons for this pattern[1][2]. The issue seems to be worse in the human ovary, where rep1 also failed the TSS enrichment threshold(13)[2]. These interpretations are in line with the course QC guidelines, the ATAC-seq technique, and the literature that demonstrates that high background, decreased TSS enrichment, and poor nucleosomal patterning frequently indicate subpar sample quality or excessive transposase[2].

The adrenal gland datasets, on the other hand, demonstrated more TSS enrichment and clearer periodicity, making them comparatively more dependable for further investigation. Nevertheless, since none of the datasets met the suggested NRF standard, suggesting that library complexity was not perfect overall, adrenal gland should still be characterized as comparatively superior rather than ideal[3].

## Reference

[1] J. D. Buenrostro, P. G. Giresi, L. C. Zaba, H. Y. Chang, and W. J. Greenleaf, “Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position,” Nature Methods, vol. 10, no. 12, pp. 1213–1218, Dec. 2013, doi: 10.1038/nmeth.2688.

[2] M. R. Corces et al., “An improved ATAC-seq protocol reduces background and enables interrogation of frozen tissues,” Nature Methods, vol. 14, no. 10, pp. 959–962, Oct. 2017, doi: 10.1038/nmeth.4396.

[3] ENCODE Project Consortium, “ATAC-seq Data Standards and Processing Pipeline,” ENCODE, accessed Apr. 6, 2026.
