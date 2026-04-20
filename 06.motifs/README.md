# Motif Enrichment

This folder contains Task 5 motif enrichment analyses using HOMER.

Main workflow:

1. Task 4 promoter/enhancer peak sets are staged into `downstream/`.
2. `downstream/build_centered_beds.sh` creates 200 bp summit-centered peak windows.
3. `downstream/run_all.sh` runs the HOMER analyses for human and mouse, centered and non-centered windows, and all-peak backgrounds.
4. `downstream/summarize_known_motifs.py` summarizes HOMER `knownResults.txt` files into `downstream/results_summary/`.

The HOMER installation itself is not stored in the repository. The full pipeline can link an external HOMER installation with the `--homer /path/to/homer` argument.

See `HOMER_Installation_Log_and_Setup_Guide.md` for setup notes.

The Results can be found below from summarize_known_motifs.py:

# HOMER Known Motif Summary — All 20 Experiments
# Filter: q-value < 0.05, delta_pct > 5.0
# Top 10 shown per experiment per ranking

## out_allHumanEnh_hg38_defaultBg
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 26

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 42.6% | 13.5% | +29.1 |
| 2 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 38.9% | 11.3% | +27.6 |
| 3 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 44.1% | 14.5% | +29.6 |
| 4 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 46.5% | 16.6% | +29.9 |
| 5 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0000 | 41.5% | 13.5% | +28.1 |
| 6 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 45.4% | 16.3% | +29.1 |
| 7 | Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer | 0.0000 | 30.0% | 7.4% | +22.6 |
| 8 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 48.2% | 18.7% | +29.4 |
| 9 | Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer | 0.0000 | 23.9% | 5.2% | +18.6 |
| 10 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 14.4% | 1.9% | +12.6 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 46.5% | 16.6% | +29.9 |
| 2 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 44.1% | 14.5% | +29.6 |
| 3 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 48.2% | 18.7% | +29.4 |
| 4 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 45.4% | 16.3% | +29.1 |
| 5 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 42.6% | 13.5% | +29.1 |
| 6 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0000 | 41.5% | 13.5% | +28.1 |
| 7 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 38.9% | 11.3% | +27.6 |
| 8 | Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer | 0.0000 | 30.0% | 7.4% | +22.6 |
| 9 | Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer | 0.0000 | 37.1% | 18.1% | +18.9 |
| 10 | Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer | 0.0000 | 23.9% | 5.2% | +18.6 |

---

## out_allHumanProm_hg38_defaultBg
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 39

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 27.4% | 13.9% | +13.4 |
| 2 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 24.0% | 11.7% | +12.2 |
| 3 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 29.4% | 15.9% | +13.5 |
| 4 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 26.0% | 13.5% | +12.4 |
| 5 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0000 | 25.8% | 13.5% | +12.4 |
| 6 | Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer | 0.0000 | 18.6% | 8.3% | +10.2 |
| 7 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 31.4% | 18.0% | +13.4 |
| 8 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 28.1% | 15.5% | +12.6 |
| 9 | Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer | 0.0000 | 14.6% | 6.0% | +8.5 |
| 10 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 13.4% | 5.5% | +7.9 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 29.4% | 15.9% | +13.5 |
| 2 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 27.4% | 13.9% | +13.4 |
| 3 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 31.4% | 18.0% | +13.4 |
| 4 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 28.1% | 15.5% | +12.6 |
| 5 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 26.0% | 13.5% | +12.4 |
| 6 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0000 | 25.8% | 13.5% | +12.4 |
| 7 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 24.0% | 11.7% | +12.2 |
| 8 | Elk4(ETS)/Hela-Elk4-ChIP-Seq(GSE31477)/Homer | 0.0000 | 37.5% | 26.3% | +11.2 |
| 9 | Sp1(Zf)/Promoter/Homer | 0.0000 | 36.0% | 24.8% | +11.2 |
| 10 | Elk1(ETS)/Hela-Elk1-ChIP-Seq(GSE31477)/Homer | 0.0000 | 37.3% | 26.4% | +10.8 |

---

## out_allMouseEnh_mm10_defaultBg
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 51

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 21.1% | 3.6% | +17.5 |
| 2 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 51.0% | 24.4% | +26.6 |
| 3 | BORIS(Zf)/K562-CTCFL-ChIP-Seq(GSE32465)/Homer | 0.0000 | 22.9% | 5.8% | +17.2 |
| 4 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 54.8% | 28.6% | +26.3 |
| 5 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 52.8% | 27.4% | +25.4 |
| 6 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 66.9% | 40.5% | +26.3 |
| 7 | EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Hom | 0.0000 | 38.6% | 16.6% | +22.0 |
| 8 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 62.0% | 36.7% | +25.2 |
| 9 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 47.7% | 24.4% | +23.3 |
| 10 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 54.6% | 30.4% | +24.2 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 51.0% | 24.4% | +26.6 |
| 2 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 66.9% | 40.5% | +26.3 |
| 3 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 54.8% | 28.6% | +26.3 |
| 4 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 52.8% | 27.4% | +25.4 |
| 5 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 62.0% | 36.7% | +25.2 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 54.6% | 30.4% | +24.2 |
| 7 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 47.7% | 24.4% | +23.3 |
| 8 | EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Hom | 0.0000 | 38.6% | 16.6% | +22.0 |
| 9 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 49.0% | 27.6% | +21.4 |
| 10 | EHF(ETS)/LoVo-EHF-ChIP-Seq(GSE49402)/Homer | 0.0000 | 53.1% | 32.7% | +20.4 |

---

## out_allMouseProm_mm10_defaultBg
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 56

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Elk4(ETS)/Hela-Elk4-ChIP-Seq(GSE31477)/Homer | 0.0000 | 49.2% | 25.4% | +23.8 |
| 2 | Sp1(Zf)/Promoter/Homer | 0.0000 | 49.8% | 26.0% | +23.8 |
| 3 | Elk1(ETS)/Hela-Elk1-ChIP-Seq(GSE31477)/Homer | 0.0000 | 48.1% | 25.4% | +22.8 |
| 4 | ETS(ETS)/Promoter/Homer | 0.0000 | 31.1% | 12.7% | +18.5 |
| 5 | ELF1(ETS)/Jurkat-ELF1-ChIP-Seq(SRA014231)/Homer | 0.0000 | 44.4% | 22.5% | +21.8 |
| 6 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 60.7% | 38.8% | +22.0 |
| 7 | KLF17(Zf)/2cell-Klf17-CutnTag(GSE211845)/Homer | 0.0000 | 64.5% | 42.6% | +21.9 |
| 8 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 50.6% | 29.6% | +21.1 |
| 9 | Ronin(THAP)/ES-Thap11-ChIP-Seq(GSE51522)/Homer | 0.0000 | 6.7% | 0.7% | +6.0 |
| 10 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 62.8% | 42.4% | +20.4 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Elk4(ETS)/Hela-Elk4-ChIP-Seq(GSE31477)/Homer | 0.0000 | 49.2% | 25.4% | +23.8 |
| 2 | Sp1(Zf)/Promoter/Homer | 0.0000 | 49.8% | 26.0% | +23.8 |
| 3 | Elk1(ETS)/Hela-Elk1-ChIP-Seq(GSE31477)/Homer | 0.0000 | 48.1% | 25.4% | +22.8 |
| 4 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 60.7% | 38.8% | +22.0 |
| 5 | KLF17(Zf)/2cell-Klf17-CutnTag(GSE211845)/Homer | 0.0000 | 64.5% | 42.6% | +21.9 |
| 6 | ELF1(ETS)/Jurkat-ELF1-ChIP-Seq(SRA014231)/Homer | 0.0000 | 44.4% | 22.5% | +21.8 |
| 7 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 50.6% | 29.6% | +21.1 |
| 8 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 62.8% | 42.4% | +20.4 |
| 9 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 50.6% | 31.9% | +18.7 |
| 10 | ETS(ETS)/Promoter/Homer | 0.0000 | 31.1% | 12.7% | +18.5 |

---

## out_conservedEnh_centered_vs_mouseSpec
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 2

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 27.3% | 15.9% | +11.4 |
| 2 | BORIS(Zf)/K562-CTCFL-ChIP-Seq(GSE32465)/Homer | 0.0000 | 28.4% | 17.3% | +11.1 |
| 3 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 15.8% | 12.1% | +3.7 |
| 4 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 16.5% | 12.7% | +3.8 |
| 5 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 15.5% | 11.9% | +3.6 |
| 6 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 16.9% | 13.1% | +3.8 |
| 7 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0000 | 15.3% | 11.8% | +3.5 |
| 8 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 17.7% | 13.9% | +3.7 |
| 9 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 14.6% | 11.2% | +3.4 |
| 10 | Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer | 0.0000 | 10.5% | 7.8% | +2.8 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 27.3% | 15.9% | +11.4 |
| 2 | BORIS(Zf)/K562-CTCFL-ChIP-Seq(GSE32465)/Homer | 0.0000 | 28.4% | 17.3% | +11.1 |

---

## out_conservedEnh_hg38_centered_vs_humanSpec
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 11

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 23.5% | 15.7% | +7.8 |
| 2 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 22.3% | 14.8% | +7.5 |
| 3 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 27.8% | 19.7% | +8.1 |
| 4 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 20.6% | 13.5% | +7.1 |
| 5 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 24.8% | 17.4% | +7.4 |
| 6 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 18.2% | 12.0% | +6.2 |
| 7 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 22.4% | 15.7% | +6.8 |
| 8 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 17.6% | 11.6% | +6.0 |
| 9 | BORIS(Zf)/K562-CTCFL-ChIP-Seq(GSE32465)/Homer | 0.0000 | 24.5% | 17.9% | +6.7 |
| 10 | EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Hom | 0.0000 | 12.7% | 8.1% | +4.7 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 27.8% | 19.7% | +8.1 |
| 2 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 23.5% | 15.7% | +7.8 |
| 3 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 22.3% | 14.8% | +7.5 |
| 4 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 24.8% | 17.4% | +7.4 |
| 5 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 20.6% | 13.5% | +7.1 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 22.4% | 15.7% | +6.8 |
| 7 | BORIS(Zf)/K562-CTCFL-ChIP-Seq(GSE32465)/Homer | 0.0000 | 24.5% | 17.9% | +6.7 |
| 8 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 18.2% | 12.0% | +6.2 |
| 9 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 17.6% | 11.6% | +6.0 |
| 10 | EHF(ETS)/LoVo-EHF-ChIP-Seq(GSE49402)/Homer | 0.0000 | 19.9% | 14.6% | +5.3 |

---

## out_conservedEnh_hg38_vs_humanSpec
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 167

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 63.2% | 47.2% | +16.0 |
| 2 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 45.6% | 30.8% | +14.8 |
| 3 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 49.1% | 34.1% | +15.0 |
| 4 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 51.5% | 36.4% | +15.1 |
| 5 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 57.4% | 42.2% | +15.2 |
| 6 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 43.4% | 29.3% | +14.1 |
| 7 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 51.9% | 37.2% | +14.6 |
| 8 | EHF(ETS)/LoVo-EHF-ChIP-Seq(GSE49402)/Homer | 0.0000 | 51.4% | 37.4% | +14.0 |
| 9 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 45.1% | 31.9% | +13.2 |
| 10 | EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Hom | 0.0000 | 34.6% | 22.6% | +12.0 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 63.2% | 47.2% | +16.0 |
| 2 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 57.4% | 42.2% | +15.2 |
| 3 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 51.5% | 36.4% | +15.1 |
| 4 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 49.1% | 34.1% | +15.0 |
| 5 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 45.6% | 30.8% | +14.8 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 51.9% | 37.2% | +14.6 |
| 7 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 43.4% | 29.3% | +14.1 |
| 8 | EHF(ETS)/LoVo-EHF-ChIP-Seq(GSE49402)/Homer | 0.0000 | 51.4% | 37.4% | +14.0 |
| 9 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 45.1% | 31.9% | +13.2 |
| 10 | SPDEF(ETS)/VCaP-SPDEF-ChIP-Seq(SRA014231)/Homer | 0.0000 | 42.3% | 30.1% | +12.2 |

---

## out_conservedEnh_vs_mouseSpec
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 13

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 30.4% | 17.7% | +12.7 |
| 2 | BORIS(Zf)/K562-CTCFL-ChIP-Seq(GSE32465)/Homer | 0.0000 | 33.3% | 20.9% | +12.3 |
| 3 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 25.0% | 18.3% | +6.8 |
| 4 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 25.7% | 18.9% | +6.8 |
| 5 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 23.4% | 17.0% | +6.4 |
| 6 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 27.5% | 20.6% | +6.8 |
| 7 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0000 | 24.6% | 18.3% | +6.3 |
| 8 | Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer | 0.0000 | 16.1% | 11.0% | +5.1 |
| 9 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 26.7% | 20.4% | +6.3 |
| 10 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 29.2% | 22.7% | +6.5 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 30.4% | 17.7% | +12.7 |
| 2 | BORIS(Zf)/K562-CTCFL-ChIP-Seq(GSE32465)/Homer | 0.0000 | 33.3% | 20.9% | +12.3 |
| 3 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 27.5% | 20.6% | +6.8 |
| 4 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 25.7% | 18.9% | +6.8 |
| 5 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 25.0% | 18.3% | +6.8 |
| 6 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 29.2% | 22.7% | +6.5 |
| 7 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 23.4% | 17.0% | +6.4 |
| 8 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0000 | 24.6% | 18.3% | +6.3 |
| 9 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 26.7% | 20.4% | +6.3 |
| 10 | GLIS3(Zf)/Thyroid-Glis3.GFP-ChIP-Seq(GSE103297)/Homer | 0.0000 | 49.4% | 43.6% | +5.8 |

---

## out_conservedProm_centered_vs_mouseSpec
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 1

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | NFY(CCAAT)/Promoter/Homer | 0.0000 | 23.8% | 18.3% | +5.5 |
| 2 | Ronin(THAP)/ES-Thap11-ChIP-Seq(GSE51522)/Homer | 0.0000 | 5.3% | 3.1% | +2.2 |
| 3 | GFY-Staf(?,Zf)/Promoter/Homer | 0.0000 | 5.6% | 3.5% | +2.1 |
| 4 | GFY(?)/Promoter/Homer | 0.0000 | 5.0% | 3.3% | +1.7 |
| 5 | En1(Homeobox)/SUM149-EN1-ChIP-Seq(GSE120957)/Homer | 0.0000 | 12.3% | 10.1% | +2.2 |
| 6 | REST-NRSF(Zf)/Jurkat-NRSF-ChIP-Seq/Homer | 0.0000 | 0.2% | 0.0% | +0.2 |
| 7 | Tgif1(Homeobox)/mES-Tgif1-ChIP-Seq(GSE55404)/Homer | 0.0000 | 17.3% | 14.8% | +2.4 |
| 8 | Tgif2(Homeobox)/mES-Tgif2-ChIP-Seq(GSE55404)/Homer | 0.0000 | 20.0% | 17.4% | +2.5 |
| 9 | Hoxd12(Homeobox)/ChickenMSG-Hoxd12.Flag-ChIP-Seq(GSE86088)/H | 0.0000 | 5.0% | 3.7% | +1.3 |
| 10 | YY1(Zf)/Promoter/Homer | 0.0000 | 4.3% | 3.1% | +1.2 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | NFY(CCAAT)/Promoter/Homer | 0.0000 | 23.8% | 18.3% | +5.5 |

---

## out_conservedProm_hg38_centered_vs_humanSpec
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 31

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | NFY(CCAAT)/Promoter/Homer | 0.0000 | 22.5% | 9.4% | +13.1 |
| 2 | Elk4(ETS)/Hela-Elk4-ChIP-Seq(GSE31477)/Homer | 0.0000 | 24.3% | 11.4% | +12.9 |
| 3 | ETS(ETS)/Promoter/Homer | 0.0000 | 15.6% | 5.7% | +9.9 |
| 4 | Elk1(ETS)/Hela-Elk1-ChIP-Seq(GSE31477)/Homer | 0.0000 | 24.3% | 11.8% | +12.5 |
| 5 | ELF1(ETS)/Jurkat-ELF1-ChIP-Seq(SRA014231)/Homer | 0.0000 | 22.2% | 10.3% | +11.9 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 31.2% | 17.1% | +14.1 |
| 7 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 23.4% | 11.6% | +11.8 |
| 8 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 28.4% | 15.6% | +12.8 |
| 9 | Klf15(Zf)/Liver-Klf15-ChIP-Seq(GSE166083)/Homer | 0.0000 | 53.1% | 37.1% | +16.0 |
| 10 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 21.1% | 10.2% | +10.9 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Klf15(Zf)/Liver-Klf15-ChIP-Seq(GSE166083)/Homer | 0.0000 | 53.1% | 37.1% | +16.0 |
| 2 | Sp5(Zf)/mES-Sp5.Flag-ChIP-Seq(GSE72989)/Homer | 0.0000 | 54.8% | 39.9% | +14.9 |
| 3 | KLF17(Zf)/2cell-Klf17-CutnTag(GSE211845)/Homer | 0.0000 | 46.5% | 31.6% | +14.9 |
| 4 | KLF1(Zf)/HUDEP2-KLF1-CutnRun(GSE136251)/Homer | 0.0000 | 48.1% | 33.9% | +14.2 |
| 5 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 31.2% | 17.1% | +14.1 |
| 6 | Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer | 0.0000 | 62.2% | 48.7% | +13.5 |
| 7 | Sp1(Zf)/Promoter/Homer | 0.0000 | 35.9% | 22.5% | +13.4 |
| 8 | NFY(CCAAT)/Promoter/Homer | 0.0000 | 22.5% | 9.4% | +13.1 |
| 9 | Elk4(ETS)/Hela-Elk4-ChIP-Seq(GSE31477)/Homer | 0.0000 | 24.3% | 11.4% | +12.9 |
| 10 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 28.4% | 15.6% | +12.8 |

---

## out_conservedProm_hg38_vs_humanSpec
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 269

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Elk4(ETS)/Hela-Elk4-ChIP-Seq(GSE31477)/Homer | 0.0000 | 57.5% | 31.8% | +25.8 |
| 2 | ELF1(ETS)/Jurkat-ELF1-ChIP-Seq(SRA014231)/Homer | 0.0000 | 52.5% | 27.7% | +24.8 |
| 3 | Elk1(ETS)/Hela-Elk1-ChIP-Seq(GSE31477)/Homer | 0.0000 | 56.9% | 31.8% | +25.1 |
| 4 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 70.2% | 45.2% | +25.1 |
| 5 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 72.2% | 47.2% | +25.0 |
| 6 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 58.4% | 33.6% | +24.8 |
| 7 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 70.4% | 46.0% | +24.5 |
| 8 | Sp1(Zf)/Promoter/Homer | 0.0000 | 60.7% | 36.6% | +24.1 |
| 9 | KLF17(Zf)/2cell-Klf17-CutnTag(GSE211845)/Homer | 0.0000 | 77.2% | 54.1% | +23.1 |
| 10 | Klf15(Zf)/Liver-Klf15-ChIP-Seq(GSE166083)/Homer | 0.0000 | 85.4% | 64.0% | +21.4 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Elk4(ETS)/Hela-Elk4-ChIP-Seq(GSE31477)/Homer | 0.0000 | 57.5% | 31.8% | +25.8 |
| 2 | Elk1(ETS)/Hela-Elk1-ChIP-Seq(GSE31477)/Homer | 0.0000 | 56.9% | 31.8% | +25.1 |
| 3 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 70.2% | 45.2% | +25.1 |
| 4 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 72.2% | 47.2% | +25.0 |
| 5 | ELF1(ETS)/Jurkat-ELF1-ChIP-Seq(SRA014231)/Homer | 0.0000 | 52.5% | 27.7% | +24.8 |
| 6 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 58.4% | 33.6% | +24.8 |
| 7 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 70.4% | 46.0% | +24.5 |
| 8 | Sp1(Zf)/Promoter/Homer | 0.0000 | 60.7% | 36.6% | +24.1 |
| 9 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 56.9% | 33.7% | +23.2 |
| 10 | KLF17(Zf)/2cell-Klf17-CutnTag(GSE211845)/Homer | 0.0000 | 77.2% | 54.1% | +23.1 |

---

## out_conservedProm_vs_mouseSpec
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 44

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Tgif2(Homeobox)/mES-Tgif2-ChIP-Seq(GSE55404)/Homer | 0.0000 | 63.3% | 54.1% | +9.2 |
| 2 | E2F3(E2F)/MEF-E2F3-ChIP-Seq(GSE71376)/Homer | 0.0000 | 56.7% | 47.5% | +9.2 |
| 3 | ZNF91(Zf)/HEK-ZNF91.HA-ChIP-Seq(GSE162571)/Homer | 0.0000 | 56.6% | 47.7% | +9.0 |
| 4 | LRF(Zf)/Erythroblasts-ZBTB7A-ChIP-Seq(GSE74977)/Homer | 0.0000 | 71.2% | 63.1% | +8.1 |
| 5 | Rfx6(HTH)/Min6b1-Rfx6.HA-ChIP-Seq(GSE62844)/Homer | 0.0000 | 43.6% | 35.5% | +8.2 |
| 6 | Zac1(Zf)/Neuro2A-Plagl1-ChIP-Seq(GSE75942)/Homer | 0.0000 | 79.6% | 72.7% | +6.9 |
| 7 | Tbx5(T-box)/HL1-Tbx5.biotin-ChIP-Seq(GSE21529)/Homer | 0.0000 | 70.6% | 63.0% | +7.6 |
| 8 | STAT4(Stat)/CD4-Stat4-ChIP-Seq(GSE22104)/Homer | 0.0000 | 22.7% | 16.5% | +6.2 |
| 9 | Tgif1(Homeobox)/mES-Tgif1-ChIP-Seq(GSE55404)/Homer | 0.0000 | 57.6% | 49.7% | +7.9 |
| 10 | Meis1(Homeobox)/MastCells-Meis1-ChIP-Seq(GSE48085)/Homer | 0.0000 | 47.3% | 39.6% | +7.7 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Tgif2(Homeobox)/mES-Tgif2-ChIP-Seq(GSE55404)/Homer | 0.0000 | 63.3% | 54.1% | +9.2 |
| 2 | E2F3(E2F)/MEF-E2F3-ChIP-Seq(GSE71376)/Homer | 0.0000 | 56.7% | 47.5% | +9.2 |
| 3 | ZNF91(Zf)/HEK-ZNF91.HA-ChIP-Seq(GSE162571)/Homer | 0.0000 | 56.6% | 47.7% | +9.0 |
| 4 | Rfx6(HTH)/Min6b1-Rfx6.HA-ChIP-Seq(GSE62844)/Homer | 0.0000 | 43.6% | 35.5% | +8.2 |
| 5 | LRF(Zf)/Erythroblasts-ZBTB7A-ChIP-Seq(GSE74977)/Homer | 0.0000 | 71.2% | 63.1% | +8.1 |
| 6 | Tgif1(Homeobox)/mES-Tgif1-ChIP-Seq(GSE55404)/Homer | 0.0000 | 57.6% | 49.7% | +7.9 |
| 7 | Meis1(Homeobox)/MastCells-Meis1-ChIP-Seq(GSE48085)/Homer | 0.0000 | 47.3% | 39.6% | +7.7 |
| 8 | Tbx5(T-box)/HL1-Tbx5.biotin-ChIP-Seq(GSE21529)/Homer | 0.0000 | 70.6% | 63.0% | +7.6 |
| 9 | E2F6(E2F)/Hela-E2F6-ChIP-Seq(GSE31477)/Homer | 0.0000 | 50.7% | 43.1% | +7.5 |
| 10 | Bapx1(Homeobox)/VertebralCol-Bapx1-ChIP-Seq(GSE36672)/Homer | 0.0000 | 51.9% | 44.4% | +7.5 |

---

## out_humanSpecEnh_hg38_centered_bidir
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 5

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer | 0.0000 | 23.1% | 14.3% | +8.8 |
| 2 | Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer | 0.0000 | 18.3% | 11.1% | +7.1 |
| 3 | SF1(NR)/H295R-Nr5a1-ChIP-Seq(GSE44220)/Homer | 0.0000 | 15.4% | 9.5% | +5.9 |
| 4 | Esrrb(NR)/mES-Esrrb-ChIP-Seq(GSE11431)/Homer | 0.0000 | 18.3% | 12.3% | +6.0 |
| 5 | ERRg(NR)/Kidney-ESRRG-ChIP-Seq(GSE104905)/Homer | 0.0000 | 20.9% | 15.3% | +5.6 |
| 6 | RAR:RXR(NR),DR0/ES-RAR-ChIP-Seq(GSE56893)/Homer | 0.0000 | 3.3% | 1.4% | +1.9 |
| 7 | REST-NRSF(Zf)/Jurkat-NRSF-ChIP-Seq/Homer | 0.0000 | 0.3% | 0.0% | +0.3 |
| 8 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 36.4% | 31.5% | +4.8 |
| 9 | RARg(NR)/ES-RARg-ChIP-Seq(GSE30538)/Homer | 0.0000 | 1.1% | 0.4% | +0.7 |
| 10 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 33.4% | 29.2% | +4.2 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer | 0.0000 | 23.1% | 14.3% | +8.8 |
| 2 | Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer | 0.0000 | 18.3% | 11.1% | +7.1 |
| 3 | Esrrb(NR)/mES-Esrrb-ChIP-Seq(GSE11431)/Homer | 0.0000 | 18.3% | 12.3% | +6.0 |
| 4 | SF1(NR)/H295R-Nr5a1-ChIP-Seq(GSE44220)/Homer | 0.0000 | 15.4% | 9.5% | +5.9 |
| 5 | ERRg(NR)/Kidney-ESRRG-ChIP-Seq(GSE104905)/Homer | 0.0000 | 20.9% | 15.3% | +5.6 |

---

## out_humanSpecEnh_hg38_vs_conserved
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 2

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer | 0.0000 | 38.4% | 32.0% | +6.4 |
| 2 | Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer | 0.0000 | 31.1% | 25.8% | +5.3 |
| 3 | SF1(NR)/H295R-Nr5a1-ChIP-Seq(GSE44220)/Homer | 0.0000 | 27.0% | 22.0% | +4.9 |
| 4 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 47.9% | 43.4% | +4.5 |
| 5 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 45.3% | 41.3% | +4.0 |
| 6 | Esrrb(NR)/mES-Esrrb-ChIP-Seq(GSE11431)/Homer | 0.0000 | 32.7% | 29.0% | +3.7 |
| 7 | RARg(NR)/ES-RARg-ChIP-Seq(GSE30538)/Homer | 0.0000 | 2.2% | 1.2% | +1.0 |
| 8 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 49.5% | 45.6% | +3.9 |
| 9 | ERRg(NR)/Kidney-ESRRG-ChIP-Seq(GSE104905)/Homer | 0.0000 | 37.5% | 33.9% | +3.6 |
| 10 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 46.9% | 43.1% | +3.7 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer | 0.0000 | 38.4% | 32.0% | +6.4 |
| 2 | Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer | 0.0000 | 31.1% | 25.8% | +5.3 |

---

## out_humanSpecProm_hg38_centered_bidir
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 12

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer | 0.0000 | 15.9% | 8.9% | +7.1 |
| 2 | Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer | 0.0000 | 12.5% | 6.8% | +5.7 |
| 3 | REST-NRSF(Zf)/Jurkat-NRSF-ChIP-Seq/Homer | 0.0000 | 1.2% | 0.1% | +1.1 |
| 4 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 24.6% | 17.5% | +7.1 |
| 5 | SF1(NR)/H295R-Nr5a1-ChIP-Seq(GSE44220)/Homer | 0.0000 | 10.8% | 6.1% | +4.7 |
| 6 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 22.1% | 15.8% | +6.3 |
| 7 | ERRg(NR)/Kidney-ESRRG-ChIP-Seq(GSE104905)/Homer | 0.0000 | 14.8% | 9.6% | +5.2 |
| 8 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 23.4% | 17.0% | +6.3 |
| 9 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 21.9% | 16.1% | +5.9 |
| 10 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 20.0% | 14.6% | +5.5 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0000 | 24.6% | 17.5% | +7.1 |
| 2 | Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer | 0.0000 | 15.9% | 8.9% | +7.1 |
| 3 | Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer | 0.0000 | 23.4% | 17.0% | +6.3 |
| 4 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0000 | 22.1% | 15.8% | +6.3 |
| 5 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0000 | 21.9% | 16.1% | +5.9 |
| 6 | Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer | 0.0000 | 12.5% | 6.8% | +5.7 |
| 7 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 21.0% | 15.4% | +5.6 |
| 8 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 20.0% | 14.6% | +5.5 |
| 9 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0000 | 20.6% | 15.1% | +5.5 |
| 10 | RARa(NR)/K562-RARa-ChIP-Seq(Encode)/Homer | 0.0000 | 32.0% | 26.6% | +5.4 |

---

## out_humanSpecProm_hg38_vs_conserved
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 0

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | REST-NRSF(Zf)/Jurkat-NRSF-ChIP-Seq/Homer | 0.0000 | 1.7% | 0.3% | +1.3 |
| 2 | Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer | 0.0000 | 28.2% | 25.1% | +3.1 |
| 3 | Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer | 0.0000 | 30.4% | 27.5% | +2.9 |
| 4 | Nrf2(bZIP)/Lymphoblast-Nrf2-ChIP-Seq(GSE37589)/Homer | 0.0002 | 3.5% | 2.6% | +0.9 |
| 5 | JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer | 0.0003 | 30.1% | 27.7% | +2.5 |
| 6 | Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer | 0.0007 | 31.9% | 29.6% | +2.4 |
| 7 | BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer | 0.0028 | 32.3% | 30.1% | +2.2 |
| 8 | Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer | 0.0047 | 17.6% | 15.9% | +1.7 |
| 9 | AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer | 0.0058 | 35.9% | 33.7% | +2.1 |
| 10 | Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer | 0.0062 | 22.3% | 20.5% | +1.8 |

No motifs passed both q<0.05 and delta>5.0.

---

## out_mouseSpecEnh_centered_bidir
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 13

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 56.5% | 46.0% | +10.5 |
| 2 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 41.9% | 31.9% | +9.9 |
| 3 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 42.5% | 32.9% | +9.7 |
| 4 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 44.0% | 34.3% | +9.7 |
| 5 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 50.5% | 40.6% | +9.9 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 43.4% | 33.9% | +9.5 |
| 7 | EHF(ETS)/LoVo-EHF-ChIP-Seq(GSE49402)/Homer | 0.0000 | 40.2% | 32.2% | +8.0 |
| 8 | EWS:FLI1-fusion(ETS)/SK_N_MC-EWS:FLI1-ChIP-Seq(SRA014231)/Ho | 0.0000 | 27.0% | 20.1% | +6.9 |
| 9 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 37.0% | 29.4% | +7.6 |
| 10 | EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Hom | 0.0000 | 29.6% | 22.7% | +6.9 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 56.5% | 46.0% | +10.5 |
| 2 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 41.9% | 31.9% | +9.9 |
| 3 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 50.5% | 40.6% | +9.9 |
| 4 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 44.0% | 34.3% | +9.7 |
| 5 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 42.5% | 32.9% | +9.7 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 43.4% | 33.9% | +9.5 |
| 7 | EHF(ETS)/LoVo-EHF-ChIP-Seq(GSE49402)/Homer | 0.0000 | 40.2% | 32.2% | +8.0 |
| 8 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 37.0% | 29.4% | +7.6 |
| 9 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 37.4% | 30.4% | +7.0 |
| 10 | EWS:FLI1-fusion(ETS)/SK_N_MC-EWS:FLI1-ChIP-Seq(SRA014231)/Ho | 0.0000 | 27.0% | 20.1% | +6.9 |

---

## out_mouseSpecEnh_vs_conserved_bidir
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 14

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 54.8% | 45.9% | +8.9 |
| 2 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 56.3% | 47.6% | +8.8 |
| 3 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 70.5% | 62.1% | +8.4 |
| 4 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 57.5% | 49.2% | +8.3 |
| 5 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 65.0% | 57.1% | +7.9 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 56.7% | 48.9% | +7.8 |
| 7 | EHF(ETS)/LoVo-EHF-ChIP-Seq(GSE49402)/Homer | 0.0000 | 56.4% | 49.2% | +7.2 |
| 8 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 52.1% | 45.1% | +7.0 |
| 9 | EWS:FLI1-fusion(ETS)/SK_N_MC-EWS:FLI1-ChIP-Seq(SRA014231)/Ho | 0.0000 | 37.6% | 31.6% | +6.1 |
| 10 | EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Hom | 0.0000 | 42.0% | 36.0% | +6.0 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 54.8% | 45.9% | +8.9 |
| 2 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 56.3% | 47.6% | +8.8 |
| 3 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 70.5% | 62.1% | +8.4 |
| 4 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 57.5% | 49.2% | +8.3 |
| 5 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 65.0% | 57.1% | +7.9 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 56.7% | 48.9% | +7.8 |
| 7 | EHF(ETS)/LoVo-EHF-ChIP-Seq(GSE49402)/Homer | 0.0000 | 56.4% | 49.2% | +7.2 |
| 8 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 52.1% | 45.1% | +7.0 |
| 9 | EWS:FLI1-fusion(ETS)/SK_N_MC-EWS:FLI1-ChIP-Seq(SRA014231)/Ho | 0.0000 | 37.6% | 31.6% | +6.1 |
| 10 | EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Hom | 0.0000 | 42.0% | 36.0% | +6.0 |

---

## out_mouseSpecProm_centered_bidir
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 7

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 41.6% | 34.8% | +6.9 |
| 2 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 29.9% | 24.0% | +5.9 |
| 3 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 39.9% | 33.6% | +6.2 |
| 4 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 41.9% | 35.7% | +6.2 |
| 5 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 34.1% | 28.3% | +5.9 |
| 6 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 40.3% | 34.3% | +6.1 |
| 7 | Ets1-distal(ETS)/CD4+-PolII-ChIP-Seq(Barski_et_al.)/Homer | 0.0000 | 10.4% | 7.1% | +3.3 |
| 8 | EWS:ERG-fusion(ETS)/CADO_ES1-EWS:ERG-ChIP-Seq(SRA014231)/Hom | 0.0000 | 16.0% | 12.0% | +4.1 |
| 9 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 32.8% | 27.7% | +5.1 |
| 10 | Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer | 0.0000 | 31.4% | 26.9% | +4.5 |

### Top 10 by delta_pct (q<0.05)
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.0000 | 41.6% | 34.8% | +6.9 |
| 2 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 0.0000 | 39.9% | 33.6% | +6.2 |
| 3 | ETV1(ETS)/GIST48-ETV1-ChIP-Seq(GSE22441)/Homer | 0.0000 | 41.9% | 35.7% | +6.2 |
| 4 | ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer | 0.0000 | 40.3% | 34.3% | +6.1 |
| 5 | Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer | 0.0000 | 29.9% | 24.0% | +5.9 |
| 6 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 0.0000 | 34.1% | 28.3% | +5.9 |
| 7 | GABPA(ETS)/Jurkat-GABPa-ChIP-Seq(GSE17954)/Homer | 0.0000 | 32.8% | 27.7% | +5.1 |

---

## out_mouseSpecProm_vs_conserved_bidir
Total motifs tested: 472
Significant (q<0.05 & delta>5.0): 0

### Top 10 by q-value
| Rank | Motif | q-value | %Target | %Bg | delta |
|------|-------|---------|---------|-----|-------|
| 1 | T1ISRE(IRF)/ThioMac-Ifnb-Expression/Homer | 0.0859 | 0.7% | 0.3% | +0.4 |
| 2 | Ets1-distal(ETS)/CD4+-PolII-ChIP-Seq(Barski_et_al.)/Homer | 0.0859 | 16.5% | 14.4% | +2.1 |
| 3 | EWS:FLI1-fusion(ETS)/SK_N_MC-EWS:FLI1-ChIP-Seq(SRA014231)/Ho | 0.1647 | 33.7% | 31.2% | +2.4 |
| 4 | ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer | 0.2445 | 62.9% | 60.5% | +2.4 |
| 5 | Pitx1:Ebox(Homeobox,bHLH)/Hindlimb-Pitx1-ChIP-Seq(GSE41591)/ | 0.4188 | 2.4% | 1.8% | +0.6 |
| 6 | Oct4:Sox17(POU,Homeobox,HMG)/F9-Sox17-ChIP-Seq(GSE44553)/Hom | 0.9837 | 1.9% | 1.4% | +0.5 |
| 7 | PGR(NR)/EndoStromal-PGR-ChIP-Seq(GSE69539)/Homer | 1.0000 | 5.2% | 4.5% | +0.8 |
| 8 | ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer | 1.0000 | 52.3% | 50.7% | +1.6 |
| 9 | Fli1(ETS)/CD8-FLI-ChIP-Seq(GSE20898)/Homer | 1.0000 | 59.5% | 57.9% | +1.6 |
| 10 | GATA(Zf),IR4/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer | 1.0000 | 1.0% | 0.7% | +0.3 |

No motifs passed both q<0.05 and delta>5.0.

---

