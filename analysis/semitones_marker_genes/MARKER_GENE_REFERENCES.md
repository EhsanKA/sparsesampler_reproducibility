# Marker Gene References for SEMITONES Analysis

This document contains citable references for marker genes used as ground truth in the SEMITONES marker gene analysis comparing SPS vs Random sampling methods.

---

## Primary Marker Gene Databases

### 1. CellMarker 2.0 (Recommended - Most Comprehensive)

**Citation:**
> Hu C., Li T., Xu Y., Zhang X., Li F., Bai J., Chen J., Jiang W., Yang K., Ou Q., Li X., Wang P., Zhang Y. (2023). CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data. *Nucleic Acids Research*, 51(D1), D870–D876. https://doi.org/10.1093/nar/gkac947

**Key Features:**
- 83,361 tissue–cell type–marker entries
- 656 tissues and 2,578 cell types
- 26,915 unique marker genes
- Incorporates marker gene information from 48 single-cell sequencing technologies
- URL: http://bio-bigdata.hrbmu.edu.cn/CellMarker/

### 2. PanglaoDB

**Citation:**
> Franzén O., Gan L.M., Björkegren J.L.M. (2019). PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data. *Database*, 2019, baz046. https://doi.org/10.1093/database/baz046

**Key Features:**
- ~6,600 marker gene-to-cell-type associations
- 155 distinct cell types across human and mouse
- Manually curated from thousands of published articles
- Includes sensitivity and specificity scores for markers
- URL: https://panglaodb.se/

---

## Marker Genes by Cell Type

### Osteoblast

**Markers:** RUNX2, SP7 (Osterix), ALPL, COL1A1, BGLAP (Osteocalcin), SPARC, BSP, OPN (SPP1)

**References:**
- CellMarker 2.0 (Hu et al., 2023)
- PanglaoDB (Franzén et al., 2019)
- PMC9583853 - Bone cell markers review

**Marker Descriptions:**
- RUNX2: Master transcription factor for osteoblast differentiation
- SP7: Essential for bone formation, follows RUNX2
- ALPL: Early osteoblast function marker (Alkaline Phosphatase)
- COL1A1: Major matrix protein synthesized by osteoblasts
- BGLAP: Late marker, secreted by mature osteoblasts

---

### Chondrocyte

**Markers:** SOX9, SOX5, SOX6, COL2A1, ACAN, COL10A1, COL11A1

**References:**
- CellMarker 2.0 (Hu et al., 2023)
- PanglaoDB (Franzén et al., 2019)
- Ji Q., et al. (2019). Single-cell RNA-seq analysis reveals the progression of human osteoarthritis. *Annals of the Rheumatic Diseases*, 78(1), 100-110. (PMC7724342)

**Marker Descriptions:**
- SOX9: Key transcription factor for chondrogenic lineage
- COL2A1: Collagen type II, cartilage matrix component
- ACAN: Aggrecan, major proteoglycan in cartilage
- COL10A1: Marker of hypertrophic chondrocytes

---

### Fibroblast

**Markers:** VIM, PDGFRA, PDGFRB, COL1A1, COL1A2, THY1, FAP, S100A4, ACTA2

**References:**
- CellMarker 2.0 (Hu et al., 2023)
- PanglaoDB (Franzén et al., 2019)

**Marker Descriptions:**
- VIM (Vimentin): Cytoskeletal marker of mesenchymal origin
- PDGFRA/PDGFRB: Growth factor receptors in fibroblasts
- COL1A1/COL1A2: Collagen type I production
- THY1: Surface marker
- FAP: Fibroblast activation protein

---

### Mesodermal Cell

**Markers:** TBXT (Brachyury/T), MIXL1, MESP1, EOMES, KDR, GSC

**References:**
- CellMarker 2.0 (Hu et al., 2023)
- PMID:19134196 - Mesoderm specification markers

**Marker Descriptions:**
- TBXT: Classic mesoderm specification marker in early development
- MIXL1: Early mesendoderm marker
- MESP1: Cardiac mesoderm specification
- KDR (VEGFR2): Marks mesodermal precursors

---

### Lateral Mesodermal Cell

**Markers:** HAND1, FOXF1, GATA4, BMP4, WT1, TBX18

**References:**
- Loh K.M., et al. (2016). Mapping the Pairwise Choices Leading from Pluripotency to Human Bone, Heart, and Other Mesoderm Cell Types. *Cell*, 166(2), 451-467. https://doi.org/10.1016/j.cell.2016.06.011 (PMC5474394)
- Reactome Pathway: R-HSA-9758920

**Marker Descriptions:**
- HAND1: Transcription factor strongly expressed in LPM vs paraxial mesoderm (~98-100% of LPM cells)
- FOXF1: Central marker of LPM, required for lateral mesoderm development (~98-100% of LPM cells)
- GATA4: Acts downstream of BMP4 and FOXF1
- BMP4: Ligand in regulatory loop with FOXF1

**Note:** Negative markers to distinguish from other mesoderm: MSGN1, DLL3 (paraxial mesoderm), PAX2, PAX8 (intermediate mesoderm)

---

## Summary Table

| Cell Type | Key Markers | Primary Database |
|-----------|-------------|------------------|
| Osteoblast | RUNX2, SP7, ALPL, COL1A1, BGLAP, SPARC | CellMarker 2.0 |
| Chondrocyte | SOX9, COL2A1, ACAN, COL10A1, SOX5, SOX6 | CellMarker 2.0 |
| Fibroblast | VIM, PDGFRA, PDGFRB, COL1A1, COL1A2, THY1, FAP | CellMarker 2.0, PanglaoDB |
| Mesodermal cell | TBXT, MIXL1, MESP1, KDR, EOMES | CellMarker 2.0 |
| Lateral mesodermal cell | HAND1, FOXF1, GATA4, BMP4 | PMC5474394 |

---

## BibTeX References

```bibtex
@article{hu2023cellmarker,
  title={CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data},
  author={Hu, Congxue and Li, Tengyue and Xu, Yingqi and Zhang, Xinxin and Li, Feng and Bai, Jing and Chen, Jing and Jiang, Wenjin and Yang, Kaiyue and Ou, Qi and Li, Xia and Wang, Peng and Zhang, Yunpeng},
  journal={Nucleic Acids Research},
  volume={51},
  number={D1},
  pages={D870--D876},
  year={2023},
  doi={10.1093/nar/gkac947}
}

@article{franzen2019panglaodb,
  title={PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data},
  author={Franz{\'e}n, Oscar and Gan, Li-Ming and Bj{\"o}rkegren, Johan LM},
  journal={Database},
  volume={2019},
  pages={baz046},
  year={2019},
  doi={10.1093/database/baz046}
}

@article{loh2016mapping,
  title={Mapping the pairwise choices leading from pluripotency to human bone, heart, and other mesoderm cell types},
  author={Loh, Kyle M and Chen, Angela and Koh, Lay Teng and Deng, Ting and Sinha, Rajan and Tsai, Jonathan M and Barber, Avis A and Cai, Winnie and Somers, Amira and Delaney, John P and others},
  journal={Cell},
  volume={166},
  number={2},
  pages={451--467},
  year={2016},
  doi={10.1016/j.cell.2016.06.011}
}
```

---

## Notes for Analysis

1. **Ground Truth Definition:** Literature-based marker genes from CellMarker 2.0 will serve as ground truth for comparing SPS vs Random sampling methods.

2. **Evaluation Metrics:**
   - Recall: % of literature markers recovered in top N enriched genes
   - Rank: Average rank of literature markers in SEMITONES enrichment results
   - Enrichment score: Mean enrichment of literature markers

3. **Expected Outcome:** SPS sampling should show better recovery of literature markers, especially for rare cell types (osteoblast ~1% of dataset).

---

*Last updated: January 2026*
