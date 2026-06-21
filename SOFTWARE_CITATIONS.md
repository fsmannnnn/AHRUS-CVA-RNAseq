# Software & database citations

Versions used in this study, with formal references and DOIs. Tools without a
peer-reviewed DOI are cited as software (URL + version). The first block lists
upstream tools (run by the sequencing provider); the second lists the
downstream R/Bioconductor packages used in this repository.

---

## Upstream — bulk RNA-seq

1. **Skewer v0.2.2** — Jiang H, Lei R, Ding S-W, Zhu S. Skewer: a fast and
   accurate adapter trimmer for next-generation sequencing paired-end reads.
   *BMC Bioinformatics.* 2014;15:182. https://doi.org/10.1186/1471-2105-15-182

2. **FastQC v0.11.5** — Andrews S. *FastQC: a quality control tool for high
   throughput sequence data.* Babraham Bioinformatics, 2010. (Software, no DOI.)
   https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

3. **STAR v2.5.3** — Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast
   universal RNA-seq aligner. *Bioinformatics.* 2013;29(1):15–21.
   https://doi.org/10.1093/bioinformatics/bts635

4. **StringTie v1.3.1** — Pertea M, Pertea GM, Antonescu CM, et al. StringTie
   enables improved reconstruction of a transcriptome from RNA-seq reads.
   *Nat Biotechnol.* 2015;33(3):290–295. https://doi.org/10.1038/nbt.3122

---

## Upstream — 16S rRNA amplicon

5. **fastp v0.19.6** — Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast
   all-in-one FASTQ preprocessor. *Bioinformatics.* 2018;34(17):i884–i890.
   https://doi.org/10.1093/bioinformatics/bty560

6. **FLASH v1.2.11** — Magoč T, Salzberg SL. FLASH: fast length adjustment of
   short reads to improve genome assemblies. *Bioinformatics.*
   2011;27(21):2957–2963. https://doi.org/10.1093/bioinformatics/btr507

7. **QIIME 2 v2020.2** — Bolyen E, Rideout JR, Dillon MR, et al. Reproducible,
   interactive, scalable and extensible microbiome data science using QIIME 2.
   *Nat Biotechnol.* 2019;37(8):852–857.
   https://doi.org/10.1038/s41587-019-0209-9

8. **DADA2** — Callahan BJ, McMurdie PJ, Rosen MJ, et al. DADA2: high-resolution
   sample inference from Illumina amplicon data. *Nat Methods.*
   2016;13(7):581–583. https://doi.org/10.1038/nmeth.3869

9. **SILVA v138** — Quast C, Pruesse E, Yilmaz P, et al. The SILVA ribosomal RNA
   gene database project: improved data processing and web-based tools.
   *Nucleic Acids Res.* 2013;41(D1):D590–D596.
   https://doi.org/10.1093/nar/gks1219

---

## Downstream — R / Bioconductor

10. **DESeq2 v1.44.0** — Love MI, Huber W, Anders S. Moderated estimation of fold
    change and dispersion for RNA-seq data with DESeq2. *Genome Biol.*
    2014;15(12):550. https://doi.org/10.1186/s13059-014-0550-8

11. **edgeR (filterByExpr)** — Robinson MD, McCarthy DJ, Smyth GK. edgeR: a
    Bioconductor package for differential expression analysis of digital gene
    expression data. *Bioinformatics.* 2010;26(1):139–140.
    https://doi.org/10.1093/bioinformatics/btp616
    — and: Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways…
    *F1000Research.* 2016;5:1438. https://doi.org/10.12688/f1000research.8987.2

12. **IOBR** — Zeng D, Ye Z, Shen R, et al. IOBR: multi-omics immuno-oncology
    biological research to decode tumor microenvironment and signatures.
    *Front Immunol.* 2021;12:687975. https://doi.org/10.3389/fimmu.2021.687975

13. **CIBERSORT** — Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell
    subsets from tissue expression profiles. *Nat Methods.* 2015;12(5):453–457.
    https://doi.org/10.1038/nmeth.3337

14. **xCell** — Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue
    cellular heterogeneity landscape. *Genome Biol.* 2017;18(1):220.
    https://doi.org/10.1186/s13059-017-1349-1

15. **clusterProfiler** — Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: a
    universal enrichment tool for interpreting omics data. *The Innovation.*
    2021;2(3):100141. https://doi.org/10.1016/j.xinn.2021.100141
    — original: Yu G, Wang L-G, Han Y, He Q-Y. clusterProfiler: an R package for
    comparing biological themes among gene clusters. *OMICS.*
    2012;16(5):284–287. https://doi.org/10.1089/omi.2011.0118

16. **GSEA** — Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment
    analysis: a knowledge-based approach for interpreting genome-wide expression
    profiles. *Proc Natl Acad Sci USA.* 2005;102(43):15545–15550.
    https://doi.org/10.1073/pnas.0506580102

17. **MSigDB** — Liberzon A, Subramanian A, Pinchback R, et al. Molecular
    Signatures Database (MSigDB) 3.0. *Bioinformatics.* 2011;27(12):1739–1740.
    https://doi.org/10.1093/bioinformatics/btr260
    — Hallmark collection: Liberzon A, Birger C, Thorvaldsdóttir H, et al. The
    Molecular Signatures Database hallmark gene set collection. *Cell Syst.*
    2015;1(6):417–425. https://doi.org/10.1016/j.cels.2015.12.004

18. **msigdbr** — Dolgalev I. *msigdbr: MSigDB gene sets for multiple organisms
    in a tidy data format.* R package. (Software, no DOI.)
    https://cran.r-project.org/package=msigdbr

19. **fgsea** — Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment
    analysis. *bioRxiv.* 2021. https://doi.org/10.1101/060012

20. **EnhancedVolcano** — Blighe K, Rana S, Lewis M. *EnhancedVolcano:
    publication-ready volcano plots with enhanced colouring and labeling.*
    R/Bioconductor package. (Software, no DOI.)
    https://bioconductor.org/packages/EnhancedVolcano

21. **ComplexHeatmap** — Gu Z, Eils R, Schlesner M. Complex heatmaps reveal
    patterns and correlations in multidimensional genomic data. *Bioinformatics.*
    2016;32(18):2847–2849. https://doi.org/10.1093/bioinformatics/btw313

22. **pheatmap** — Kolde R. *pheatmap: pretty heatmaps.* R package.
    (Software, no DOI.) https://cran.r-project.org/package=pheatmap

23. **ggplot2** — Wickham H. *ggplot2: Elegant Graphics for Data Analysis.*
    Springer-Verlag New York, 2016. https://doi.org/10.1007/978-3-319-24277-4

24. **factoextra** — Kassambara A, Mundt F. *factoextra: extract and visualize
    the results of multivariate data analyses.* R package. (Software, no DOI.)
    https://cran.r-project.org/package=factoextra

25. **Seurat** — Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of
    multimodal single-cell data. *Cell.* 2021;184(13):3573–3587.e29.
    https://doi.org/10.1016/j.cell.2021.04.048

26. **CIBERSORTx** — Newman AM, Steen CB, Liu CL, et al. Determining cell type
    abundance and expression from bulk tissues with digital cytometry.
    *Nat Biotechnol.* 2019;37(7):773–782.
    https://doi.org/10.1038/s41587-019-0114-2

---

### Notes
- Confirm the **exact installed versions** of IOBR, clusterProfiler, Seurat and
  msigdbr from `results/sessionInfo.txt`; the numbers above reflect the versions
  used during analysis and should match your environment.
- `org.Mm.eg.db` is an annotation package (Bioconductor); cite its version from
  `sessionInfo()` if requested by the journal.
