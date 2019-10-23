Data analysis of single cell RNA-seq of thymus
================
Jihed Chouaref
22nd october 2019

The Goal of this notebook is to recapitulate all the data I have
gathered so far on the Thymus single-cell RNA-seq performed on the
homozygous sample.

# Dataset

The experiment was performed by Veronica and Kelly, on E18.5 days mouse
embryos. Only CD45+ cells were isolated and then prepare for single cell
RNA library preparation with 10X 3’.

Load the require packages for the analysis:

``` r
suppressMessages(require(Seurat))
suppressMessages(require(ggplot2))
suppressMessages(require(dplyr))
suppressMessages(require(scater))
suppressMessages(require(scran))
suppressMessages(require(Matrix))
```

First I import the dataset and create the seurat object.

``` r
sample1.data <- Read10X(data.dir="/Volumes/Elements/scRNA-seq Thymus/sample_1/filtered_feature_bc_matrix/")
sample1.seurat <- CreateSeuratObject(counts = sample1.data, project = "Thymus_homo_5582", min.cells = 3, min.features = 200)
rm(sample1.data) #save some memory 
```

## 1\. Quality control

### Calculate mitochondrial proportion

``` r
sample1.seurat[["percent.mt"]] <- PercentageFeatureSet(sample1.seurat, pattern = "^mt-")

head(sample1.seurat@meta.data)
```

    ##                        orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCCAAGAGTATAC Thymus_homo_5582      21685         3056  0.3965875
    ## AAACCCAAGATTGGGC Thymus_homo_5582      17267         3728  3.0404818
    ## AAACCCAAGGTAGTAT Thymus_homo_5582      14334         3471  1.9394447
    ## AAACCCACAATTGAAG Thymus_homo_5582      21514         4257  2.6540857
    ## AAACCCACACACACGC Thymus_homo_5582       3182         1048  0.1571339
    ## AAACCCACAGAGATTA Thymus_homo_5582       1192          616 10.4865772

``` r
VlnPlot(sample1.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol=3),pt.size = 0.1,cols = "red")
```

    ## Warning in FetchData(object = object, vars = features, slot = slot): The
    ## following requested variables were not found: 3

![](sample1.complete.analysis_files/figure-gfm/QC-mitochondria-plot-1.png)<!-- -->

### Calculate ribosomal proportion

``` r
sample1.seurat[["percent.rp"]]<-PercentageFeatureSet(sample1.seurat, pattern = "Rp[sl][[:digit:]]")
head(sample1.seurat@meta.data)
```

    ##                        orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCCAAGAGTATAC Thymus_homo_5582      21685         3056  0.3965875
    ## AAACCCAAGATTGGGC Thymus_homo_5582      17267         3728  3.0404818
    ## AAACCCAAGGTAGTAT Thymus_homo_5582      14334         3471  1.9394447
    ## AAACCCACAATTGAAG Thymus_homo_5582      21514         4257  2.6540857
    ## AAACCCACACACACGC Thymus_homo_5582       3182         1048  0.1571339
    ## AAACCCACAGAGATTA Thymus_homo_5582       1192          616 10.4865772
    ##                  percent.rp
    ## AAACCCAAGAGTATAC   1.918377
    ## AAACCCAAGATTGGGC  26.680952
    ## AAACCCAAGGTAGTAT  22.185015
    ## AAACCCACAATTGAAG  24.500325
    ## AAACCCACACACACGC   2.231301
    ## AAACCCACAGAGATTA   7.382550

``` r
a<- sample1.seurat[["percent.rp"]]
mean(a$percent.rp)
```

    ## [1] 14.20766

Cells contains on average 20% ribosomal RNA.

#### Plot QC

I plot all the QC as violin plots

``` r
VlnPlot(sample1.seurat, features = "nFeature_RNA", pt.size = 0.1, cols = "red") + NoLegend()
```

![](sample1.complete.analysis_files/figure-gfm/plot_nGenes-1.png)<!-- -->

``` r
VlnPlot(sample1.seurat, features = "nCount_RNA", pt.size = 0.1,  cols = "red") + NoLegend()
```

![](sample1.complete.analysis_files/figure-gfm/plot_nUMIs-1.png)<!-- -->

``` r
VlnPlot(sample1.seurat, features = "percent.mt", pt.size = 0.1,  cols = "red", y.max = 20) + NoLegend()
```

    ## Warning: Removed 29 rows containing non-finite values (stat_ydensity).

    ## Warning: Removed 29 rows containing missing values (geom_point).

![](sample1.complete.analysis_files/figure-gfm/plot_mito-1.png)<!-- -->

``` r
VlnPlot(sample1.seurat, features = "percent.rp", pt.size = 0.1,  cols = "red") + NoLegend()
```

![](sample1.complete.analysis_files/figure-gfm/plot_ribo-1.png)<!-- -->

To help me to choose the values for the filtering I am plotting the QC
information as scatter plots

``` r
FeatureScatter(sample1.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "red")
```

![](sample1.complete.analysis_files/figure-gfm/plot_scatter1-1.png)<!-- -->

``` r
FeatureScatter(sample1.seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt", cols = "red")
```

![](sample1.complete.analysis_files/figure-gfm/plot_scatter2-1.png)<!-- -->

``` r
FeatureScatter(sample1.seurat, feature1="nFeature_RNA", feature2="percent.rp", cols = "red")
```

![](sample1.complete.analysis_files/figure-gfm/plot_scatter3-1.png)<!-- -->
This plot seems to indicate that there is a correlation between the
quantity of ribosomal RNA and the total quantity of RNA detected in a
cell. This makes sense but is it due to biology that some cells have
more ribosomal RNA? Are they cycling?

``` r
plots<- list(FeatureScatter(sample1.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "red"))
plots[[2]]<- FeatureScatter(sample1.seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt", cols = "red")
plots[[3]]<- FeatureScatter(sample1.seurat, feature1="nFeature_RNA", feature2="percent.rp", cols = "red")
CombinePlots(plots = plots, ncol = 2, label_size = 12)
```

![](sample1.complete.analysis_files/figure-gfm/plot_scatterall-1.png)<!-- -->

The violin plots indicates that most cells contains less than 5000
transcript per cells. I have then decided to plot select the cells with
a minimum of 200 and a maximum of 5000 transcripts per cells, and a
percentage of mt RNA under 5%. I have decided to not filter the cells
based on the percentage of ribosomal RNA because it would imply
filtering cells that are in a different step of the cell cycle and could
be a biological effect.

An interesting observation is that both violin plot of mt-RNA and rp-RNA
displays 2 groups of populations.

### 2.Filtering

#### Before filtering

``` r
sample1.seurat
```

    ## An object of class Seurat 
    ## 15365 features across 5582 samples within 1 assay 
    ## Active assay: RNA (15365 features)

``` r
sample1.seurat<-subset(sample1.seurat,subset= nFeature_RNA> 200 & nFeature_RNA<5000 & percent.mt <5 )
sample1.seurat
```

    ## An object of class Seurat 
    ## 15365 features across 5122 samples within 1 assay 
    ## Active assay: RNA (15365 features)

With those thresholds, 453 cells were filtered.

#### Plot QC-stats again

Lets plot the same qc-stats another time.

``` r
VlnPlot(sample1.seurat, features = "nFeature_RNA", pt.size = 0.1, cols = "red") + NoLegend()
```

![](sample1.complete.analysis_files/figure-gfm/vln.plot2-1.png)<!-- -->

``` r
VlnPlot(sample1.seurat, features = "nCount_RNA", pt.size = 0.1, cols = "red") + NoLegend()
```

![](sample1.complete.analysis_files/figure-gfm/vln.plot2-2.png)<!-- -->

``` r
VlnPlot(sample1.seurat, features = "percent.mt", pt.size = 0.1, cols = "red") + NoLegend()
```

![](sample1.complete.analysis_files/figure-gfm/vln.plot2-3.png)<!-- -->

``` r
VlnPlot(sample1.seurat, features = "percent.rp", pt.size = 0.1, cols = "red") + NoLegend()
```

![](sample1.complete.analysis_files/figure-gfm/vln.plot2-4.png)<!-- -->

``` r
# and check the number of cells per sample before and after filtering
table(Idents(sample1.seurat))
```

    ## 
    ## Thymus_homo_5582 
    ##             5122

### 3\. Normalization: Log

The normalization is done to compensate for the variability of the
number of RNA per cell. This normalizes the feature expression
measurements for each cell by the total expression, multiplies this by a
scale factor (10,000 by default), and log-transforms the result.
Normalized values are stored in

``` r
sample1.seurat <- NormalizeData(sample1.seurat, normalize.method = "LogNormalize", scale.factor = 10000)
```

    ## Warning: The following arguments are not used: normalize.method

The information about the normalization are stored in the
`sample1.seurat@assays$RNA@data`.

### 4\. Feature Selection

The default method in Seurat 3 is variance-stabilizing transformation. A
trend is fitted to to predict the variance of each gene as a function of
its mean. For each gene, the variance of standardized values is computed
across all cells and used to rank the features. By default, 2000 top
genes are returned.

``` r
sample1.seurat <- FindVariableFeatures(sample1.seurat, selection.method = "vst")
top10 <- head(VariableFeatures(sample1.seurat), 10)
vplot <- VariableFeaturePlot(sample1.seurat)
LabelPoints(plot = vplot, points = top10, repel = TRUE)
```

    ## Warning: Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
    ## Please use `as_label()` or `as_name()` instead.
    ## This warning is displayed once per session.

    ## When using repel, set xnudge and ynudge to 0 for optimal results

    ## Warning: Transformation introduced infinite values in continuous x-axis

![](sample1.complete.analysis_files/figure-gfm/hvg%20identification-1.png)<!-- -->

Those genes are the genes that shows the most variation in the dataset.
They will be used for the calculation of the principal components.

### 5\. Dimensionality Reduction

#### Scaling the data

Next, we apply a linear transformation (‘scaling’) that is a standard
pre-processing step prior to dimensional reduction techniques like PCA.
The ScaleData function:

  - Shifts the expression of each gene, so that the mean expression
    across cells is 0
  - Scales the expression of each gene, so that the variance across
    cells is 1 This step gives equal weight in downstream analyses, so
    that highly-expressed genes do not dominate The results of this are
    stored in `sample1.seurat[["RNA"]]@scale.data`

<!-- end list -->

``` r
all.genes <- rownames(sample1.seurat)
sample1.seurat <- ScaleData(sample1.seurat, features = all.genes)
```

    ## Centering and scaling data matrix

#### Principal component analysis

``` r
sample1.seurat <- RunPCA(sample1.seurat, features = VariableFeatures(object = sample1.seurat))
```

    ## PC_ 1 
    ## Positive:  Gm5483, Hp, Stfa2l1, BC100530, Spi1, Alox5ap, Wfdc21, Stfa2, Pglyrp1, Lyz2 
    ##     Anxa1, Gm5416, Mcemp1, Prdx5, Asprv1, Mmp8, Olfm4, F630028O10Rik, Chil1, Cebpd 
    ##     Retnlg, Slfn1, Stfa3, Raet1e, Mgst1, Cd3d, Mmp9, Cxcr2, Igsf6, Mxd1 
    ## Negative:  Ptma, Ran, Npm1, Rpl3, Hint1, Rpl12, Rpl14, Nme1, Rpl4, Set 
    ##     Snrpd1, Erh, Tmsb10, Eef1d, Rps17, Ranbp1, Lgals1, Snrpg, Rps2, Ptprcap 
    ##     Rpl5, Rps6, Serbp1, Selenoh, Rpl36a, Snrpf, Prdx1, mt-Nd1, Srsf3, Stmn1 
    ## PC_ 2 
    ## Positive:  Cpa3, Cyp11a1, Gata2, Slc7a8, Kit, Gm11697, Cma1, Hdc, Gstp1, Ms4a2 
    ##     Oaf, Hs3st1, Tph1, Mcpt4, Clec12b, Tpsb2, Rapsn, Chst1, Homer2, Fam105a 
    ##     Slc45a3, Smarca1, Lxn, Smpx, Rgs10, Mcpt2, Rab4a, Hpgds, Cd200r3, Car8 
    ## Negative:  Nusap1, Prc1, Knl1, Kif11, Mki67, Ccna2, Cenpf, Ube2c, Cdk1, Top2a 
    ##     Cdca8, Tpx2, Cdca3, H2afx, Birc5, Aurkb, Hist1h3c, Hist1h1b, Hmmr, Plk1 
    ##     Hist1h2ap, Sgo2a, Spc25, Mxd3, Ckap2l, Cenpe, Ndc80, Kif15, Incenp, Mis18bp1 
    ## PC_ 3 
    ## Positive:  Cpa3, Lat2, Gata2, Cyp11a1, Cma1, Hdc, Itgb1, Ms4a2, Atp1b1, Cd200r4 
    ##     Clec12b, Mcpt4, Fam105a, Cst7, Tpsb2, Gm11697, Chst1, Oaf, Cd63, Tubb2a 
    ##     Slc6a4, Smpx, Slc45a3, Nkg7, P2ry14, Car8, Tpsab1, Kif23, Smarca1, Mcpt2 
    ## Negative:  Tcrg-C1, Cysltr1, Il2ra, Icos, Myo6, Rxrg, Il13, Mdfic, Neb, Bcl11b 
    ##     Stc2, Ly6a, Rbp1, Gata3, Il17rb, 1500009L16Rik, Ptpn13, Arl5a, Xlr4a, Pcca 
    ##     Tcrg-C4, Tnfrsf18, Tnfsf11, Tnfrsf4, Ahcyl2, Nampt, Avpi1, Akr1c13, Fabp5, Ccr6 
    ## PC_ 4 
    ## Positive:  Cd81, Nusap1, Serpinb1a, Prc1, Tpx2, Aurka, Hmmr, Cdca8, Plk1, Ube2c 
    ##     Fcnb, Camp, Top2a, Ccnb1, Chil3, Cenpf, Ngp, Ccna2, Ckap2l, Kif11 
    ##     Cdk1, Cdca3, Spc25, Ltf, Homer2, H2afx, Hist1h3c, Cenpe, Kpna2, Aurkb 
    ## Negative:  Klrd1, Gzma, Ccl5, Ctsw, Ncr1, Gimap4, AW112010, Slamf7, Cd160, Klrk1 
    ##     Serpinb6b, Car2, Cxcr3, Lck, Klrb1a, Ifng, Klrb1c, Cd2, Cd244, Nkg7 
    ##     Klrc1, Il2rb, Tbx21, Klrb1b, Tmem37, Gimap6, Lbhd2, Cd7, Samd3, Prf1 
    ## PC_ 5 
    ## Positive:  Ccl9, 1110032F04Rik, Grm6, Tbc1d4, Rnase4, Wdr95, Itgb7, C3ar1, Adgre1, Hgf 
    ##     Alox15, Fut8, Mboat1, Cdh1, Itga1, Bmp4, Ccl3, Pim1, Cysltr1, Fez1 
    ##     Inpp4b, Itgb3, Asb2, Stk19, Lilr4b, Ms4a2, Neb, Dnase2b, Il7r, Adora2b 
    ## Negative:  Serpinb1a, Cma1, Mcpt4, Mt1, Ramp1, Clec12b, Igfbp4, Tpsb2, Myb, Slc45a3 
    ##     Chst1, Mcpt2, Camp, Syce2, Cd81, Cdca7, Smpx, Nme4, Ngp, Smarca1 
    ##     Krt4, Ltf, Tpsab1, Chil3, Fcnb, Gchfr, Hpgds, Fdx1, Ung, Basp1

We can now plot the first two components using the `DimPlot` function.
The first argument here is the Seurat data object `sample1.seurat`. By
providing the `Reduction = "pca"` argument, `DimPlot` looks in the
object for the PCA we created and assigned above.

``` r
DimPlot(sample1.seurat, reduction = "pca", cols = "red")
```

![](sample1.complete.analysis_files/figure-gfm/plot_pca-1.png)<!-- -->

The `DimHeatmap` function allows us to have a look at the selected PCs
`dim =1:dims` for determined `n` number of cells . Each component will
produce one heatmap, the cells are the columns and the top genes of each
component are the rows. `balanced` arguments plot an equal number of
features as positive and negatives

``` r
dims = 6
n_cells = 500
DimHeatmap(sample1.seurat, dims = 1:dims, cells = n_cells, balanced = TRUE)
```

![](sample1.complete.analysis_files/figure-gfm/plot_pc_heatmaps-1.png)<!-- -->

An easy method to determine the number of PCs that represent most of the
variation in the data set is to plot an Elbowplot with the function
`Elbowplot`.

``` r
ElbowPlot(sample1.seurat, ndims = 100)
```

    ## Warning in ElbowPlot(sample1.seurat, ndims = 100): The object only has
    ## information for 50 reductions

![](sample1.complete.analysis_files/figure-gfm/plot_pc_std_dev-1.png)<!-- -->

In the sample1 (homozygous), it seems that 20PCs covers most of the
variation in the data.

#### t-SNE

##### Determination of the good number of PCs

``` r
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca")

DimPlot(smartseq2, reduction = "tsne")
```

![](sample1.complete.analysis_files/figure-gfm/dim_red_tSNE-1.png)<!-- -->
This is the t-SNE representation with the default parameters, I want now
to have a look if differents number of PCs influences the representation
I get and compare it to the selected `dims` (20). For this I will
iterate through different number of PCs.

``` r
# PC_1 to PC_5
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:5)
tsne_mult <- list(DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("5 PCs"))

# PC_1 to PC_10
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:10)
tsne_mult[[2]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("10 PCs")

# PC_1 to PC_30
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20)
tsne_mult[[3]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("20 PCs")

# PC_1 to PC_100
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50)
tsne_mult[[4]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50 PCs")

smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:100)
```

    ## Warning in `[[.DimReduc`(args$object, cells, args$dims): The following
    ## embeddings are not present: NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    ## NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    ## NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA

``` r
tsne_mult[[5]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("100 PCs")

smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:150)
```

    ## Warning in `[[.DimReduc`(args$object, cells, args$dims): The following
    ## embeddings are not present: NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    ## NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    ## NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    ## NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    ## NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    ## NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA

``` r
tsne_mult[[6]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("150 PCs")

CombinePlots(plots = tsne_mult)
```

![](sample1.complete.analysis_files/figure-gfm/dim_red_tSNE_multPCA-1.png)<!-- -->

There are not many differences between 50PCs and more. But surprisingly
the former `dims =20` doesn’t seems to be appropriate for this dataset.

##### Determination of the good number of iterations

t-SNE has a few hyper-parameters that can be tuned for better
visualization. There is an [excellent
tutorial](https://distill.pub/2016/misread-tsne/). The main parameter is
the perplexity, basically indicating how many neighbors to look at. I
will run different perplexities to see the effect, the default
perplexity is 30.

``` r
# Perplexity 3
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20, perplexity = 3)
tsne_mult <- list(DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("20PCs, Perplexity 3"))

# Perplexity 10
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20, perplexity = 10)
tsne_mult[[2]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("20PCs, Perplexity 10")

# Perplexity 30
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20, perplexity = 30)
tsne_mult[[3]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("20PCs, Perplexity 30")

# Perplexity 200
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20, perplexity = 200)
tsne_mult[[4]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("20PCs, Perplexity 200")

CombinePlots(plots = tsne_mult)
```

![](sample1.complete.analysis_files/figure-gfm/dim_red_tSNE_perplexity-1.png)<!-- -->

``` r
# Perplexity 3
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50, perplexity = 3)
tsne_mult <- list(DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50PCs, Perplexity 3"))

# Perplexity 10
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50, perplexity = 10)
tsne_mult[[2]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50PCs, Perplexity 10")

# Perplexity 30
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50, perplexity = 30)
tsne_mult[[3]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50PCs, Perplexity 30")

# Perplexity 500
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50, perplexity = 200)
tsne_mult[[4]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50PCs, Perplexity 200")

CombinePlots(plots = tsne_mult)
```

![](sample1.complete.analysis_files/figure-gfm/dim_red_tSNE_perplexity%2050-1.png)<!-- -->
In both case the defaults value of 30 for the perplexity seems to works
nicely, I will maintain `perplexity = 30`.

##### Determining the number of iterations

``` r
# 100 iterations
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20, max_iter = 100)
tsne_mult <- list(DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("100 iterations"))

# 500 iterations
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20, max_iter = 500)
tsne_mult[[2]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("500 iterations")

# 1000 iterations
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20, max_iter = 1000)
tsne_mult[[3]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("1000 iterations")

# 2000 iterations
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:20, max_iter = 2000)
tsne_mult[[4]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("2000 iterations")

CombinePlots(plots = tsne_mult)
```

![](sample1.complete.analysis_files/figure-gfm/dim_red_tSNE_iterations%2020-1.png)<!-- -->

``` r
# 100 iterations
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50, max_iter = 100)
tsne_mult <- list(DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50PCs,100 iterations"))

# 500 iterations
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50, max_iter = 500)
tsne_mult[[2]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50PCs,500 iterations")

# 1000 iterations
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50, max_iter = 1000)
tsne_mult[[3]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50PCs,1000 iterations")

# 2000 iterations
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:50, max_iter = 2000)
tsne_mult[[4]] <- DimPlot(smartseq2, reduction = "tsne") + NoLegend() + ggtitle("50PCs,2000 iterations")

CombinePlots(plots = tsne_mult)
```

![](sample1.complete.analysis_files/figure-gfm/dim_red_tSNE_iterations-1.png)<!-- -->

As we see after 100 iterations the main structure becomes apparent, but
there is very little detail. 500 to 2000 iterations all look very
similar, with 500 still a bit more loose than 1000, indicating that the
optimization converges somewhere between 500 and 1000 iterations. It’s
not necessary to run more than 500, `max_iter=1000` for `PCs =50` and
`PCs =20`.

The final parameters for the t-SNE, I will uses are :

``` r
dims=50
max_iter=1000
perplexity = 30
sample1.seurat <- FindNeighbors(sample1.seurat, dims = 1:dims)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
sample1.seurat <- FindClusters(sample1.seurat, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5122
    ## Number of edges: 218684
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8898
    ## Number of communities: 11
    ## Elapsed time: 4 seconds

``` r
sample1.seurat <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:dims, max_iter = max_iter, perplexity=perplexity)
DimPlot(sample1.seurat, reduction = "tsne") + ggtitle("50PCs, 1000 iterations, 30perplexity Heteorzygous")
```

![](sample1.complete.analysis_files/figure-gfm/parameters%20t-SNE-1.png)<!-- -->

##### Determining the resolution of the clustering

The `FindCluster` function offers a parameter `resolution` which
determine the number of cluster identified by Seurat. Sajita’s lab claim
that for 3000 the value of the resolution varies between 0.4 to 1.2.
Optimal resolution usually increases for larger datasets.

``` r
dims=50
max_iter=1000
perplexity = 30
sample1.seurat <- FindNeighbors(sample1.seurat, dims = 1:dims)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
sample1.seurat <- FindClusters(sample1.seurat, resolution = 1)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5122
    ## Number of edges: 218684
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8363
    ## Number of communities: 17
    ## Elapsed time: 4 seconds

``` r
sample1.seurat <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:dims, max_iter = max_iter, perplexity=perplexity)
DimPlot(sample1.seurat, reduction = "tsne") + ggtitle("50PCs, 1000 iterations, 30perplexity, homozygous")
```

![](sample1.complete.analysis_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# Resolution 0.4
smartseq2 <- FindClusters(sample1.seurat, resolution = 0.4)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5122
    ## Number of edges: 218684
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9033
    ## Number of communities: 11
    ## Elapsed time: 3 seconds

``` r
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:dims, max_iter = max_iter, perplexity=perplexity)
tsne_mult <- list(DimPlot(smartseq2, reduction = "tsne") + ggtitle("Selected paramters, resolution =0.4"))

# Resolution 0.6
smartseq2 <- FindClusters(sample1.seurat, resolution = 0.6)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5122
    ## Number of edges: 218684
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8777
    ## Number of communities: 11
    ## Elapsed time: 3 seconds

``` r
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:dims, max_iter = max_iter, perplexity=perplexity)
tsne_mult[[2]] <- DimPlot(smartseq2, reduction = "tsne") + ggtitle("Selected paramters, resolution =0.6")

# Resolution 0.8
smartseq2 <- FindClusters(sample1.seurat, resolution = 0.8)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5122
    ## Number of edges: 218684
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8554
    ## Number of communities: 15
    ## Elapsed time: 4 seconds

``` r
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:dims, max_iter = max_iter, perplexity=perplexity)
tsne_mult[[3]] <- DimPlot(smartseq2, reduction = "tsne") + ggtitle("Selected paramters, resolution =0.8")

# Resolution 1
smartseq2 <- FindClusters(sample1.seurat, resolution = 1)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5122
    ## Number of edges: 218684
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8363
    ## Number of communities: 17
    ## Elapsed time: 4 seconds

``` r
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:dims, max_iter = max_iter, perplexity=perplexity)
tsne_mult[[4]] <- DimPlot(smartseq2, reduction = "tsne") + ggtitle("Selected paramters, resolution =1")

# Resolution 1.2
smartseq2 <- FindClusters(sample1.seurat, resolution = 1.2)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5122
    ## Number of edges: 218684
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8204
    ## Number of communities: 18
    ## Elapsed time: 4 seconds

``` r
smartseq2 <- RunTSNE(sample1.seurat, reduction = "pca", dims = 1:dims, max_iter = max_iter, perplexity=perplexity)
tsne_mult[[5]] <- DimPlot(smartseq2, reduction = "tsne")+ ggtitle("Selected paramters, resolution =1.2")

a<-CombinePlots(plots = tsne_mult, ncol = 2)

ggsave(a, filename = "~/Desktop/scRNA-seq/October reanalysis/determining_resolution_sample1.pdf", height = 20, width = 20)

a
```

![](sample1.complete.analysis_files/figure-gfm/resolution%20test-1.png)<!-- -->

The resolution paramters doesn’t change anything to the represenation of
the clusters. I will maintain `resolution=0.5`.

#### UMAP

``` r
smartseq2 <- RunUMAP(sample1.seurat, dims = 1:20)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 03:24:48 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 03:24:48 Read 5122 rows and found 20 numeric columns

    ## 03:24:48 Using Annoy for neighbor search, n_neighbors = 30

    ## 03:24:48 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 03:24:51 Writing NN index file to temp file /var/folders/wx/8gvk2hs93mg75l4qvlr032d00000gn/T//RtmpHBrRRk/file7f693b9b1f36
    ## 03:24:51 Searching Annoy index using 1 thread, search_k = 3000
    ## 03:24:58 Annoy recall = 100%
    ## 03:24:58 Commencing smooth kNN distance calibration using 1 thread
    ## 03:24:59 Initializing from normalized Laplacian + noise
    ## 03:25:05 Commencing optimization for 500 epochs, with 204476 positive edges
    ## 03:25:24 Optimization finished

``` r
DimPlot(smartseq2, reduction = "umap")
```

![](sample1.complete.analysis_files/figure-gfm/dim_red_UMAP-1.png)<!-- -->

``` r
smartseq2 <- RunUMAP(sample1.seurat, dims = 1:50)
```

    ## 03:25:25 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 03:25:25 Read 5122 rows and found 50 numeric columns

    ## 03:25:25 Using Annoy for neighbor search, n_neighbors = 30

    ## 03:25:25 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 03:25:28 Writing NN index file to temp file /var/folders/wx/8gvk2hs93mg75l4qvlr032d00000gn/T//RtmpHBrRRk/file7f694115c487
    ## 03:25:28 Searching Annoy index using 1 thread, search_k = 3000
    ## 03:25:33 Annoy recall = 100%
    ## 03:25:33 Commencing smooth kNN distance calibration using 1 thread
    ## 03:25:34 Initializing from normalized Laplacian + noise
    ## 03:25:40 Commencing optimization for 500 epochs, with 218866 positive edges
    ## 03:26:00 Optimization finished

``` r
DimPlot(smartseq2, reduction = "umap")
```

![](sample1.complete.analysis_files/figure-gfm/dim_red_UMAP%2050-1.png)<!-- -->

``` r
# PC_1 to PC_5
#smartseq2 <- RunUMAP(sample1.seurat, reduction = "pca", dims = 1:5)
#umap_mult <- list(DimPlot(smartseq2, reduction = "umap") + NoLegend() + ggtitle("5 PCs"))

# PC_1 to PC_10
#smartseq2 <- RunUMAP(sample1.seurat, reduction = "pca", dims = 1:10)
#umap_mult[[2]] <- DimPlot(smartseq2, reduction = "umap") + NoLegend() + ggtitle("10 PCs")

# PC_1 to PC_30
#smartseq2 <- RunUMAP(sample1.seurat, reduction = "pca", dims = 1:20)
#umap_mult[[3]] <- DimPlot(smartseq2, reduction = "umap") + NoLegend() + ggtitle("20 PCs")

# PC_1 to PC_100
#smartseq2 <- RunUMAP(sample1.seurat, reduction = "pca", dims = 1:50)
#umap_mult[[4]] <- DimPlot(smartseq2, reduction = "umap") + NoLegend() + ggtitle("50 PCs")

#smartseq2 <- RunUMAP(sample1.seurat, reduction = "pca", dims = 1:100)
#umap_mult[[5]] <- DimPlot(smartseq2, reduction = "umap") + NoLegend() + ggtitle("100 PCs")
#CombinePlots(plots = umap_mult)
```

Most of the identification of the cluster is done, now I want to assign
cell identity to each cluster. I can save the RDS file for the
homozygous sample without cell cycle regression.

``` r
saveRDS(sample1.seurat, file="~/Desktop/scRNA-seq/October reanalysis/sample1.seurat.wo.cell_cycle.rds")
```

### 5\. Assign cell identity

#### Finding differentially expressed feature

``` r
sample1.seurat.markers <- FindAllMarkers(sample1.seurat, min.pct = 0.25, logfc.threshold = 0.25)
```

    ## Calculating cluster 0

    ## Calculating cluster 1

    ## Calculating cluster 2

    ## Calculating cluster 3

    ## Calculating cluster 4

    ## Calculating cluster 5

    ## Calculating cluster 6

    ## Calculating cluster 7

    ## Calculating cluster 8

    ## Calculating cluster 9

    ## Calculating cluster 10

    ## Calculating cluster 11

    ## Calculating cluster 12

    ## Calculating cluster 13

    ## Calculating cluster 14

    ## Calculating cluster 15

    ## Calculating cluster 16

``` r
sample1.seurat.markers %>% group_by(cluster) %>% top_n(n=2, wt = avg_logFC)
```

    ## # A tibble: 34 x 7
    ## # Groups:   cluster [17]
    ##        p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene   
    ##        <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
    ##  1 0.             2.05 0.612 0.099 0.        0       Il1b   
    ##  2 0.             1.89 0.687 0.147 0.        0       Bcl2a1a
    ##  3 1.62e-206      1.45 0.998 0.501 2.48e-202 1       Retnlg 
    ##  4 1.45e-184      1.11 1     0.499 2.23e-180 1       Stfa2  
    ##  5 7.36e-261      1.71 1     0.501 1.13e-256 2       Retnlg 
    ##  6 5.07e-242      1.40 1     0.436 7.79e-238 2       Stfa3  
    ##  7 2.58e-266      2.02 0.935 0.27  3.96e-262 3       Gzma   
    ##  8 1.90e-246      2.52 0.845 0.225 2.93e-242 3       Ccl5   
    ##  9 0.             2.19 0.99  0.082 0.        4       Tcrg-C1
    ## 10 0.             2.04 0.939 0.134 0.        4       Ly6a   
    ## # … with 24 more rows

![](sample1.complete.analysis_files/figure-gfm/Heatmap%20markers-1.png)<!-- -->
Here I plot the expression of the top markers on the cell clusters, we
can see that some genes, like Rag1 are shared between the clusters and
indicate similarity between the cells.

``` r
top_markers <-top10$gene


init=1
end=init
pdf("~/Desktop/scRNA-seq/October reanalysis/sample1_feature.plot.markers.top10markers.pdf",width=20,height=20) #change name or will overwrite
for(index in seq(9,length(top_markers),5)){
  end=index
  plot1 <- VlnPlot(
  object = sample1.seurat,
  features = as.character(top_markers[init:end]),
  ncol = 2,
  pt.size = 0.1
)
plot2 <- FeaturePlot(
  object = sample1.seurat,
  features = as.character(top_markers[init:end]),
  ncol = 2,
  pt.size = 0.1,
  max.cutoff = 'q95'
)
init=end+1
p=CombinePlots(list(plot1,plot2))
print(p)
}

dev.off()
```

    ## quartz_off_screen 
    ##                 2

Here I plot the Tconv markers identified by Kernfeld et al.2018 to show
that in the homozygous most of the cells display the latest Tconv
markers (Tconv5 and Tconv4). This indicates that the homozygous samples
develop mature T-cells

Observation of the homozygous sample with Cerebro Shiny App showed that
the cells cluster also according to the stage of the cell cycle they
belong to. In order to check if this observation is due to an artifact
of the cell cycle or due to biology I performed a cell cycle regression.

## Analysis with Cell cycle regression

``` r
sample1.seurat
```

    ## An object of class Seurat 
    ## 15365 features across 5122 samples within 1 assay 
    ## Active assay: RNA (15365 features)
    ##  2 dimensional reductions calculated: pca, tsne

This is the set of genes, Kernfeld and al.2018 used to normalize their
data for the cell cycle.

``` r
S_genes <- c("Abcc5","Abhd10","Ankrd18a","Asf1b","Atad2","Bbs2","Bivm","Blm","Bmi1","Brca1","Brip1","C5orf42","C11orf82","Cald1","Calm2","Casp2","Ccdc14","Ccdc84","Ccdc150","Cdc7","Cdc45","Cdca5","Cdkn2aip","Cenpm","Cenpq","Cers6","Chml","Coq9","Cpne8","Crebzf","Crls1","Dcaf16","Depdc7","Dhfr","Dna2","Dnajb4","Donson","Dscc1","Dync1li2","E2f8","Eif4ebp2","Enosf1","Esco2","Exo1","Ezh2","Fam178a","Fanca","Fanci","Fen1","Gclm","Golga8a","Golga8b","H1f0","Hells","Hist1h2ac","Hist1h4c","Ints7","Kat2a","Kat2b","Kdelc1","Kiaa1598","Lmo4","Lyrm7","Man1a2","Map3k2","Mastl","Mbd4","Mcm8","Mlf1ip","Mycbp2","Nab1","Neat1","Nfe2l2","Nrd1","Nsun3","Nt5dc1","Nup160","Ogt","Orc3","Osgin2","Phip","Phtf1","Phtf2","Pkmyt1","Pola1","Prim1","Ptar1","Rad18","Rad51","Rad51ap1","Rbbp8","Reep1","Rfc2","Rhobtb3","Rmi1","Rpa2","Rrm1","Rrm2","Rsrc2","Sap30bp","Slc38a2","Sp1","Srsf5","Svip","Top2a","Ttc31","Ttll7","Tyms","Ube2t","Ubl3","Usp1","Zbed5","Zwint")

G2_M_genes <- c("Anln","Ap3d1","Arhgap19","Arl4a","Armc1","Asxl1","Atl2","Aurkb","Bclaf1","Bora","Brd8","Bub3","C2orf69","C14orf80","Casp3","Cbx5","Ccdc107","Ccna2","Ccnf","Cdc16","Cdc25c","Cdca2","Cdca3","Cdca8","Cdk1","Cdkn1b","Cdkn2c","Cdr2","Cenpl","Cep350","Cfd","Cflar","Chek2","Ckap2","Ckap2l","Cyth2","Dcaf7","Dhx8","Dnajb1","Entpd5","Espl1","Fadd","Fam83d","Fan1","Fancd2","G2e3","Gabpb1","Gas1","Gas2l3","H2afx","Haus8","Hint3","Hipk2","Hjurp","Hmgb2","Hn1","Hp1bp3","Hrsp12","Ifnar1","Iqgap3","Katna1","Kctd9","Kdm4a","Kiaa1524","Kif5b","Kif11","Kif20b","Kif22","Kif23","Kifc1","Klf6","Kpna2","Lbr","Lix1l","Lmnb1","Mad2l1","Malat1","Melk","Mgat2","Mid1","Mis18bp1","Mnd1","Ncapd3","Ncaph","Ncoa5","Ndc80","Neil3","Nfic","Nipbl","Nmb","Nr3c1","Nucks1","Numa1","Nusap1","Pif1","Pknox1","Polq","Ppp1r2","Psmd11","Psrc1","Rangap1","Rccd1","Rdh11","Rnf141","Sap30","Ska3","Smc4","Stat1","Stil","Stk17b","Suclg2","Tfap2a","Timp1","Tmem99","Tmpo","Tnpo2","Top2a","Traip","Trim59","Trmt2a","Ttf2","Tuba1a","Tubb","Tubb2a","Tubb4b","Tubd1","Uaca","Ube2c","Vps25","Vta1","Wsb1","Znf587","Znhit2")
```

``` r
sample1.seurat <- CellCycleScoring(sample1.seurat, s.features = S_genes, g2m.features = G2_M_genes, set.ident = TRUE)
```

    ## Warning: The following features are not present in the object: Ankrd18a,
    ## C5orf42, C11orf82, Dcaf16, Enosf1, Fam178a, Golga8a, Golga8b, Kiaa1598,
    ## Mlf1ip, Ttc31, Ttll7, not searching for symbol synonyms

    ## Warning: The following features are not present in the object: C2orf69,
    ## C14orf80, Cfd, Gas1, Hn1, Hrsp12, Kiaa1524, Timp1, Tmem99, Tubb, Znf587,
    ## not searching for symbol synonyms

``` r
head(sample1.seurat[[]])
```

    ##                        orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCCAAGAGTATAC Thymus_homo_5582      21685         3056 0.39658750
    ## AAACCCAAGATTGGGC Thymus_homo_5582      17267         3728 3.04048184
    ## AAACCCAAGGTAGTAT Thymus_homo_5582      14334         3471 1.93944468
    ## AAACCCACAATTGAAG Thymus_homo_5582      21514         4257 2.65408571
    ## AAACCCACACACACGC Thymus_homo_5582       3182         1048 0.15713388
    ## AAACCCAGTGGTTTAC Thymus_homo_5582       3013          930 0.09956854
    ##                  percent.rp RNA_snn_res.0.5 seurat_clusters RNA_snn_res.1
    ## AAACCCAAGAGTATAC   1.918377               4               5             5
    ## AAACCCAAGATTGGGC  26.680952               5               4             4
    ## AAACCCAAGGTAGTAT  22.185015               3               6             6
    ## AAACCCACAATTGAAG  24.500325               3               6             6
    ## AAACCCACACACACGC   2.231301               1               1             1
    ## AAACCCAGTGGTTTAC   1.526718               0               0             0
    ##                      S.Score    G2M.Score Phase old.ident
    ## AAACCCAAGAGTATAC -0.04957815 -0.014357925    G1         5
    ## AAACCCAAGATTGGGC  0.03583419 -0.023230254     S         4
    ## AAACCCAAGGTAGTAT  0.01464740  0.246256700   G2M         6
    ## AAACCCACAATTGAAG  0.06726391  0.276110819   G2M         6
    ## AAACCCACACACACGC  0.02694748 -0.008335882     S         1
    ## AAACCCAGTGGTTTAC  0.01362137  0.019976468   G2M         0

Now each cell is assigned to a cell cycle stage. I proceed to the
regression of the cell cycle score

``` r
test <- RunPCA(sample1.seurat, features = c(S_genes, G2_M_genes))
```

    ## Warning in PrepDR(object = object, features = features, verbose = verbose):
    ## The following 23 features requested have not been scaled (running reduction
    ## without them): Ankrd18a, C5orf42, C11orf82, Dcaf16, Enosf1, Fam178a,
    ## Golga8a, Golga8b, Kiaa1598, Mlf1ip, Ttc31, Ttll7, C2orf69, C14orf80, Cfd,
    ## Gas1, Hn1, Hrsp12, Kiaa1524, Timp1, Tmem99, Tubb, Znf587

    ## PC_ 1 
    ## Positive:  Malat1, Neat1, Svip, Klf6, Ubl3, Ogt, Cers6, Ppp1r2, Stk17b, Asxl1 
    ##     Stat1, Lmo4, Nfe2l2, Cdkn1b, Nfic, Rnf141, Sap30, Abcc5, Nsun3, Srsf5 
    ##     Cflar, Wsb1, Nmb, Hipk2, Nr3c1, Kdm4a, Tubb2a, Cald1, Tfap2a, Kat2b 
    ## Negative:  Top2a, H2afx, Cdca8, Cdk1, Nucks1, Asf1b, Ccna2, Rrm2, Nusap1, Aurkb 
    ##     Ezh2, Smc4, Kif11, Cdca3, Tyms, Ube2c, Rrm1, Cbx5, Ndc80, Esco2 
    ##     Ckap2l, Kif20b, Mis18bp1, Kif22, Lmnb1, Tubb4b, Kpna2, Rangap1, Cdca2, Cenpm 
    ## PC_ 2 
    ## Positive:  Ube2c, Nusap1, Malat1, Ckap2l, Hmgb2, Kpna2, Kif22, Mis18bp1, Kif11, Arhgap19 
    ##     Pif1, Cenpl, Ccna2, Kif20b, Neat1, Cdca8, Cdca3, Tubb4b, Klf6, Cdc25c 
    ##     Kifc1, Kif23, Gas2l3, Ckap2, Ccnf, Ndc80, H2afx, Calm2, Anln, Cdca2 
    ## Negative:  Hells, Rpa2, Dhfr, Prim1, Rfc2, Dscc1, Fen1, Atad2, Exo1, Pola1 
    ##     Mgat2, Tyms, Cbx5, Donson, Arl4a, Usp1, Brca1, Suclg2, Psmd11, Kat2a 
    ##     Rad51, Zwint, Rrm2, Armc1, Trmt2a, Rrm1, Cdc45, Bclaf1, Ccdc107, Mycbp2 
    ## PC_ 3 
    ## Positive:  Hmgb2, Neat1, Esco2, Ppp1r2, Svip, Rrm2, Brca1, Rad51ap1, Cers6, E2f8 
    ##     Cdc45, Tmpo, Sap30, Rad51, Klf6, Stat1, Neil3, Tyms, Exo1, Asxl1 
    ##     Malat1, Ubl3, Dhfr, Brip1, Stk17b, Asf1b, Pola1, Cdca5, Atad2, Ogt 
    ## Negative:  Tubb2a, Nmb, Ccdc107, Arl4a, Hp1bp3, H1f0, Nucks1, Mgat2, Ube2c, Hint3 
    ##     Psmd11, G2e3, Mycbp2, Armc1, Psrc1, Ckap2, Tuba1a, Kif23, Coq9, Tnpo2 
    ##     Gas2l3, Atl2, Kif5b, Kpna2, Vta1, Rfc2, Rdh11, Kat2a, Katna1, Cdc25c 
    ## PC_ 4 
    ## Positive:  Stk17b, Malat1, Ogt, Neat1, Cdkn1b, Stat1, Nrd1, Calm2, Slc38a2, Srsf5 
    ##     Rsrc2, Mycbp2, Cers6, Nipbl, Klf6, Tubb2a, Hist1h2ac, Ifnar1, Numa1, Eif4ebp2 
    ##     Cep350, Wsb1, Bclaf1, Phip, Tmpo, Kif23, Sp1, Blm, Casp3, Nab1 
    ## Negative:  Hipk2, Abcc5, H1f0, Hmgb2, Lmo4, Arhgap19, Cdc45, Abhd10, Cdca8, Cdca3 
    ##     Rrm2, Aurkb, Iqgap3, H2afx, Anln, Ube2c, Cdk1, Ccna2, Asf1b, Katna1 
    ##     Kif22, Donson, Ccnf, Entpd5, Pkmyt1, Crls1, Esco2, Rrm1, Fanci, Kpna2 
    ## PC_ 5 
    ## Positive:  Hmgb2, Srsf5, Hp1bp3, Lmnb1, Stat1, Tuba1a, Nucks1, Gclm, Rangap1, Hells 
    ##     Ppp1r2, Entpd5, Ubl3, Ckap2, Tmpo, Tubb4b, Rdh11, Rsrc2, Ube2c, Ints7 
    ##     Cpne8, Mad2l1, Rpa2, Cenpq, Phip, Cdca8, Sap30, Tnpo2, Cep350, Ccdc84 
    ## Negative:  Nmb, Tubb2a, Calm2, Malat1, Esco2, Kif23, H1f0, Hist1h2ac, Rfc2, E2f8 
    ##     Blm, Cdca5, Neil3, Casp3, Ncaph, Aurkb, Rad51ap1, Brca1, G2e3, Anln 
    ##     Stil, Brip1, Ccnf, Tyms, Rrm2, Nrd1, Top2a, Mis18bp1, Ncapd3, Atl2

``` r
DimPlot(test)
```

![](sample1.complete.analysis_files/figure-gfm/PCA%20cell%20cycle-1.png)<!-- -->
As clearly shown on the t-SNE plot above the cells cluster according to
their cell cycle stage.

``` r
sample1.seurat<- RunPCA(sample1.seurat, features = VariableFeatures(object=sample1.seurat), npcs =50)
```

    ## PC_ 1 
    ## Positive:  Gm5483, Hp, Stfa2l1, BC100530, Spi1, Pglyrp1, Wfdc21, Stfa2, Lyz2, Alox5ap 
    ##     Gm5416, Anxa1, Mcemp1, Asprv1, Prdx5, F630028O10Rik, Olfm4, Chil1, Mmp8, Cebpd 
    ##     Slfn1, Retnlg, Raet1e, Stfa3, Cd3d, Mgst1, Cxcr2, Mmp9, Igsf6, Mxd1 
    ## Negative:  Ptma, Npm1, Rpl3, Rpl12, Rpl14, Ran, Rpl4, Hint1, Lgals1, Rps2 
    ##     Tmsb10, Rpl36a, Nme1, Rpl5, Eef1d, Rps17, Rps6, Set, Ptprcap, Snrpg 
    ##     Erh, mt-Nd1, Hspe1, Rpl31, Serbp1, Eef1b2, Ranbp1, Snrpd1, Rpl22l1, Atp5g1 
    ## PC_ 2 
    ## Positive:  Cpa3, Cyp11a1, Gata2, Hdc, Cma1, Gm11697, Kit, Mcpt4, Oaf, Ms4a2 
    ##     Clec12b, Tpsb2, Slc7a8, Chst1, Slc45a3, Smpx, Smarca1, Mcpt2, Cd63, Hpgds 
    ##     Hs3st1, Fam105a, Slc6a4, Atp8b5, Tpsab1, Krt4, Tubb2a, Gchfr, Basp1, Lxn 
    ## Negative:  Ctsw, Klrd1, AW112010, Cd160, Il2rb, Cd3g, S100a4, Ncr1, Gzma, Lck 
    ##     Slamf7, Gimap4, Cd2, Cxcr3, Ly6c2, Ccl5, Car2, Klrk1, Serpinb6b, Tmem37 
    ##     Txk, Ifitm3, Sh2d2a, Ifngr1, Ifng, Klrc1, Dusp2, Klrb1a, Trbc2, Evl 
    ## PC_ 3 
    ## Positive:  Klrd1, Ctsw, Gzma, Ncr1, Nkg7, Car2, Serpinb6b, Klrc1, Ccl5, AW112010 
    ##     Cxcr3, Slamf7, Cd2, Klrk1, Klrb1a, Itgb1, Gimap4, Cd160, Serpinb9, Cd244 
    ##     Klrb1c, Klrb1b, Samd3, Lbhd2, Klrc2, Tmem37, Tbx21, Lck, Ifng, Atp1b1 
    ## Negative:  Tcrg-C1, Cysltr1, Il2ra, Icos, Myo6, Ptpn13, Il13, Rxrg, Neb, Ly6a 
    ##     1500009L16Rik, Stc2, Bcl11b, Il17rb, Mdfic, Rbp1, Rnf128, Akr1c13, Gata3, Homer2 
    ##     Kcnn4, Arl5a, Tnfrsf4, Il7r, Tnfsf11, Pcca, Ahcyl2, Fam213a, 1700061F12Rik, Ccr6 
    ## PC_ 4 
    ## Positive:  Camp, Serpinb1a, Ltf, Ngp, Chil3, Fcnb, Cybb, BC117090, Lcn2, Lbp 
    ##     Anxa3, Hmgn2, Ifitm6, Dstn, Adpgk, Cebpe, Cd177, Orm1, Igfbp4, Ms4a3 
    ##     Ceacam1, Mgst2, Lta4h, Tmem216, Aldh2, Golim4, Gca, St3gal5, Stfa3, Cd81 
    ## Negative:  Ccl6, Cyp4f18, Pla2g7, Emilin2, Rnase4, Ccl9, Apoe, 1110032F04Rik, Tbc1d4, Itgb7 
    ##     Wdr95, Grm6, Csf1, Ifitm1, Bcl2a1a, Il1b, Pim1, Ccl3, Stat1, C3ar1 
    ##     Adgre1, Cxcr2, Bcl2a1b, Marcks, Epcam, Csf3r, Alox15, Fut8, Aqp9, S100a10 
    ## PC_ 5 
    ## Positive:  Ccl9, 1110032F04Rik, Grm6, Tbc1d4, Hgf, C3ar1, Alox15, Wdr95, Adgre1, Ms4a2 
    ##     Mcpt8, Cdh1, Adora2b, Bmp4, Lilr4b, P2rx1, Fez1, Syne1, Dnase2b, Ccl3 
    ##     Ctse, Itga1, Plac8, Stk19, F13a1, Csf2rb2, Cd200r3, Igfbp7, Cd200r4, Cela3b 
    ## Negative:  Cma1, Mcpt4, Clec12b, Tpsb2, Chst1, Tpsab1, Rab27b, Slc45a3, Smpx, Krt4 
    ##     Hpgds, Mcpt2, Smarca1, Basp1, Kit, Gchfr, Gnai1, Atp8b5, Lxn, Tph1 
    ##     Clec10a, Rgs13, Fam198b, Cited4, Ddah2, Rgs10, Maob, Fam105a, Fjx1, Pde3a

``` r
sample1.seurat <- FindNeighbors(sample1.seurat, dims = 1:dims)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
sample1.seurat <- FindClusters(sample1.seurat, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5122
    ## Number of edges: 205255
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8862
    ## Number of communities: 10
    ## Elapsed time: 3 seconds

#### T-SNE representation after cell cycle regression

``` r
sample1.seurat <- RunTSNE(sample1.seurat, dims = 1:dims)
DimPlot(sample1.seurat, reduction = "tsne", label = TRUE)
```

![](sample1.complete.analysis_files/figure-gfm/t-SNE%20representation%20after%20cell%20cycle%20regression-1.png)<!-- -->

``` r
test <- RunPCA(sample1.seurat, features = c(S_genes, G2_M_genes))
```

    ## Warning in PrepDR(object = object, features = features, verbose = verbose):
    ## The following 23 features requested have not been scaled (running reduction
    ## without them): Ankrd18a, C5orf42, C11orf82, Dcaf16, Enosf1, Fam178a,
    ## Golga8a, Golga8b, Kiaa1598, Mlf1ip, Ttc31, Ttll7, C2orf69, C14orf80, Cfd,
    ## Gas1, Hn1, Hrsp12, Kiaa1524, Timp1, Tmem99, Tubb, Znf587

    ## PC_ 1 
    ## Positive:  Malat1, Neat1, Klf6, Svip, Cers6, Ppp1r2, Ubl3, Ogt, Stk17b, Cdkn1b 
    ##     Asxl1, Stat1, Sap30, Nfe2l2, Lmo4, Nfic, Rnf141, Hmgb2, Srsf5, Wsb1 
    ##     Nr3c1, Lbr, Abcc5, Slc38a2, Cflar, Cep350, Nsun3, Kdm4a, Kat2b, Hipk2 
    ## Negative:  Nucks1, Top2a, Rrm2, Ezh2, Tyms, Cdca8, Asf1b, Hells, H2afx, Cbx5 
    ##     Rpa2, Fen1, Ccna2, Atad2, Prim1, Rrm1, Cdk1, Smc4, Rangap1, Ube2c 
    ##     Nusap1, Usp1, Cdca3, Aurkb, Rfc2, Dhfr, Psmd11, Cenpm, Ncapd3, Rad51 
    ## PC_ 2 
    ## Positive:  Nusap1, Calm2, Ube2c, Tubb2a, H1f0, Kif23, Mis18bp1, Ckap2l, Nmb, Kpna2 
    ##     Kif11, G2e3, Kif22, Pif1, Aurkb, Gas2l3, Ndc80, Ccnf, Kif20b, Cdca2 
    ##     Cdk1, Hist1h2ac, H2afx, Ccna2, Fam83d, Cdc25c, Anln, Nrd1, Malat1, Top2a 
    ## Negative:  Hells, Hmgb2, Dhfr, Tmpo, Rpa2, Lmnb1, Stat1, Prim1, Fen1, Ppp1r2 
    ##     Dscc1, Rrm2, Exo1, Atad2, Donson, Cbx5, Brca1, Stk17b, Rad51, Sap30 
    ##     Tyms, Cdc45, Usp1, Rad51ap1, Nr3c1, Asxl1, Ncapd3, Pola1, Klf6, Trmt2a 
    ## PC_ 3 
    ## Positive:  Esco2, Hmgb2, Aurkb, Rrm2, Top2a, Rad51ap1, Ccnf, Neil3, Asf1b, Cdca5 
    ##     Nusap1, Cdk1, Cdc45, Ncaph, Ndc80, Anln, Tyms, Kif11, Rad51, Hipk2 
    ##     E2f8, Abcc5, Cdca2, Cdkn2c, H2afx, Cdca3, Brca1, Mis18bp1, Rrm1, Sap30 
    ## Negative:  Tubb2a, Nrd1, Mycbp2, Hp1bp3, Arl4a, Calm2, Ccdc107, Rsrc2, Bclaf1, Ogt 
    ##     Srsf5, Nipbl, Mgat2, Stk17b, Eif4ebp2, Numa1, Tuba1a, Malat1, Armc1, Stat1 
    ##     Kat2a, Nmb, Kif23, Rfc2, Cdkn1b, Zwint, Psmd11, Ifnar1, Tnpo2, Slc38a2 
    ## PC_ 4 
    ## Positive:  Hmgb2, Ube2c, Cdca8, Tubb4b, Ccna2, Srsf5, Ckap2, Cdca3, Nucks1, Kif20b 
    ##     Arhgap19, Cdc25c, Lmnb1, Hp1bp3, Ubl3, Tuba1a, Mad2l1, Lmo4, Gas2l3, Kpna2 
    ##     Gclm, Ckap2l, Kif22, Rsrc2, Rangap1, Trim59, Nfe2l2, Cenpq, Svip, Psrc1 
    ## Negative:  Nmb, Tubb2a, Esco2, Malat1, Kif23, Casp3, Calm2, Ncaph, Blm, Rfc2 
    ##     Brca1, Rad51ap1, Ncapd3, Nr3c1, E2f8, Tyms, Brip1, Neil3, Cdca5, Aurkb 
    ##     Hist1h2ac, Nipbl, Atad2, Pola1, Rrm2, Stil, Suclg2, Stk17b, Prim1, Top2a 
    ## PC_ 5 
    ## Positive:  Calm2, Ccna2, Tuba1a, Cdca5, Cenpm, Lmnb1, Melk, Malat1, Ndc80, Slc38a2 
    ##     Nab1, Cdca3, Ckap2, Kif23, Cdca8, Mad2l1, Kif20b, Ogt, Neat1, Cdkn1b 
    ##     Kif22, Tubb4b, Tmpo, Neil3, Traip, Hp1bp3, Fancd2, Ttf2, Ska3, Asf1b 
    ## Negative:  Lmo4, Hipk2, Anln, Abcc5, H1f0, Kpna2, Kif5b, Bub3, Cdk1, Arhgap19 
    ##     Entpd5, Hjurp, Gclm, Rnf141, Arl4a, Iqgap3, Cdc45, Sap30, Cflar, Vta1 
    ##     Donson, Dna2, Cpne8, Ubl3, Coq9, Gabpb1, Hint3, G2e3, Ints7, Ap3d1

``` r
DimPlot(test)
```

![](sample1.complete.analysis_files/figure-gfm/t-SNE%20representation%20after%20cell%20cycle%20regression%20with%20markers-1.png)<!-- -->

#### UMAP representation after cell cycle regression

``` r
sample1.seurat <- RunUMAP(sample1.seurat, dims = 1:dims)
```

    ## 04:16:22 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 04:16:22 Read 5122 rows and found 50 numeric columns

    ## 04:16:22 Using Annoy for neighbor search, n_neighbors = 30

    ## 04:16:22 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 04:16:25 Writing NN index file to temp file /var/folders/wx/8gvk2hs93mg75l4qvlr032d00000gn/T//RtmpHBrRRk/file7f695952834a
    ## 04:16:25 Searching Annoy index using 1 thread, search_k = 3000
    ## 04:16:30 Annoy recall = 100%
    ## 04:16:30 Commencing smooth kNN distance calibration using 1 thread
    ## 04:16:31 Initializing from normalized Laplacian + noise
    ## 04:16:36 Commencing optimization for 500 epochs, with 222538 positive edges
    ## 04:16:55 Optimization finished

``` r
DimPlot(sample1.seurat, reduction = "umap", label = TRUE)
```

![](sample1.complete.analysis_files/figure-gfm/UMAP%20representation%20after%20cell%20cycle%20regression-1.png)<!-- -->

``` r
test <- RunPCA(sample1.seurat, features = c(S_genes, G2_M_genes))
```

    ## Warning in PrepDR(object = object, features = features, verbose = verbose):
    ## The following 23 features requested have not been scaled (running reduction
    ## without them): Ankrd18a, C5orf42, C11orf82, Dcaf16, Enosf1, Fam178a,
    ## Golga8a, Golga8b, Kiaa1598, Mlf1ip, Ttc31, Ttll7, C2orf69, C14orf80, Cfd,
    ## Gas1, Hn1, Hrsp12, Kiaa1524, Timp1, Tmem99, Tubb, Znf587

    ## PC_ 1 
    ## Positive:  Malat1, Neat1, Klf6, Svip, Cers6, Ppp1r2, Ubl3, Ogt, Stk17b, Cdkn1b 
    ##     Asxl1, Stat1, Sap30, Nfe2l2, Lmo4, Nfic, Rnf141, Hmgb2, Srsf5, Wsb1 
    ##     Nr3c1, Lbr, Abcc5, Slc38a2, Cflar, Cep350, Nsun3, Kdm4a, Kat2b, Hipk2 
    ## Negative:  Nucks1, Top2a, Rrm2, Ezh2, Tyms, Cdca8, Asf1b, Hells, H2afx, Cbx5 
    ##     Rpa2, Fen1, Ccna2, Atad2, Prim1, Rrm1, Cdk1, Smc4, Rangap1, Ube2c 
    ##     Nusap1, Usp1, Cdca3, Aurkb, Rfc2, Dhfr, Psmd11, Cenpm, Ncapd3, Rad51 
    ## PC_ 2 
    ## Positive:  Nusap1, Calm2, Ube2c, Tubb2a, H1f0, Kif23, Mis18bp1, Ckap2l, Nmb, Kpna2 
    ##     Kif11, G2e3, Kif22, Pif1, Aurkb, Gas2l3, Ndc80, Ccnf, Kif20b, Cdca2 
    ##     Cdk1, Hist1h2ac, H2afx, Ccna2, Fam83d, Cdc25c, Anln, Nrd1, Malat1, Top2a 
    ## Negative:  Hells, Hmgb2, Dhfr, Tmpo, Rpa2, Lmnb1, Stat1, Prim1, Fen1, Ppp1r2 
    ##     Dscc1, Rrm2, Exo1, Atad2, Donson, Cbx5, Brca1, Stk17b, Rad51, Sap30 
    ##     Tyms, Cdc45, Usp1, Rad51ap1, Nr3c1, Asxl1, Ncapd3, Pola1, Klf6, Trmt2a 
    ## PC_ 3 
    ## Positive:  Esco2, Hmgb2, Aurkb, Rrm2, Top2a, Rad51ap1, Ccnf, Neil3, Asf1b, Cdca5 
    ##     Nusap1, Cdk1, Cdc45, Ncaph, Ndc80, Anln, Tyms, Kif11, Rad51, Hipk2 
    ##     E2f8, Abcc5, Cdca2, Cdkn2c, H2afx, Cdca3, Brca1, Mis18bp1, Rrm1, Sap30 
    ## Negative:  Tubb2a, Nrd1, Mycbp2, Hp1bp3, Arl4a, Calm2, Ccdc107, Rsrc2, Bclaf1, Ogt 
    ##     Srsf5, Nipbl, Mgat2, Stk17b, Eif4ebp2, Numa1, Tuba1a, Malat1, Armc1, Stat1 
    ##     Kat2a, Nmb, Kif23, Rfc2, Cdkn1b, Zwint, Psmd11, Ifnar1, Tnpo2, Slc38a2 
    ## PC_ 4 
    ## Positive:  Hmgb2, Ube2c, Cdca8, Tubb4b, Ccna2, Srsf5, Ckap2, Cdca3, Nucks1, Kif20b 
    ##     Arhgap19, Cdc25c, Lmnb1, Hp1bp3, Ubl3, Tuba1a, Mad2l1, Lmo4, Gas2l3, Kpna2 
    ##     Gclm, Ckap2l, Kif22, Rsrc2, Rangap1, Trim59, Nfe2l2, Cenpq, Svip, Psrc1 
    ## Negative:  Nmb, Tubb2a, Esco2, Malat1, Kif23, Casp3, Calm2, Ncaph, Blm, Rfc2 
    ##     Brca1, Rad51ap1, Ncapd3, Nr3c1, E2f8, Tyms, Brip1, Neil3, Cdca5, Aurkb 
    ##     Hist1h2ac, Nipbl, Atad2, Pola1, Rrm2, Stil, Suclg2, Stk17b, Prim1, Top2a 
    ## PC_ 5 
    ## Positive:  Calm2, Ccna2, Tuba1a, Cdca5, Cenpm, Lmnb1, Melk, Malat1, Ndc80, Slc38a2 
    ##     Nab1, Cdca3, Ckap2, Kif23, Cdca8, Mad2l1, Kif20b, Ogt, Neat1, Cdkn1b 
    ##     Kif22, Tubb4b, Tmpo, Neil3, Traip, Hp1bp3, Fancd2, Ttf2, Ska3, Asf1b 
    ## Negative:  Lmo4, Hipk2, Anln, Abcc5, H1f0, Kpna2, Kif5b, Bub3, Cdk1, Arhgap19 
    ##     Entpd5, Hjurp, Gclm, Rnf141, Arl4a, Iqgap3, Cdc45, Sap30, Cflar, Vta1 
    ##     Donson, Dna2, Cpne8, Ubl3, Coq9, Gabpb1, Hint3, G2e3, Ints7, Ap3d1

``` r
DimPlot(test, reduction = "umap")
```

![](sample1.complete.analysis_files/figure-gfm/UMAP%20representation%20after%20cell%20cycle%20regression%20with%20markers-1.png)<!-- -->

``` r
saveRDS(sample1.seurat, file="~/Desktop/scRNA-seq/October reanalysis/sample1.seurat.cell_cycle.regressed.rds")
```

#### Finding differentially expressed feature

``` r
sample1.seurat.markers <- FindAllMarkers(sample1.seurat, min.pct = 0.25, logfc.threshold = 0.25)
```

    ## Calculating cluster 0

    ## Calculating cluster 1

    ## Calculating cluster 2

    ## Calculating cluster 3

    ## Calculating cluster 4

    ## Calculating cluster 5

    ## Calculating cluster 6

    ## Calculating cluster 7

    ## Calculating cluster 8

    ## Calculating cluster 9

``` r
sample1.seurat.markers %>% group_by(cluster) %>% top_n(n=2, wt = avg_logFC)
```

    ## # A tibble: 20 x 7
    ## # Groups:   cluster [10]
    ##        p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene     
    ##        <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
    ##  1 0.             1.92 0.997 0.651 0.        0       Msrb1    
    ##  2 0.             1.92 0.506 0.065 0.        0       Il1b     
    ##  3 0.             2.53 0.899 0.148 0.        1       Ltf      
    ##  4 0.             2.43 0.987 0.387 0.        1       Ngp      
    ##  5 0.             2.48 0.845 0.199 0.        2       Ccl5     
    ##  6 0.             2.15 0.94  0.241 0.        2       Gzma     
    ##  7 0.             2.38 0.968 0.289 0.        3       Hist1h2ae
    ##  8 0.             2.36 0.954 0.239 0.        3       Hist1h1b 
    ##  9 0.             2.19 0.996 0.083 0.        4       Tcrg-C1  
    ## 10 0.             2.00 0.94  0.135 0.        4       Ly6a     
    ## 11 0.             2.84 0.592 0.024 0.        5       Mcpt8    
    ## 12 0.             2.63 0.989 0.024 0.        5       Ccl9     
    ## 13 5.97e-169      1.29 0.587 0.07  9.17e-165 6       Lta      
    ## 14 9.68e-101      1.34 0.884 0.364 1.49e- 96 6       S100a4   
    ## 15 0.             5.57 0.911 0.039 0.        7       Mcpt4    
    ## 16 0.             3.70 0.745 0.006 0.        7       Mcpt2    
    ## 17 0.             2.55 0.588 0.018 0.        8       Prtn3    
    ## 18 9.08e-197      2.40 0.296 0.008 1.39e-192 8       Elane    
    ## 19 1.54e-216      4.43 1     0.043 2.37e-212 9       Tpsb2    
    ## 20 6.61e-135      4.18 1     0.079 1.01e-130 9       Cma1

![](sample1.complete.analysis_files/figure-gfm/Heatmap%20markers%20cell%20cycle%20regression-1.png)<!-- -->
Here I plot the expression of the top markers on the cell clusters, we
can see that some genes, like Rag1 are shared between the clusters and
indicate similarity between the cells.

``` r
top_markers <-top10$gene


init=1
end=init
pdf("~/Desktop/scRNA-seq/October reanalysis/sample1_feature.cellcycleregressed.plot.markers.top10markers.pdf",width=20,height=20) #change name or will overwrite
for(index in seq(9,length(top_markers),5)){
  end=index
  plot1 <- VlnPlot(
  object = sample1.seurat,
  features = as.character(top_markers[init:end]),
  ncol = 2,
  pt.size = 0.1
)
plot2 <- FeaturePlot(
  object = sample1.seurat,
  features = as.character(top_markers[init:end]),
  ncol = 2,
  pt.size = 0.1,
  max.cutoff = 'q95'
)
init=end+1
p=CombinePlots(list(plot1,plot2))
print(p)
}

dev.off()
```

    ## quartz_off_screen 
    ##                 2

Here I plot the Tconv markers identified by Kernfeld et al.2018 to show
that in the homozygous most of the cells display the latest Tconv
markers (Tconv5 and Tconv4). This indicates that the homozygous samples
develop mature T-cells
