README file for: Assessment of pharmacogenomic agreement

Zhaleh Safikhani, Nehme El-Hachem, Rene Quevedo, Petr Smirnov, Anna Goldenberg, Nicolai Juul Birkbak, Christopher Mason, Christos Hatzis, Leming Shi, Hugo JWL Aerts, John Quackenbush, Benjamin Haibe-Kains

1. Copy Data_File_1.R and Data_File_2.csv in one directory
2. Open an R session and set the directory to the one containing Data Files.
3. Run Data_File_1.R
4. Note that it requires installation of PharmcoGx package which is coded in the R file and your R session should be as following after loading PharmacoGx package into your R session:
 sessionInfo()
R version 3.2.2 (2015-08-14)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.11.1 (El Capitan)

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] Biobase_2.30.0       BiocGenerics_0.16.1  PharmacoGx_1.1.4     BiocInstaller_1.20.1

loaded via a namespace (and not attached):
 [1] mclust_5.1             Rcpp_0.12.3            lattice_0.20-33        relations_0.6-6       
 [5] class_7.3-13           piano_1.10.2           gtools_3.5.0           digest_0.6.9          
 [9] slam_0.1-32            plyr_1.8.3             acepack_1.3-3.3        stats4_3.2.2          
[13] RSQLite_1.0.0          e1071_1.6-7            ggplot2_2.1.0          gplots_2.17.0         
[17] iC10_1.1.3             gdata_2.17.0           S4Vectors_0.8.3        rpart_4.1-10          
[21] magicaxis_1.9.4        splines_3.2.2          sets_1.0-16            epibasix_1.3          
[25] genefu_2.2.0           downloader_0.4         foreign_0.8-65         igraph_1.0.1          
[29] RCurl_1.95-4.7         biomaRt_2.26.1         munsell_0.4.3          survcomp_1.20.0       
[33] marray_1.48.0          AIMS_1.2.0             nnet_7.3-10            gridExtra_2.2.1       
[37] prodlim_1.5.5          Hmisc_3.17-2           IRanges_2.4.4          XML_3.98-1.3          
[41] MASS_7.3-43            bitops_1.0-6           SuppDists_1.1-9.1      grid_3.2.2            
[45] gtable_0.2.0           DBI_0.3.1              magrittr_1.5           scales_0.4.0          
[49] KernSmooth_2.23-15     amap_0.8-14            iC10TrainingData_1.0.1 sm_2.2-5.4            
[53] limma_3.26.3           latticeExtra_0.6-28    Formula_1.2-1          rmeta_2.16            
[57] lava_1.4.1             RColorBrewer_1.1-2     tools_3.2.2            pamr_1.55             
[61] plotrix_3.6-1          survival_2.38-3        AnnotationDbi_1.32.0   colorspace_1.2-6      
[65] cluster_2.0.3          caTools_1.17.1         survivalROC_1.0.3      bootstrap_2015.2     