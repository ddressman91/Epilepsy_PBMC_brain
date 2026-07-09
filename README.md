# Epilepsy_PBMC_brain
Code used to analyze data from plasma, PBMCs, and brain tissue from epilepsy patients and healthy controls

# R version and package version info:
R version 4.5.1 (2025-06-13 ucrt)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scRepertoire_2.6.2     janitor_2.2.1          coloc_5.2.3           
 [4] Seurat_5.4.0           SeuratObject_5.3.0     sp_2.2-1              
 [7] openxlsx_4.2.8.1       purrr_1.2.1            rlang_1.1.7           
[10] readxl_1.4.5           scales_1.4.0           EnhancedVolcano_1.28.2
[13] ggrepel_0.9.8          ggplot2_4.0.2          stringr_1.6.0         
[16] clipr_0.8.0            org.Hs.eg.db_3.22.0    AnnotationDbi_1.72.0  
[19] IRanges_2.44.0         S4Vectors_0.48.1       Biobase_2.70.0        
[22] BiocGenerics_0.56.0    generics_0.1.4         clusterProfiler_4.18.4
[25] dplyr_1.2.0           

loaded via a namespace (and not attached):
  [1] fs_2.1.0                    matrixStats_1.5.0          
  [3] spatstat.sparse_3.1-0       enrichplot_1.30.5          
  [5] lubridate_1.9.5             httr_1.4.8                 
  [7] RColorBrewer_1.1-3          tools_4.5.1                
  [9] sctransform_0.4.3           utf8_1.2.6                 
 [11] R6_2.6.1                    lazyeval_0.2.2             
 [13] uwot_0.2.4                  withr_3.0.2                
 [15] gridExtra_2.3               progressr_0.19.0           
 [17] quantreg_6.1                cli_3.6.5                  
 [19] spatstat.explore_3.8-0      fastDummies_1.7.5          
 [21] iNEXT_3.0.2                 scatterpie_0.2.6           
 [23] labeling_0.4.3              S7_0.2.1                   
 [25] spatstat.data_3.1-9         ggridges_0.5.7             
 [27] pbapply_1.7-4               mixsqp_0.3-54              
 [29] systemfonts_1.3.2           yulab.utils_0.2.4          
 [31] gson_0.1.0                  DOSE_4.4.0                 
 [33] R.utils_2.13.0              parallelly_1.46.1          
 [35] susieR_0.14.2               rstudioapi_0.18.0          
 [37] RSQLite_2.4.6               gridGraphics_0.5-1         
 [39] ica_1.0-3                   spatstat.random_3.4-5      
 [41] zip_2.3.3                   GO.db_3.22.0               
 [43] Matrix_1.7-4                abind_1.4-8                
 [45] R.methodsS3_1.8.2           lifecycle_1.0.5            
 [47] snakecase_0.11.1            SummarizedExperiment_1.40.0
 [49] SparseArray_1.10.10         qvalue_2.42.0              
 [51] Rtsne_0.17                  grid_4.5.1                 
 [53] blob_1.3.0                  promises_1.5.0             
 [55] crayon_1.5.3                miniUI_0.1.2               
 [57] ggtangle_0.1.1              lattice_0.22-9             
 [59] cowplot_1.2.0               KEGGREST_1.50.0            
 [61] pillar_1.11.1               GenomicRanges_1.62.1       
 [63] fgsea_1.36.2                rjson_0.2.23               
 [65] future.apply_1.20.2         codetools_0.2-20           
 [67] fastmatch_1.1-8             glue_1.8.0                 
 [69] ggiraph_0.9.6               ggfun_0.2.0                
 [71] spatstat.univar_3.1-7       fontLiberation_0.1.0       
 [73] data.table_1.18.2.1         vctrs_0.7.1                
 [75] png_0.1-9                   treeio_1.34.0              
 [77] spam_2.11-3                 cellranger_1.1.0           
 [79] gtable_0.3.6                cachem_1.1.0               
 [81] S4Arrays_1.10.1             mime_0.13                  
 [83] tidygraph_1.3.1             Seqinfo_1.0.0              
 [85] survival_3.8-6              SingleCellExperiment_1.32.0
 [87] iterators_1.0.14            fitdistrplus_1.2-6         
 [89] ROCR_1.0-12                 nlme_3.1-168               
 [91] ggtree_4.0.5                bit64_4.6.0-1              
 [93] fontquiver_0.2.1            RcppAnnoy_0.0.23           
 [95] irlba_2.3.7                 KernSmooth_2.23-26         
 [97] otel_0.2.0                  DBI_1.3.0                  
 [99] tidyselect_1.2.1            bit_4.6.0                  
[101] compiler_4.5.1              rvest_1.0.5                
[103] SparseM_1.84-2              xml2_1.5.2                 
[105] ggdendro_0.2.0              DelayedArray_0.36.1        
[107] fontBitstreamVera_0.1.1     plotly_4.12.0              
[109] lmtest_0.9-40               rappdirs_0.3.4             
[111] digest_0.6.39               goftest_1.2-3              
[113] spatstat.utils_3.2-2        XVector_0.50.0             
[115] turboGliph_0.99.2           htmltools_0.5.9            
[117] pkgconfig_2.0.3             MatrixGenerics_1.22.0      
[119] fastmap_1.2.0               htmlwidgets_1.6.4          
[121] shiny_1.13.0                immApex_1.4.3              
[123] farver_2.1.2                zoo_1.8-15                 
[125] jsonlite_2.0.0              BiocParallel_1.44.0        
[127] GOSemSim_2.36.0             R.oo_1.27.1                
[129] magrittr_2.0.4              ggplotify_0.1.3            
[131] dotCall64_1.2               patchwork_1.3.2            
[133] evmix_2.12                  Rcpp_1.1.1                 
[135] ape_5.8-1                   ggnewscale_0.5.2           
[137] viridis_0.6.5               gdtools_0.5.0              
[139] reticulate_1.45.0           stringi_1.8.7              
[141] ggraph_2.2.2                ggalluvial_0.12.6          
[143] MASS_7.3-65                 plyr_1.8.9                 
[145] parallel_4.5.1              listenv_0.10.1             
[147] deldir_2.0-4                graphlayouts_1.2.3         
[149] Biostrings_2.78.0           splines_4.5.1              
[151] hash_2.2.6.4                tensor_1.5.1               
[153] igraph_2.2.2                spatstat.geom_3.7-3        
[155] RcppHNSW_0.6.0              reshape2_1.4.5             
[157] foreach_1.5.2               tweenr_2.0.3               
[159] httpuv_1.6.17               MatrixModels_0.5-4         
[161] RANN_2.6.2                  tidyr_1.3.2                
[163] polyclip_1.10-7             future_1.70.0              
[165] reshape_0.8.10              scattermore_1.2            
[167] ggforce_0.5.0               xtable_1.8-8               
[169] RSpectra_0.16-2             tidytree_0.4.7             
[171] tidydr_0.0.6                later_1.4.8                
[173] viridisLite_0.4.3           gsl_2.1-9                  
[175] snow_0.4-4                  tibble_3.3.1               
[177] aplot_0.2.9                 memoise_2.0.1              
[179] cluster_2.1.8.2             timechange_0.4.0           
[181] globals_0.19.1
