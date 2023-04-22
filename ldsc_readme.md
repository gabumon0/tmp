##### step1
```
####to generate celltype and disease program
/data/COLOR/sclinker_test/script/generateGenePrograms_step1.sh
```

* input
  * /data/mydata6T1/users/R/R_packges.bak/expr.h5ad, h5ad file
* output
  * information from different gene analysis in sc.tl.rank_genes_groups function, results like below:
  * /data/COLOR/sclinker_test/data/test_celltype_genescores.csv
  * /data/COLOR/sclinker_test/data/test_disease_genescores.csv
  * /data/COLOR/sclinker_test/data/test_celltype_logfold.csv
  * /data/COLOR/sclinker_test/data/test_disease_logfold.csv
  * /data/COLOR/sclinker_test/data/test_celltype_pval.csv
  * /data/COLOR/sclinker_test/data/test_disease_pval.csv
  * /data/COLOR/sclinker_test/data/test_celltype_score.csv
  * /data/COLOR/sclinker_test/data/test_disease_score.csv

#### step2
```
####to generate cellular process
/data/COLOR/sclinker_test/script/generateNMFmodule_step2.py
```

* input 
  * /data/mydata6T1/users/R/R_packges.bak/expr.h5ad, h5ad file
* output
  * NMF dimensionality reduction
  * /data/COLOR/sclinker_test/data/nmf/test_cellprograms.csv
  * /data/COLOR/sclinker_test/data/nmf/test_correlation_cellprograms.csv
  * /data/COLOR/sclinker_test/data/nmf/test_geneprograms.csv

#### step3
```
####to generate a Gene Program
/data/COLOR/sclinker_test/script/run_process_sclinker_output_step3.sh (process_sclinker_output.R)
```

* input
  * /data/COLOR/sclinker_test/data/test_celltype_score.csv

* output
  * one file each celltype
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/ALL.txt
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal.txt
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/endothelial.mural.txt
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/interneurons.txt
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/microglia.txt
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/oligodendrocytes.txt
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/pyramidal.CA1.txt
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/pyramidal.SS.txt

#### step4
```
####Gene sets to bedgraph format files using region-gene linking
/data/COLOR/sclinker_test/script/run_geneset_to_bed_sclinker_step4.sh (geneset_to_bed_sclinker.R)
```

* input
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal.txt
  * /data/COLOR/sclinker_test/processed_data/Dey_Enhancer_MasterReg/data/
* output
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/ABC_Road_ALL.bed
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/100kb.bed

#### step5
```
####cleaning the bedgraph format files
/data/COLOR/sclinker_test/script/run_clean_bed_step5.sh (clean_bed.sh)
```

* input
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/ABC_Road_ALL.bed
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/100kb.bed
* output
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/ABC_Road_ALL.bed (has cleaned)
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/100kb.bed (has cleaned)

#### step6
```
####from bedgraph to annotation files
/data/COLOR/sclinker_test/script/run_bed_to_annot_step6.sh (bed_to_annot.sh, bedgraph_to_annot.py)
```

* input
  * /data/COLOR/sclinker_test/processed_data/Dey_Enhancer_MasterReg/data/BIMS
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/ABC_Road_ALL.bed & 100kb.bed
* output
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype_annot/astrocytes_ependymal/100kb/*annot.gz
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype_annot/astrocytes_ependymal/ABC_Road_ALL/*annot.gz

#### step7
```
####create LD-scores from a reference panel
/data/COLOR/sclinker_test/script/create_ldscores.sh
```

* input 
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/astrocytes_ependymal*.annot.gz  (copy from /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype_annot/astrocytes_ependymal/ABC_Road_ALL/*annot.gz)
* output
  * /data/COLOR/sclinker_test/data/sclinker_genescores/test_celltype/astrocytes_ependymal/astrocytes_ependymal*l2.ldscore.gz

#### step8
```
####run regression model in S-LDSC 
/data/COLOR/sclinker_test/script/run_ldsc_reg.sh   (--ref-ld-chr参数需要调整)
```

* input

* output
  * /data/COLOR/sclinker_test/data/sclinker_out/baselineLD_v2.2/astrocytes_ependymal/*results
