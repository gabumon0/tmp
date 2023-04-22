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
/data/COLOR/sclinker_test/script/run_process_sclinker_output_step3.sh
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
