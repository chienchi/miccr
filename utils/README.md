## Utility Scripts to generate the classification result. 

### Usage:

Two steps:
 
1. Add taxonomy lineage to the result tsv file

```
./add_lineage.py dbPath test/TEST.ctg.tsv  > TEST.ctg.tsv.lineage
```

2. Run the plot R script with coverage information table file.

```
Rscript classification_plot.R TEST.ctg.tsv.lineage out_prefix coverage.table
```


The coverage.table file is a tab-delimited text file with 6 columns or more:
ex:

| ID  | Length|  GC%  | Avg_fold | Fold_std | Base_Coverage% |
|-----|-------|-------|----------|----------|----------------|
|seq00| 60366 | 34.86 |   21.71  |   5.60   |  100.0000      |
|seq01| 36715 | 59.60 |   63.79  |  32.01   |  100.0000      |
|seq02| 29799 | 67.30 |   43.94  |  18.12   |   99.9933      | 



#### Alternative 

Run the plot R script without coverage information table file.

```
Rscript classification_plot.R TEST.ctg.tsv.lineage out_prefix
```


