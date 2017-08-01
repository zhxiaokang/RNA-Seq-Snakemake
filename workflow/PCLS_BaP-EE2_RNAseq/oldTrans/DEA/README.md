The files named {condition}_VS_control are the raw counts of each pair (a condition versus the control), after removing the all-zero genes.

The files named dea_{condition} are the results of differential expression analysis (DEA) from DESeq2 based on paired sample test.

The files named deg_{condition} only include the differentially expressed genes (DEG) with padj < 0.05
