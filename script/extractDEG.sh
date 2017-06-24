# extract the genes with padj < 0.05
awk -F ',' '$NF<0.05{print $0}' $1