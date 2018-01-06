# extract the genes with padj < 0.1
awk -F ',' 'NR==1{next} $NF<0.15{print $1}' $1 | cut -f 2 -d '"'
