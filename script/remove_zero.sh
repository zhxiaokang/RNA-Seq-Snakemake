#remove all-zero lines from a paired countFile (control and one treated group)

paste $1 $2 | awk 'NR==1{print $0} NR>1{s=0; for (i=2; i<=NF; i++) s+=$i} s>0{print $0}'
