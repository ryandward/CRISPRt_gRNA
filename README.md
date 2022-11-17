# CRISPRt_gRNA

Run as usual, then

```
awk 'NR == 1 {print; next} $2 <= -95 && $9 == "sense" && $10 == 39 {print $0 | "sort -k1,1V -k2,2nr"}' sgRNA_outputs.tsv | 
awk 'NR == 1 {print; next} {guides[$1]++} guides[$1] <= 10 {print}' > top_ten_sgRNA_outputs.tsv #ranked
```
