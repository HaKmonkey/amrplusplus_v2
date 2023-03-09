f_line=$(head -n 1 $1);
echo "gene_accession,"$f_line > $2;
tail -n +2 $1 >> $2;