#!/bin/bash
aakbar init_config_file .
for dir in `cat dirs.lst`; do
  name=`cat "${dir}/common_name.txt"`
  aakbar define_set ${dir} ${dir}
  aakbar label_set ${dir} "${name}"
done
aakbar set_complexity_function --function runlength
#aakbar --progress calculate_peptide_kmers -k 10 protein.faa protein-10mer all
#aakbar --progress filter_peptide_kmers --runlength 2 protein-10mer_terms.tsv protein-10mer_terms-cutoff-2 all
aakbar define_summary intersect "Intersection"
aakbar merge_peptide_kmers protein-10mer_terms-cutoff-2_terms.tsv all
