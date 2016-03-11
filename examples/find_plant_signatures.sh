#!/bin/bash
# This tests every aakbar command.
aakbar
aakbar init_config_file .
aakbar show_context_object
for dir in `cat dirs.lst`; do
  name=`cat "${dir}/scientific_name.txt"`
  aakbar define_set ${dir} ${dir}
  aakbar label_set ${dir} "${name}"
done
aakbar define_summary intersect "Intersection"
aakbar set_plot_type
aakbar set_plot_type png
aakbar set_simplicity_object
aakbar set_simplicity_object runlength
aakbar show_configuration
aakbar demo_simplicity
aakbar --progress peptide_simplicity_mask --plot --cutoff 3 protein.faa protein-runlength-3 all
aakbar --progress calculate_peptide_terms -k 10 protein-runlength-3.faa protein-runlength-3-k-10 all
aakbar intersect_peptide_terms protein-runlength-3-k-10.tsv all
