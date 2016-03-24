#!/bin/bash
# This tests every aakbar command.
aakbar
aakbar init_config_file .
aakbar show_context_object
for dir in `cat dirs.lst`; do
  label=`cat "${dir}/scientific_name.txt"`
  identifier=`cat "${dir}/abbreviation.txt"`
  aakbar define_set ${identifier} ${dir}
  aakbar label_set ${identifier} "${label}"
done
aakbar define_summary intersect "Intersection"
aakbar set_plot_type
aakbar set_plot_type png
aakbar set_letterfreq_window 12
aakbar set_simplicity_object
aakbar show_configuration
#
# loop over all simplicity functions except null
#
simplicity_functions=`FS=":" aakbar set_simplicity_object | grep ":" | awk '{print substr($1,0,length($1)-1)}' | grep -v null`
for func in $simplicity_functions; do
  aakbar set_simplicity_object $func
  aakbar demo_simplicity
  for cut in 3 2; do
     aakbar --progress peptide_simplicity_mask --plot --cutoff ${cut} protein.faa protein-${func}-${cut} all
     aakbar --progress calculate_peptide_terms -k 10 protein-${func}-${cut}.faa protein-${func}-${cut}_k-10 all
     aakbar intersect_peptide_terms protein-${func}-${cut}_k-10 all
done
