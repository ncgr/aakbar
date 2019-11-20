#!/bin/bash
set -e
error_exit() {
   >&2 echo "ERROR--unexpected exit from ${BASH_SOURCE} script at line:"
   >&2 echo "   $BASH_COMMAND"
   >&2 echo "Directory is \"${PWD}\"."
}
trap error_exit EXIT
docstring="
 Use aakbar to calculate a set of signatures from a list of directories containing
 identically-names FASTA files of called protein sequences.

Usage:
       calculate_signatures.sh [-k K] [-f FUNCTION] [-c CUT] [-i INNAME] [-o OUTNAME] [-s SCORE] [-p PLOT_TYPE] DIRLIST

         where
           K is the k-mer size in amino-acid residues (default is 10)
    FUNCTION is the simplicity function to use (default is letterfreq12)
         CUT is the simplicity cutoff value (default is 6)
      INNAME is the name of the input FASTA file (default is \"protein.faa\")
     OUTNAME is the name of the resulting signature files (default is \"signatures\")
       SCORE is the simplicity cutoff score (default is 0.1)
   PLOT_TYPE is the output plot type (default is PNG)
     DIRLIST is the list of directories containing input files.
"
# default values
k="10"
inname="protein.faa"
func="letterfreq12"
cut="5"
score="0.1"
plot_type="png"
outname="signatures"
# handle options
while getopts "hk:f:c:i:o:s:p:" opt; do
  case $opt in
     h)
       echo "${docstring}"
       trap - EXIT
       exit 0
       ;;
     k)
       k="${OPTARG}"
       ;;
     f)
      func="${OPTARG}"
      ;;
     c)
      cut="${OPTARG}"
      ;;
     i)
       inname="${OPTARG}"
       ;;
     o)
       outname="${OPTARG}"
       ;;
     s)
      score="${OPTARG}"
      ;;
     p)
      plot_type="${OPTARG}"
      ;;
    \?)
       echo "Invalid option -$OPTARG" >&2
       echo "${docstring}"
       trap - EXIT
       exit 1
       ;;
     :)
       echo "ERROR-option -$OPTARG requires an argument." >&2
       trap - EXIT
       exit 1
       ;;
  esac
done
# handle positional arguments
shift $(($OPTIND - 1))
ndirs=${#@}
if [ $ndirs -eq 0 ] ; then
  echo "Must specify a list of directories." >&2
  echo "${docstring}"
  trap - EXIT
  exit 1
fi
echo "Signatures will be calculated from ${ndirs} genomes."
# echo arguments
inext="${inname##*.}"
instem="${inname%.*}"
echo "Input file name is \"${instem}.${inext}\"."
echo "k-mer length will be ${k}."
echo "Simplicity object will be ${func} with cutoff of ${cut}."
echo "Simplicity filter score cutoff will be ${score}."
echo "Plot type will be ${plot_type}."
echo "Output signature directory and name will be \"${outname}\"."
#
echo
echo "Initializing aakbar configuration:"
aakbar init-config-file .
for dir in $@; do
  while read key value; do
    if [[ "${key}" == "scientific_name" ]]; then
      label=$value
    elif [[ "${key}"  ==  "code" ]]; then
      identifier=$value
    fi
  done <${dir}/metadata.tsv
  if [ "$labels" ]; then
    labels="${labels}+${identifier}"
  else
    labels="${identifier}"
  fi
  aakbar define-set ${identifier} ${dir}
  aakbar label-set ${identifier} "${label}"
  echo 
done
aakbar define-summary ${outname} "${labels}"
if [[ $func == letterfreq* ]]; then
  aakbar set-letterfreq-window ${func##letterfreq}
fi
aakbar set-simplicity-object $func
aakbar set-plot-type $plot_type
aakbar show-config
#
echo
echo "Demo $func simplicity with cutoff of $cut: "
aakbar demo-simplicity -k $k --cutoff $cut
#
echo
echo "Calculate simplicity masks:"
aakbar --progress peptide-simplicity-mask --plot --cutoff $cut ${inname} ${instem}_${func}-${cut} all
#
echo
echo "Calculate peptide terms: "
aakbar --progress calculate-peptide-terms -k ${k} \
  ${instem}_${func}-${cut}.${inext} ${instem}_${func}-${cut}_k-${k} all
#
echo
echo "Compute terms in common across ${ndirs} genomes:"
aakbar intersect-peptide-terms ${instem}_${func}-${cut}_k-${k} all
#
echo
echo "Filter terms with score of ${score}: "
aakbar filter-peptide-terms --cutoff ${score} ${instem}_${func}-${cut}_k-${k} \
       ${outname}
#
echo
echo "Back-search for occurrances of terms:"
aakbar --progress search-peptide-occurrances  ${inname} \
         ${outname} all
#
echo
echo "Done with ${outname} signature calculations"
trap - EXIT
exit 0
