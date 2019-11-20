#!/bin/bash
set -e
error_exit() {
   >&2 echo "ERROR--unexpected exit from ${BASH_SOURCE} script at line:"
   >&2 echo "   $BASH_COMMAND"
   >&2 echo "Directory is \"${PWD}\"."
}
trap error_exit EXIT
#
# Calculate signatures for a set of 10 Streptococci
#
download_failed="Downloader returned an error, fix the problem."
#
# Download data files for signature calculation.
#
echo "Downloading data files from GenBank for signature calculations"
echo "If downloading fails, this script may be restarted without harm"
./genbank_downloader.sh -q -c strep bacteria/Streptococcus_pneumoniae/reference/GCA_000007045.1_ASM704v1
./genbank_downloader.sh -q -c S.oral bacteria/Streptococcus_oralis/latest_assembly_versions/GCA_000148565.2_ASM14856v1
./genbank_downloader.sh -q -c S.infant bacteria/Streptococcus_infantis/latest_assembly_versions/GCA_000148975.2_ASM14897v1
./genbank_downloader.sh -q -c S.gordon bacteria/Streptococcus_gordonii/representative/GCA_000017005.1_ASM1700v1
./genbank_downloader.sh -q -c S.suis bacteria/Streptococcus_suis/reference/GCA_000026745.1_ASM2674v1
./genbank_downloader.sh -q -c S.mutans bacteria/Streptococcus_mutans/reference/GCA_000007465.2_ASM746v2
./genbank_downloader.sh -q -c S.thermo bacteria/Streptococcus_thermophilus/representative/GCA_000253395.1_ASM25339v1
./genbank_downloader.sh -q -c S.horse bacteria/Streptococcus_equinus/latest_assembly_versions/GCA_000146405.1_ASM14640v1
./genbank_downloader.sh -q -c S.pyogenes bacteria/Streptococcus_pyogenes/reference/GCA_000006785.2_ASM678v2
./genbank_downloader.sh -q -c S.aglac bacteria/Streptococcus_agalactiae/reference/GCA_000007265.1_ASM726v1
# test genome
./genbank_downloader.sh -q -c S.pig -t genome bacteria/Streptococcus_porcinus/latest_assembly_versions/GCA_000187955.3_ASM18795v2
#
# Now calculate signatures.
#
./calculate_signatures.sh -o strep10 -p png Streptococcus_pneumoniae Streptococcus_equinus Streptococcus_gordonii Streptococcus_infantis Streptococcus_mutans Streptococcus_oralis Streptococcus_agalactiae Streptococcus_pyogenes Streptococcus_suis Streptococcus_thermophilus
#
# Test the signatures on a genome not included in the set above,
# taking 200-bp chunks of the genome as representative.
#
echo
echo "Searching simulated reads from test genome:"
./split.sh Streptococcus_porcinus/genome.fna 200
aakbar define-set S.pig Streptococcus_porcinus
aakbar label-set S.pig "Streptococcus porcinus"
aakbar --progress search-peptide-occurrances --nucleotides genome-200-bp_reads.fna strep10 S.pig
#
# Do stats on highly-conserved signatures.
#
echo
echo "Doing stats on signatures in test genome:"
aakbar conserved-signature-stats genome-200-bp_reads strep10 S.pig
echo "strep10 example finished"
trap - EXIT
exit 0
