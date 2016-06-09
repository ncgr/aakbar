#!/bin/bash
#
# Calculate signatures for a set of 10 Streptococci
#
download_failed="Downloader returned an error, fix the problem."
#
# Download data files for signature calculation.
#
./genbank_downloader.sh -q -c strep bacteria/Streptococcus_pneumoniae/reference/GCA_000007045.1_ASM704v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.oral bacteria/Streptococcus_oralis/latest_assembly_versions/GCA_000148565.2_ASM14856v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.infant bacteria/Streptococcus_infantis/latest_assembly_versions/GCA_000148975.2_ASM14897v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.gordon bacteria/Streptococcus_gordonii/representative/GCA_000017005.1_ASM1700v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.suis bacteria/Streptococcus_suis/reference/GCA_000026745.1_ASM2674v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.mutans bacteria/Streptococcus_mutans/reference/GCA_000007465.2_ASM746v2 
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.thermo bacteria/Streptococcus_thermophilus/representative/GCA_000011825.1_ASM1182v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.horse bacteria/Streptococcus_equinus/latest_assembly_versions/GCA_000146405.1_ASM14640v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.pyogenes bacteria/Streptococcus_pyogenes/reference/GCA_000006785.2_ASM678v2
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -q -c S.pig bacteria/Streptococcus_porcinus/latest_assembly_versions/GCA_000187955.3_ASM18795v2
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
# test genome
./genbank_downloader.sh -q -c S.dog -t genome bacteria/Streptococcus_canis/latest_assembly_versions/GCA_000268305.2_ASM26830v2
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
#
# Now calculate signatures.
#
./calculate_signatures.sh -o strep10 -p png Streptococcus_pneumoniae Streptococcus_equinus Streptococcus_gordonii Streptococcus_infantis Streptococcus_mutans Streptococcus_oralis Streptococcus_porcinus Streptococcus_pyogenes Streptococcus_suis Streptococcus_thermophilus
#
# Test the signatures on a genome not included in the set above.
#
./split.sh Streptococcus_canis/genome.fna 200
aakbar define_set S.dog Streptococcus_canis
aakbar label_set S.dog "Streptococcus canis"
aakbar --progress search_peptide_occurrances --nucleotides genome-200-bp_reads.fna strep10 S.dog
