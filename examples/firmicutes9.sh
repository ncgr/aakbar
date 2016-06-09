#!/bin/bash
#
#  Analyze signatures for a set of 9 firmicutes
#
download_failed="Downloader returned an error, fix the problem."
#
# Download data files for signature calculation.
#
./genbank_downloader.sh -c aerovir bacteria/Aerococcus_viridans/latest_assembly_versions/GCA_000178435.1_ASM17843v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -c b.sub bacteria/Bacillus_subtilis/reference/GCA_000009045.1_ASM904v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -c botox bacteria/Clostridium_botulinum/reference/GCA_000063585.1_ASM6358v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -c finegold bacteria/Finegoldia_magna/representative/GCA_000010185.1_ASM1018v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -c lacto bacteria/Lactobacillus_brevis/representative/GCA_000014465.1_ASM1446v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -c listeria bacteria/Listeria_monocytogenes/reference/GCA_000196035.1_ASM19603v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -c almost bacteria/Paenibacillus_alvei/representative/GCA_000293805.1_ASM29380v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -c strep bacteria/Streptococcus_pneumoniae/reference/GCA_000007045.1_ASM704v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
./genbank_downloader.sh -c thermo bacteria/Thermoanaerobacterium_saccharolyticum/latest_assembly_versions/GCA_000307585.1_ASM30758v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
# test genome
./genbank_downloader.sh -q -c kimchi -t genome bacteria/Leuconostoc_kimchii/representative/GCA_000092505.1_ASM9250v1
if [ $? -ne 0 ]; then echo "$download_failed"; exit 1; fi
#
# Now calculate signatures.
#
./calculate_signatures.sh -o firmicutes9 -p png Bacillus_subtilis Aerococcus_viridans Clostridium_botulinum Finegoldia_magna Lactobacillus_brevis Listeria_monocytogenes Paenibacillus_alvei Streptococcus_pneumoniae Thermoanaerobacterium_saccharolyticum
#
# Test the signatures on a genome not included in the set above,
# taking 200-bp chunks of the genome as representative.
#
./split.sh Leuconostoc_kimchii/genome.fna 200
aakbar define_set kimchi Leuconostoc_kimchii
aakbar label_set kimchi "Leuconostoc kimchii"
aakbar --progress search_peptide_occurrances --nucleotides genome-200-bp_reads.fna firmicutes9 kimchi
