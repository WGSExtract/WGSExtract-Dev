#!/bin/sh
#
# Script to compare Reference Genomes.
# 
# The script naively presumes any file in the folder with the correct extension is likely a FASTA file. But will skip if not.
#
# Part of the Reference Genome package in WGS Extract (https://wgsextract.github.io/)
#
# Copyright (c) 2020 Randy Harr
#

# Figure out what we have to do in the current directory; or if explicitly set in the parameter list
shopt -s nullglob
if [ $# -ne 2 ]; then
  echo "Usage: $0 ref_gen_1 ref_gen_2"
  exit
fi

LANG=POSIX		# Needed for sort command
# Throw out first line, take only columns 2 to 4, and upcase everything. Store temporarily in output file name
tail -n +2 "$1.dict" | cut -f2-4 | awk '{print toupper($0)}' | sort -d > "$1.tmp"
tail -n +2 "$2.dict" | cut -f2-4 | awk '{print toupper($0)}' | sort -d > "$2.tmp"

# Model 36, 37 and 38 chromosome length counts (one long pattern line)
chrs="LN:247249719\|LN:249250621\|LN:248956422\|LN:242951149\|LN:243199373\|LN:242193529\|LN:199501827\|LN:198022430\|LN:198295559\|LN:191273063\|LN:191154276\|LN:190214555\|LN:180857866\|LN:180915260\|LN:181538259\|LN:170899992\|LN:171115067\|LN:170805979\|LN:158821424\|LN:159138663\|LN:159345973\|LN:146274826\|LN:146364022\|LN:145138636\|LN:140273252\|LN:141213431\|LN:138394717\|LN:135374737\|LN:135534747\|LN:133797422\|LN:134452384\|LN:135006516\|LN:135086622\|LN:132349534\|LN:133851895\|LN:133275309\|LN:114142980\|LN:115169878\|LN:114364328\|LN:106368585\|LN:107349540\|LN:107043718\|LN:100338915\|LN:102531392\|LN:101991189\|LN:88827254\|LN:90354753\|LN:90338345\|LN:78774742\|LN:81195210\|LN:83257441\|LN:76117153\|LN:78077248\|LN:80373285\|LN:63811651\|LN:59128983\|LN:58617616\|LN:62435964\|LN:63025520\|LN:64444167\|LN:46944323\|LN:48129895\|LN:46709983\|LN:49691432\|LN:51304566\|LN:50818468\|LN:154913754\|LN:155270560\|LN:156040895\|LN:57772954\|LN:59373566\|LN:57227415\|LN:16569\|LN:16571"

chrsn="SN:CHR[123456789XYM]\>\|SN:CHR1[0123456789]\>\|SN:CHR2[012]\>\|SN:[123456789XY]\>\|SN:1[0123456789]\>\|SN:2[012]\>\|SN:CHRMT\>\|SN:MT\>\|SN:CM0006\|SN:J01415\|SN:NC_0000[012][0123456789]\|SN:NC_012920"    # Some alt contigs same length as chromosomes

# Capture Primary (Chromosomes) and Mitochondria key stats
grep -w -e "$chrs" "$1.tmp" | grep -e "$chrsn" > "$1.tmp1"	# Length is not enough; need to then filter by SN
grep -w -e "$chrs" "$2.tmp" | grep -e "$chrsn" > "$2.tmp1"	# Length is not enough; need to then filter by SN


# Capture Primary (Chromosomes) and Mitochondria key stats
cut -f2-3 < "$1.tmp1" > "$1.tmp2"	# Length is not enough; need to then filter by SN
cut -f2-3 < "$2.tmp1" > "$2.tmp2"	# Length is not enough; need to then filter by SN

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Difference between $1 and $2 w/out SNs"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
diff "$1.tmp2" "$2.tmp2"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Difference between $1 and $2 with SNs"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
diff "$1.tmp1" "$2.tmp1"

#rm $1.tmp $2.tmp $1.tmp1 $2.tmp1 $1.tmp2 $2.tmp2

