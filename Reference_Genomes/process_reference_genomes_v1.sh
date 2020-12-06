#!/bin/bash
#
# Script to process a folder of Reference Genome files in preparation for use in Bioinformatic tools
# 
# The script naively presumes any file in the folder with the correct extension is likely a FASTA file. But will skip if not.
#
# Part of the Reference Genome package in WGS Extract (https://wgsextract.github.io/)
#
# Copyright (c) 2020 Randy Harr
#

#
# This script works of the current directory and tries to find all Reference Genomes (actually, FASTA files). It uses each
#  identified FASTA as the base.  It can be run a second time and will simply do any needed updates. So if errors fixed,
#  simply restart.
#
#  NOTE: We have commented out the BWA Index command due to the extensive resource usage (CPU and file space)
#

echo "-----------------------------------------------------------------------------"
echo "Checking / correcting files for being FASTA, correct compression, and having indices"
echo "-----------------------------------------------------------------------------"

LANG=POSIX		# Needed for sort command
for file in *.{fa,fna,fasta,7z,gz,zip,bz,bz2} ; do

  # Check for correct file type and compression format of all files first
  [ ! -f "$file" ] && continue   # Not a valid file so simply skip; may be template if no files found for that template
  HTS=`htsfile "$file"`
  # Check if matched file with extension is even a FASTA file or maybe already compressed properly
  if grep -v "FASTA" <<< "$HTS" &> /dev/null ; then
    echo "***WARNING: Skipping a non-FASTA file with an expected FASTA extension: $file"
    continue
  elif grep "FASTA BGZF" <<< "$HTS" &> /dev/null ; then
    # echo "$file: Already a BGZF compressed FASTA file!"   # Maybe simply keep silent
    filen=$file
  else
    # OK, if here, is FASTA but not BGZF compressed
    echo "$file: Fixing Reference Genome File compression"
    ext=${file##*.}
    case "$ext" in
      gz)		# BGZF already checked above; so must be gzipped
  	# echo "Recompressing gzip'ed file using BGZF"
  	gunzip -f "$file" && bgzip -f -@ 16 "${file%.*}"
  	# This is the only format where you will end up with the same file extension; they should have used .bgz but ...
	filen=$file
  	;;

      zip)		# Likely zip'ped
	filen="${file%.*}.gz"
  	;;

      7z)		# A popular, very high compression format not normally seen in Bioinformatics
	filen="${file%.*}.gz"
  	;;

      bz | bz2)		# A popular Linux/Unix/MacOSX high compression format sometimes seen in Bioinformatics
  	# echo "$file: Recompressing BZip'ed file using BGZF"
  	bunzip2 -f "$file" && bgzip -f -@ 16 "${file%.*}"
	filen="${file%.*}.gz"
  	;;

      fn | fna | fasta)	# Appears likely not compressed; but will check
  	if grep -- "-compressed" <<< "$HTS" &> /dev/null ; then
	  echo "$HTS\n***WARNING: Is a FASTA file and compressed. But missing proper compression file name extension. Skipping."
	  continue
	fi
  	# echo "$file: Compressing using BGZF"
  	bgzip -f -@ 16 "$file"
	filen="$file.gz"
  	;;
    esac
  fi
  # All valid files must end in .gz by this point; now to create any missing or old indices (using a possibly new filename $filen

  # Check again if matched file with extension is a valid, compressed FASTA file just to make sure
  [ ! -f "$filen" ] && echo "***WARNING: Compression conversion failed. Skipping: $file" && continue
  HTS=`htsfile "$filen"`
  { grep -v "FASTA BGZF" <<< "$HTS" &> /dev/null ; } && echo "***WARNING: Compression conversion failed. Skipping: $file" && continue

  [ "$filen" -nt "$filen.dict" ] && echo "$filen: Creating DICTionary"         && samtools dict "$filen" -o "$filen.dict" 
  [ "$filen" -nt "$filen.md5c" ] && echo "$filen: Creating header MD5Sum CRAM" && tail -n +2 "$filen.dict" | cut -f2-4 | awk '{print toupper($0)}' | sort -d | md5sum | cut -d' ' -f1 > "$filen.md5c"
  [ "$filen" -nt "$filen.md5b" ] && echo "$filen: Creating header MD5Sum BAM"  && tail -n +2 "$filen.dict" | cut -f2-3 | awk '{print toupper($0)}' | sort -d | md5sum | cut -d' ' -f1 > "$filen.md5b"
  [ "$filen" -nt "$filen.fai" ]  && echo "$filen: Creating FAI Index"          && samtools faidx "$filen"
  [ "$filen" -nt "$filen.gzi" ]  && echo "$filen: Creating GZI Index"          && samtools index "$filen"

  # Those were the quick and easy ones; now the BWA Index which adds 5.5 GB of files and takes near an hour
  #   BWA Index creates .bwt (30 min, 3 GB), .pac (800MB), .ann, .amb, and .sa (10 min, 1.5GB); with 900M original, 6.4GB total!
  # [ "$filen" -nt "$filen.bwt" ]  && echo "$filen: Creating BWA Indices: 45 min, 5.5 GB" && bwa index "$filen"
  # echo "$filen: Delaying creation of the BWA Indices"   # Maybe just be silent about it
done


# Helpful coding reminder since CSH :r:f:h etc are not available!
# FILE="example.tar.gz"
# echo "${FILE%%.*}"
#example
# echo "${FILE%.*}"
#example.tar
# echo "${FILE#*.}"
#tar.gz
# echo "${FILE##*.}"
#gz
