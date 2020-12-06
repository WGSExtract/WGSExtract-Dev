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
# This script works off the current directory or a parameter list.  
#   * If parameter list is one or more files, it works off that.
#   * If a single directory specified, it works off that. 
#   * If nothing specified, then it uses the current directory
# If a directory specified, it tries to find all Reference Genomes (actually, FASTA files) using common extensions.
# It uses each identified FASTA as the base.  
# It can be run a second time and will simply do any needed updates. So if errors fixed, will simply restart.
#
#  NOTE: We have commented out the BWA Index command due to the extensive resource usage (CPU and file space)
#

# Number of Processors to use when available
np=16

# Todo Programs used; maybe parameterize for platform so $PATH not required?
# gunzip, unzip, bunzip2, 7z, bgzip, htsfile; cut, awk, ...

# Figure out what we have to do in the current directory; or if explicitly set in the parameter list
shopt -s nullglob
if [ $# -eq 1 ] && [ -d $1 ]; then
  file_list=($1/*.{fa,fna,fasta,7z,gz,zip,bz,bz2})
else
  file_list=($*)
fi
if [ $# -eq 0 ] || [ ${#file_list[@]} -eq 0 ]; then
  echo "Usage: $0 [dir]  or  $0 [file(s)]"
  echo "Use $0 . to process the current directory"
  echo "Directory searches match any FASTA file with extension: fa, fna, fasta, gz, zip, bz, bz2, or 7z"
  exit
fi
# DEBUG echo "${#file_list[@]} parameters: ${file_list[@]}"

echo "-----------------------------------------------------------------------------"
echo "Checking / correcting files for being FASTA, correct compression, and having indices"
echo "-----------------------------------------------------------------------------"

LANG=POSIX		# Needed for sort command
for file in "${file_list[@]}"; do

  # Check if even a file
  if [ ! -f "$file" ]; then
    [[ $# -gt 2 || $# -eq 1 && ! -d $1 ]] && echo "$file: ***WARNING: Skipping non-file"
    continue 	# Continue silently if processing a directory using pattern match
  fi

  # Check if file content and extension correct
  ext=${file##*.}
  HTS=`htsfile "$file"`
  if [[ "$ext" =~ ^(zip|7z|bz|bz2)$ ]]; then
    # htsfile cannot see in (to tell if FASTA) for .zip, .7z nor .bz2 archive formats (but does recognize bzip2 format)
    # So go ahead if those and check if a FASTA after recompressing; might recompress some non-bioinformatic files but ...
    # that is a risk the user takes if they do a general directory specification
    echo "$file: ***WARNING: Cannot look into .zip/.bz/.bz2/.7z files; so recompressing to peak inside."
  elif [ grep -v "FASTA" <<< "$HTS" &> /dev/null ]; then
    echo "$file: ***WARNING: Skipping a non-FASTA file"
    continue
  fi  # simply fall through as next check determines if need to recompress or not

  if grep "FASTA BGZF" <<< "$HTS" &> /dev/null ; then
    # echo "$file: Already a BGZF compressed FASTA file!"   # simply keep silent
    filen="$file"
  else
    # OK, if here, is (likely) a FASTA but is not BGZF compressed
    echo "$file: Fixing Reference Genome File compression"
    case "$ext" in
      gz)		# BGZF already checked above; so must be gzipped
  	{ gunzip -f "$file" ; bgzip -f -i -@ $np "${file%.*}" ; } && filen="$file"
  	# This is the only format where you end up with the same file extension; they should have used .bgz but ...
  	;;

      zip)		# Likely zip'ped 
	{ unzip -p "$file" | bgzip -cf@ $np > "${file%.*}.gz" ; } && rm -f "$file" && filen="${file%.*}.gz"
  	;;

      7z)		# A popular, very high compression format not normally seen in Bioinformatics
	{ 7z e -mmt$np -so "$file" | bgzip -cf@ $np > "${file%.*}.gz" ; } && rm -f "$file" && filen="${file%.*}.gz"
  	;;

      bz | bz2)		# A popular Linux/Unix/MacOSX high compression format sometimes seen in Bioinformatics
  	{ bunzip2 -cf "$file" | bgzip -cf -@ $np > "${file%.*}.gz" ;} && rm -f "$file" && filen="${file%.*}.gz"
  	;;

      fa | fna | fasta)	# Appears likely not compressed; but will check
  	if grep -- "-compressed" <<< "$HTS" &> /dev/null ; then
	  echo "$HTS\n***WARNING: Is a FASTA file and compressed. But missing proper compressed file name extension. Skipping."
	  continue
	fi
  	bgzip -f -i -@ 16 "$file"
	filen="$file.gz"
  	;;
      *)		# Well that is embarrising
	echo "$file: ***WARNING: Unknown file extension (internal error)"
	;;
    esac
    # Above does not handle .zip / .7z archive where multiple files (i.e. like tar+gzip combo from original pkzip)
    #  Singular, internal file name must the same as archive without extension (e.g. file.zip wth one entry named file)

    # Check again if matched file with extension is a valid, compressed FASTA file just to make sure
    [ ! -f "$filen" ] && echo "$file: ***WARNING: Compression conversion failed (missing file). Skipping" && continue
    HTS=`htsfile "$filen"`
    { grep -v "FASTA BGZF" <<< "$HTS" &> /dev/null ; } && echo "$file: ***WARNING: Compression conversion failed (not BGZF). Skipping" && continue
  fi
  # By this point, Ref Genome file ends in .gz, is FASTA and is BGZF compressed. Ready to create supporting indices, etc.
  #  Does allow explicitely named file through that is FASTA and BGZF compressed but without the extension. Some programs
  #  will not handle that but who are we to judge.

  # Now to (re)create any missing or old indices (using a possibly new filename $filen if compression changed)
  [ "$filen" -nt "$filen.dict" ] && echo "$filen: Creating DICTionary"         && samtools dict "$filen" -o "$filen.dict" 
  [ "$filen" -nt "$filen.fai" ]  && echo "$filen: Creating FAI Index"          && samtools faidx "$filen"
  [ "$filen" -nt "$filen.gzi" ]  && echo "$filen: Creating GZI Index"          && bgzip -r "$filen"   # Also samtools index works

  # Those were the quick and easy ones; now the BWA Index which adds 5.5 GB of files and takes near an hour
  #   BWA Index creates .bwt (30 min, 3 GB), .pac (800MB), .ann, .amb, and .sa (10 min, 1.5GB); with 900M original, 6.4GB total!
  # [ "$filen" -nt "$filen.bwt" ]  && echo "$filen: Creating BWA Indices: 45 min, 5.5 GB" && bwa index "$filen"
  # Todo add parameter to scrip to optionally turn on BWA index capability

  # Finally, lets collect all the special WGS Extract desired stats to understand the FASTA content; all from DICT file
  [ "$filen" -nt "$filen.wgse" ] && echo "$filen: Creating WGS Extract Index"  && {
    # Throw out first line, take only columns 2 to 4, and upcase everything. Store temporarily in output file name
    tail -n +2 "$filen.dict" | cut -f2-4 | awk '{print toupper($0)}' > "$filen.wgse"

    snct=`wc -l "$filen.wgse" | cut -d' ' -f1`				 # True SN count
    md5b=`cut -f1-2 < "$filen.wgse" | sort -d | md5sum | cut -d' ' -f1`  # MD5Sum BAM   (SN, LN)
    md5c=`cat         "$filen.wgse" | sort -d | md5sum | cut -d' ' -f1`  # MD5Sum CRAM  (SN, LN, MD)
    md5f=`cut -f2-3 < "$filen.wgse" | sort -d | md5sum | cut -d' ' -f1`  # MD5Sum FASTA (LN, MD)

    # Model 36, 37 and 38 chromosome length counts (one long pattern line)
    chrs="LN:247249719\|LN:249250621\|LN:248956422\|LN:242951149\|LN:243199373\|LN:242193529\|LN:199501827\|LN:198022430\|LN:198295559\|LN:191273063\|LN:191154276\|LN:190214555\|LN:180857866\|LN:180915260\|LN:181538259\|LN:170899992\|LN:171115067\|LN:170805979\|LN:158821424\|LN:159138663\|LN:159345973\|LN:146274826\|LN:146364022\|LN:145138636\|LN:140273252\|LN:141213431\|LN:138394717\|LN:135374737\|LN:135534747\|LN:133797422\|LN:134452384\|LN:135006516\|LN:135086622\|LN:132349534\|LN:133851895\|LN:133275309\|LN:114142980\|LN:115169878\|LN:114364328\|LN:106368585\|LN:107349540\|LN:107043718\|LN:100338915\|LN:102531392\|LN:101991189\|LN:88827254\|LN:90354753\|LN:90338345\|LN:78774742\|LN:81195210\|LN:83257441\|LN:76117153\|LN:78077248\|LN:80373285\|LN:63811651\|LN:59128983\|LN:58617616\|LN:62435964\|LN:63025520\|LN:64444167\|LN:46944323\|LN:48129895\|LN:46709983\|LN:49691432\|LN:51304566\|LN:50818468\|LN:154913754\|LN:155270560\|LN:156040895\|LN:57772954\|LN:59373566\|LN:57227415"
    chrsn="SN:CHR[123456789XY]\>\|SN:CHR1[0123456789]\>\|SN:CHR2[012]\>\|SN:[123456789XY]\>\|SN:1[0123456789]\>\|SN:2[012]\>\|CM0006\|SN:NC_0000[012][0123456789]"    # Some alt contigs same length as chromosomes


    # Capture Primary (Chromosomes) and Mitochondria key stats
    primary=`grep -w -e "$chrs" "$filen.wgse" | grep -e "$chrsn"`	# Length is not enough; need to then filter by SN
    pcnt=`wc -l <<< $primary | cut -d' ' -f1` 
    if [ $pcnt -ne "24" ]; then
      ycnt=`grep -w -e "LN:59373566\|LN:57227415" <<< $primary | wc -l | cut -d' ' -f1`
      if [ $ycnt -ne "1" ]; then
        # M and Y are only duplicates found so far; so give a special error if found to be Y
        echo "$filen: ***ERROR: 1 expected, $ycnt Y chromosome entries found in ref model"
        md5p="******** ERROR:Y $pcnt!=24 *********"
      else
        echo "$filen: ***ERROR: 24 expected, $pcnt chromosomes found in ref model" 
        md5p="******** ERROR:P $pcnt!=24 **********"
      fi
    else
      md5p=`cut -f2-3 <<<$primary | sort -d | md5sum | cut -d' ' -f1`  # MD5Sum Chrs-only (LN, MD)
    fi

    chrM=`grep -w -e "LN:16569\|LN:16571" < "$filen.wgse"`	# Keep SN and MD columns to help capture detail
    mcnt=`wc -l <<< $chrM | cut -d' ' -f1`
    if [ $mcnt -ne "1" ]; then
      echo "$filen: ***ERROR: 1 expected, $mcnt microchonrdrial entries found in ref model"
      chrM="SN:ERROR LN:MULTI M5:******** ERROR:M $mcnt!=1 **********"
    fi

    # Save all the values as single line, space separated text file 
    echo "$snct $md5b $md5c $md5f $md5p $chrM \"$filen\"" > "$filen.wgse"

    # Final error check -- in case previous uncaught error left intermediate .wgse file or multi-line stats (mtDNA)
    wcnt=`wc -l "$filen.wgse" | cut -d' ' -f1`
    if [ $wcnt -ne "1" ]; then
      echo "$filen: ***ERROR: Failed generating final WGSE stats file"
      rm -f "$filen.wgse"
    fi
  }

done
