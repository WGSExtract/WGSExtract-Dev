#!/bin/bash
#
# Script to process a folder of Reference Genome files in preparation for use in Bioinformatic tools
# 
# The script naively presumes any file in the folder with the correct extension is a FASTA file.
#
# Part of the Reference Genome package in WGS Extract (https://wgsextract.github.io/)
#
# Copyright (c) 2020-21 Randy Harr
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
# gunzip, unzip, bunzip2, 7z, sort, cut, (g)awk, sed, grep, tail, head, wc, md5sum, rm, bgzip, htsfile, samtools


# Figure out what we have to do in the current directory; or if explicitly set in the parameter list
shopt -s nullglob
if [ $# -eq 1 ] && [ -d $1 ] ; then
  file_list=("$1"/*.{fa,fna,fasta,7z,gz,zip,bz,bz2})
else
  file_list=("$*")
fi
if [ $# -eq 0 ] || [ ${#file_list[@]} -eq 0 ] ; then
  echo "Usage: $0 [dir]  or  $0 [file(s)]"
  echo "Use $0 . to process the current directory"
  echo "Directory searches match any FASTA file with extension: fa, fna, fasta, gz, zip, bz, bz2, or 7z"
  exit
fi
# DEBUG echo "${#file_list[@]} parameters: ${file_list[@]}"

case $OSTYPE in
  msys*|cygwin*|linux*)
    ;;
  darwin*)
    # Need to pick up Macports executables
    PATH="/opt/local/bin:$PATH"
    ;;
esac


echo "-----------------------------------------------------------------------------"
echo "Checking / correcting files for being FASTA, correct compression, and having indices"
echo "-----------------------------------------------------------------------------"

# Setup initial csv summary file header if not yet written; for when initial processing of whole directory (TSV)
if [ ! -f "WGSE.csv" ] ; then
	echo "File	Build	SN CNT	BAM (SN, LN)	CRAM (SN. LN., M5)	Fasta (LN, M5)	Chromo (LN, M5)	Mito SN	Mitp LN	Mito M5	Error" > "WGSE.csv"
fi
	
LANG=POSIX		# Needed for sort command
for file in "${file_list[@]}" ; do

  # Check if even a file
  if [ ! -f "$file" ] ; then
    [[ $# -gt 2 || $# -eq 1 && ! -d $1 ]] && echo "$file: ***WARNING: Skipping non-file"
    continue 	# Continue silently if processing a directory using pattern match
  fi

  # For common case of using current directory, strip off ./ that starts each file name
  if [ ${file:0:2} == "./" ] ; then
    file=${file:2}
  fi

  # Check if file content and extension correct
  ext=${file##*.}
  # <<< no longer working for here-documents in subsetted cygwin bootstrap used here; so change to temp file
  # HTS=$(htsfile "$file")
  htsfile "$file" > "HTS.tmp"
  if [[ "$ext" =~ ^(zip|7z|bz|bz2)$ ]]; then
    # htsfile cannot see in (to tell if FASTA) for .zip, .7z nor .bz2 archive formats (but does recognize bzip2 format)
    # So go ahead if those and check if a FASTA after recompressing; might recompress some non-bioinformatic files but ...
    # that is a risk the user takes if they do a general directory specification
    echo "$file: ***WARNING: Cannot look into .zip/.bz/.bz2/.7z files; so recompressing to peak inside."
  elif grep -v "FASTA" < "HTS.tmp" &> /dev/null ; then
    echo "$file: ***WARNING: Skipping a non-FASTA file"
    continue
  fi  # simply fall through as next check determines if need to recompress or not

  if grep "FASTA BGZF" < "HTS.tmp" &> /dev/null ; then
    # echo "$file: Already a BGZF compressed FASTA file!"   # simply keep silent
    filen="$file"
  else
    # OK, if here, is (likely) a FASTA but is not BGZF compressed
    echo "$file: Fixing Reference Genome File compression"
    case "$ext" in
      gz)		# BGZF already checked above; so must be gzipped
  { gzip -d -f "$file" ; bgzip -f -i -@ $np "${file%.*}" ; } && filen="$file"
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
  	    if grep -- "-compressed" < "HTS.tmp" &> /dev/null ; then
	echo "***WARNING: $file has a FASTA file extension only but is compressed. Missing a compressed file extension. Skipping. "
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
    if [ ! -f "$filen" ] ; then
      echo "$file: ***WARNING: Compression conversion failed (missing file). Skipping"
      continue
    fi
    #HTS=$(htsfile "$filen")
    htsfile "$filen" > "HTS.tmp"
    if grep -v "FASTA BGZF" < "HTS.tmp" &> /dev/null ; then
      echo "$file: ***WARNING: Compression conversion failed (not BGZF). Skipping"
      continue
    fi
  fi

  rm -f "HTS.tmp"
  # By this point, Ref Genome file ends in .gz, is FASTA and is BGZF compressed. Ready to create supporting indices, etc.
  #  Does allow explicitely named file through that is FASTA and BGZF compressed but without the extension. Some programs
  #  will not handle that but who are we to judge.

  # Now to (re)create any missing or old indices (using a possibly new filename $filen if compression changed)
  # GATK that uses the .dict file wants it to be named without {fasta,fa,fna}.gz extension
  # Reminder: BASH -nt is the "file newer than" operator
  filed=$(echo "$filen" | sed "s/.fasta.gz//;s/.fna.gz//;s/.fa.gz//")  # strip known extension
  [ "$filen" -nt "$filed.dict" ] && echo "$filen: Creating FA DICTionary"         && samtools dict "$filen" -o "$filed.dict"
  [ "$filen" -nt "$filen.fai" ]  && echo "$filen: Creating FA (FAI) Index"        && samtools faidx "$filen"    # Usually creates .gzi as well
  [ "$filen" -nt "$filen.gzi" ]  && echo "$filen: Creating BGZip (GZI) Index"     && bgzip -r "$filen"          # Also samtools index works

  # Those were the quick and easy ones; now the BWA Index which adds 5.5 GB of files and takes near an hour
  #   BWA Index creates .bwt (30 min, 3 GB), .pac (800MB), .ann, .amb, and .sa (10 min, 1.5GB); with 900M original, 6.4GB total!
  # [ "$filen" -nt "$filen.bwt" ]  && echo "$filen: Creating BWA Indices: 45 min, 5.5 GB" && bwa index "$filen"
  # Todo add parameter to optionally turn on BWA index capability

  # Finally, lets collect all the special WGS Extract desired stats to understand the FASTA content; all from DICT file
  [ "$filen" -nt "$filed.wgse" ] && echo "$filen: Creating WGS Extract Info"  && {
 
    # Throw out first line, take only columns 2 to 4, upcase and sort. Store temporarily
    tail -n +2 "$filed.dict" | cut -f2-4 | awk '{print toupper($0)}' | sort -d > "$filed.tmp"

    snct=$(wc -l "$filed.tmp" | cut -d' ' -f1)        				 # True SN count
    md5b=$(cut -f1-2 < "$filed.tmp" | md5sum | cut -d' ' -f1)  # MD5Sum BAM   (SN, LN)
    md5c=$(cut -f1-3 < "$filed.tmp" | md5sum | cut -d' ' -f1)  # MD5Sum CRAM  (SN, LN, MD)
    md5f=$(cut -f2-3 < "$filed.tmp" | md5sum | cut -d' ' -f1)  # MD5Sum FASTA (LN, MD)

	  if [ $snct -lt 83 ]; then
	    errp=" ***WARN: Too few SN entries (<84)***"
	  else
  	  errp=""		# Null out error report string that will become tail of entry
	  fi

    # Model 36, 37 and 38 chromosome length counts (one long pattern line)
    chrs="LN:247249719\|LN:249250621\|LN:248956422\|LN:242951149\|LN:243199373\|LN:242193529\|LN:199501827\|LN:198022430\|LN:198295559\|LN:191273063\|LN:191154276\|LN:190214555\|LN:180857866\|LN:180915260\|LN:181538259\|LN:170899992\|LN:171115067\|LN:170805979\|LN:158821424\|LN:159138663\|LN:159345973\|LN:146274826\|LN:146364022\|LN:145138636\|LN:140273252\|LN:141213431\|LN:138394717\|LN:135374737\|LN:135534747\|LN:133797422\|LN:134452384\|LN:135006516\|LN:135086622\|LN:132349534\|LN:133851895\|LN:133275309\|LN:114142980\|LN:115169878\|LN:114364328\|LN:106368585\|LN:107349540\|LN:107043718\|LN:100338915\|LN:102531392\|LN:101991189\|LN:88827254\|LN:90354753\|LN:90338345\|LN:78774742\|LN:81195210\|LN:83257441\|LN:76117153\|LN:78077248\|LN:80373285\|LN:63811651\|LN:59128983\|LN:58617616\|LN:62435964\|LN:63025520\|LN:64444167\|LN:46944323\|LN:48129895\|LN:46709983\|LN:49691432\|LN:51304566\|LN:50818468\|LN:154913754\|LN:155270560\|LN:156040895\|LN:57772954\|LN:59373566\|LN:57227415"
    # Some alt contigs same length as chromosome; so further refine by primary entry SNs
    chrsn="SN:CHR[123456789XY]\>\|SN:CHR1[0123456789]\>\|SN:CHR2[012]\>\|SN:[123456789XY]\>\|SN:1[0123456789]\>\|SN:2[012]\>\|CM0006\|SN:NC_0000[012][0123456789]"

    # Capture Primary (Chromosomes) key stats
    # <<< here-document no longer works in subsetted cygwin64 bootstrap release here; so change to temp file
    #primary=$(grep -w -e "$chrs" "$filed.tmp" | grep -e "$chrsn") 	# Length is not enough; need to then filter by SN
    grep -w -e "$chrs" "$filed.tmp" | grep -e "$chrsn" > "$filed.ptmp"  # Length is not enough; need to then filter by SN
	  md5p=$(cut -f2-3 "$filed.ptmp" | md5sum | cut -d' ' -f1)  # MD5Sum Chrs-only (LN, MD)
    pcnt=$(wc -l "$filed.ptmp" | cut -d' ' -f1)
    if [ $pcnt -ne "24" ]; then
      # Todo simplify to grep -c and drop wc and cut
      ycnt=$(grep -w -e "LN:59373566\|LN:57227415" "$filed.ptmp" | wc -l | cut -d' ' -f1)
      if [ $ycnt -ne "1" ]; then
        # M and Y are only duplicates found so far; so give a special error if found to be Y
        echo "$filen: ***ERROR: 1 expected, $ycnt Y chromosome entries found in ref model"
        errp="$errp *** ERROR:Y $pcnt!=24 ***"
      else
        echo "$filen: ***ERROR: 24 expected, $pcnt chromosomes found in ref model" 
        errp="$errp *** ERROR:P $pcnt!=24 ***"
      fi
    fi

    rm -f "$filed.ptmp"

	  # Capture MT line to get SN name, LN, and M5 hash -- not directly captured in primary above
	  # <<< here-document no longer works in subsetted cygwin64 bootstrap release here; so change to temp file
    chrM=$(grep -w -e "LN:16569\|LN:16571" < "$filed.tmp")	# Keep SN and MD columns to help capture detail
    echo "$chrM" > "$filed.mtmp"
    mcnt=$(wc -l < "$filed.mtmp" | cut -d' ' -f1)
    if [ $mcnt -ne "1" ]; then
      echo "$filen: ***ERROR: 1 expected, $mcnt microchonrdrial entries found in ref model"
      errp="$errp *** ERROR:M $mcnt!=1 ***"
	    chrM=$(head -n 1 "$filed.mtmp")		# Truncate to first line for reporting
    fi

	  # Setup trailing ID based on Sequence Naming convention
    chrMSN=$(cut -f1 < "$filed.mtmp")
    case $chrMSN in
	    SN:CHRM)		# UCSC chrN with single numeric / alphabetic naming
	      msn="h"
	      ;;
	    SN:CHRMT)		# UCSC chrN naming with oddball chrMT
	      msn="ht"
	      errp="$errp ***WARN: chrMT name is non-standard ***"
	      ;;
	    SN:MT)		  # EBI single mumeric / alphabetic naming except MT
	      msn="g"
	      ;;
	    SN:NC_012920.1)	# NCBI RefSeq Acquisition ID naming
	      msn="n"
	      ;;
	    SN:J01415.2)		# NCBI GenBank ID naming
	      msn="c"
	      ;;
	    *)			# Unrecognized naming convention; or no MT in Fasta model file?
	      msn="x"
	      echo "$filen: ***ERROR: Unregonized Sequence Naming"
	      errp="$errp ***ERROR: Unrecognized SN Name***"
	      ;;
	  esac
	
	  # Determine model Build variation based on mitochondrial model found
	  #  Use special of 19 for Yoruba Build 19/37 and 37 for rCRS
    chrMM5=$(cut -f3 < "$filed.mtmp")
	  case $chrMM5 in
	    M5:C68F52674C9FB33AEF52DCF399755519)		# rCRS
	      mbuild="37"
	      ;;
	    M5:D2ED829B8A1628D16CBEEE88E88E39EB)		# Yoruba
	      mbuild="19"
	      ;;
	    *)										# Unrecognized model; RSRS?
	      mbuild="xx"
	      echo "$filen: ***ERROR: Unregonized Mitochondrial Model"
	      errp="$errp ***ERROR: Unrecognized Mito Model***"
	      ;;
	  esac

	  rm -f "$filed.mtmp"
	
	  # Based on Chromosome-only MD5sum of LN and M5 fields; determine major/minor Build type
	  # Mid-release, the MD5Sum of the Primary chromosom-only LN and M5 fields changed for every entry (consistently).
	  #  Oddly. the MD5Sum of the total file SN, LN and M5 fields did not change though. So we simply include the
	  #  original followed by the new until we figure this out.
	  case $md5p in
	    04ecc3be916bf559dd924885b94e52e3|b05113b52031beadfb6737bc1185960b)		# HG19/37
	      build="HG${mbuild}" ;;
	    1710809e0f974613edbbf77bf161ff47|5a23f5a85bd78221010561466907bf7d)		# EBI37
	      build="EBI37"       ;;
	    1ed5687e6d6e4431e900e9c8487b7805)		# NCB38
	      build="NCB38"	      ;;
	    22698c26453cea2733359b32e379bcf1|eec5eb2eeae44c48a31eb32647cd04f6)		# EBI38
	      build="EBI38" 	    ;;
	    31ee75abc53d3b9e8b612cf1c8c43d30|7a5eb72fb45c4567431651aa6f9edfef)		# 1K19/37
	      build="1K${mbuild}" ;;
	    3edede0915ebe8c49b002cae10067e54)		# EBI37 with extra Y (error)
	      build="EBI37p"	    ;;
	    47cf1169bb18701afada92e43f5daba7)		# NCB37
	      build="NCB37" 	    ;;
	    76419821763d9ae311d30a2790f0ae5e|4bdbf8a3761d0cd03b53a398b6da026d)		# HG38
	      build="HG38"  	    ;;
	    8716c91fb4131e44549ae12196fb6b8a)		# EBI38p
	      build="EBI38p"	    ;;
	    a001f69d019f8a2c4f9591bb414fc56e)		# EBI18 (Yoruba)
	      build="EBI18"
	      mbuild="18"   	    ;;
	    a8a70bbb2a27102e7fd74506a8d1da6b|b7884451f3069579e5f2e885582b9434)		# 1K38
	      build="1K38"	      ;;
	    e629a1ffd866c84db801688ecb178174)		# 1k38p (odd Verily one)
	      build="1K38p" 	    ;;
	    f65dd5cce2a25e509dd3c5d78e475c25)		# 1K19/37 with extra Y (error)
	      build="1K${mbuild}p"   ;;
	    *)
	      build="UNK"
	      errp="$errp ***ERROR: Unrecognized Build Model***"
	      ;;
	  esac
	  build="${build}${msn}"			# Append sequence naming convention used
	
    # Save all the values as single line, tab separated text file; append to project file for directory
	  #  Note: chrM is already a 3 value, tab separated variable that we simply use to create 3 columns
    echo "\"$filen\"	$build	$snct	$md5b	$md5c	$md5f	$md5p	$chrM	\"$errp\"" > "$filed.wgse"
    cat "$filed.wgse" >> "WGSE.csv"
	  rm "$filed.tmp"

    # Final error check -- in case previous uncaught error left intermediate .wgse file or multi-line stats (mtDNA)
    wcnt=$(wc -l "$filed.wgse" | cut -d' ' -f1)
    if [ $wcnt -ne "1" ]; then
      echo "$filen: ***ERROR: Failed generating final WGSE stats file"
      rm -f "$filed.wgse"
    fi
  }
done
