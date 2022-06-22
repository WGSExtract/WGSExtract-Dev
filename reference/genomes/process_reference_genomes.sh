#!/usr/bin/env bash
#
# Script to process a folder of Reference Genome files in preparation for use in Bioinformatic tools
# 
# If a directory specified, the script naively presumes any file in the folder with the correct extension
#   is a FASTA file.
#
# Part of the Reference Genome package in WGS Extract (https://wgsextract.github.io/)
# Copyright (c) 2020-22 Randy Harr
#

#
# This script works off the current directory or a file parameter list.
#   * If parameter list is one or more files, it works off that.
#   * If a single directory specified, it works off the matching files in that directory.
# If a directory specified, it tries to find all Reference Genomes (actually, FASTA files) using common extensions.
# It uses each identified FASTA as the base.  
# It can be run a second time and will simply do any needed updates. So if errors fixed, simply rerun. Only difference
#  is the WGSE.csv file is appended to on each run if already existing.
#
# NOTE: We have commented out the BWA Index command due to the extensive resource usage (CPU and file space) for
#  doing BWA indices.  It is simply part of the main python code and run as needed.  Also different aligners have
#  different index files as it is.
#

# Number of Processors to use when available (Todo should read from system)
np=16

# Todo Programs used; maybe parameterize for platform so $PATH not required?
# gunzip, unzip, bunzip2, 7z, sort, cut, (g)awk, sed, grep, tail, head, wc, md5sum, rm, bgzip, htsfile, samtools

# Figure out what we have to do in the current directory; or if explicitly set in the parameter list
shopt -s nullglob
sumFPB="WGSE"
if [[ $# -eq 1 && "$1" = "clean" ]] ; then
  rm -f *fai *gzi *dict *wgse "${sumFPB}.csv" "${sumFPB}_dict*csv" "${sumFPB}_uniq*csv" || true   # Files created here
  exit
elif [[ $# -eq 1 && "$1" = "clean_all" ]] ; then
  rm -f *fai *gzi *dict *wgse "${sumFPB}.csv" "${sumFPB}_dict*csv" "${sumFPB}_uniq*csv" || true   # Files created here
  rm -f *amb *ann *bwt *pac *sa || true   # BWA Index files
  rm -f *fa.gz *fasta.gz *fna.gz || true  # Downloaded reference genomes; to get back to true delivery folder
  exit
fi
if [[ $# -eq 1 && -d "$1" ]] ; then
  file_list=("$1"/*.{fa,fna,fasta,7z,gz,zip,bz,bz2})
  if [[ $1 != "." ]] ; then
    sumFPB="$1"/WGSE
  fi
  # Want to process whole directory; so start with clean slate of summary files
  rm -f "${sumFPB}.csv" "${sumFPB}_dict.csv" "$1"/*wgse || true
  rm -f "${sumFPB}_uniq_build.csv" "${sumFPB}_uniq_ChrLNM5.csv" "${sumFPB}_uniq_SNcnt.csv" || true
  rm -f "${sumFPB}_dict_uniq_SNLNM5.csv" "${sumFPB}_dict_uniq_LNM5.csv" || true
  wholedir="True"
else
  file_list=("$@")
  wholedir="False"
fi
if [[ $# -eq 0 || ${#file_list[@]} -eq 0 ]] ; then
  printf "Usage: %s [dir]\n     # All files in folder specified; use dot (.) for current\n" "$0"
  printf "  or   %s [file(s)]   # List of specified files\n" "$0"
  printf "  or   %s clean       # Files created by this script; in current script folder\n" "$0"
  printf "  or   %s clean_all   # All files including downloaded reference genomes; in current script folder\n" "$0"
  printf "Folder searches use any file with extension: fa, fna, fasta, gz, zip, bz, bz2, or 7z as a FASTA file\n"
  exit
fi

# MacOS is overriding PATH and using builtin (ancient) BASH; so need absolute path override.
# Make sure we have the bioinformatic tools on the path
case $OSTYPE in
  linux*)           bashx="/usr/bin/env bash"  ;;
  msys*|cygwin*)    bashx="/bin/bash.exe"  ;;  # PATH="/usr/local/bin:/bin:${PATH}"
  darwin*)          bashx="/opt/local/bin/bash" ;  PATH="/opt/local/bin:${PATH}"  ;;
  *)  printf "*** Error: unknown OSTYPE of %s\n" "$OSTYPE" &&  exit 1  ;;
esac

printf "Processing Reference Genome file(s) ...\n"
# echo "***DEBUG:${#file_list[@]} parameters: ${file_list[@]}"

# Setup initial csv summary file header if not yet written; for when initial processing of whole directory (TSV)
# Columns are: File Name, Major/Minor build code, SN Cnt ; MD5Sums: SN/LN, SN/LN/M5, LN/M5, Chromo SN, Chromo LN/M5 ;
#              Mito: SN, LN, M5 ; Error Message(s)
# Chromo LN/M5 identifies Major Build, Full LN/M5 the Minor Class (and also SN Cnt)
if [[ ! -f "${sumFPB}.csv" ]] ; then
  printf $'File\tBuild\tSN_CNT\tBAM_(SN,_LN)\tCRAM_(SN,_LN,_M5)\tFASTA_(LN,_M5)\tChromo_(LN,_M5)\tChromo_(SN)\tMito:_SN\tLN\tM5\tError\n' > "${sumFPB}.csv"
fi
	
LANG=POSIX		# Needed for sort command
for file in "${file_list[@]}" ; do
  # echo "***DEBUG: file = $file"

  # For common case of using current directory, strip off ./ that starts each file name
  if [[ "${file:0:2}" == "./" ]] ; then
    file=${file:2}
  fi

  fpa=${file%/*}    # File path (if it exists)  (not used; but keep in case we need it later)
  fbn=${file##*/}   # File base name
  ext=${file##*.}   # File (last) extension

  # Check if even a file
  if [[ ! -f "$file" ]] ; then
    [[ $# -gt 2 || $# -eq 1 && ! -d $1 ]] && echo "$fbn: ***WARNING: Skipping non-file"
    continue 	# Continue silently if processing a directory using pattern match
  fi

  # Check if file content and extension correct
  # <<< here-documents from variable not working in subsetted cygwin bootstrap; so change to temp file
  # HTS=$(htsfile "$file")
  htsfile "$file" > "HTS.tmp"

  if [[ "$ext" =~ ^(zip|7z|bz|bz2)$ ]]; then
    # htsfile cannot see in (to tell if FASTA) for .zip, .7z nor .bz2 archive formats (but does recognize bzip2 format)
    # So go ahead if those and check if a FASTA after recompressing; might recompress some non-bioinformatic files but ...
    # that is a risk the user takes if they do a general directory specification
    echo "$fbn: ***WARNING: Cannot look into .zip/.bz/.bz2/.7z files; so recompressing to peak inside."
  elif grep -v "FASTA" < "HTS.tmp" &> /dev/null ; then
    echo "$fbn: ***WARNING: Skipping a non-FASTA file"
    rm -f "HTS.tmp" || true
    continue
  fi  # simply fall through as next check determines if you need to recompress or not

  if grep "FASTA BGZF" < "HTS.tmp" &> /dev/null ; then
    # echo "$file: Already a BGZF compressed FASTA file!"   # simply keep silent
    filen="$file"
  else
    # OK, if here, is (likely) a FASTA but is not BGZF compressed
    echo "$fbn: Fixing Reference Genome File compression"
    case "$ext" in
      gz)		# BGZF already checked above; so must be gzipped
        { gzip -d -f "$file" ; bgzip -if@ $np "${file%.*}" ; } && filen="$file"
  	    # This is the only format where you end up with the same file extension; they should have used .bgz but ...
  	    ;;

      zip)		# Likely zip'ped 
	    { unzip -p "$file" | bgzip -cf@ $np > "${file%.*}.gz" ; } && rm -f "$file" && filen="${file%.*}.gz"
  	    ;;

      7z)		# A popular, very high compression format not normally seen in Bioinformatics
	    { 7z e -mmt$np -so "$file" | bgzip -cf@ $np > "${file%.*}.gz" ; } && rm -f "$file" && filen="${file%.*}.gz"
  	    ;;

      bz | bz2)		# A popular Linux/Unix/MacOSX high compression format sometimes seen in Bioinformatics
        { bunzip2 -cf "$file" | bgzip -cf@ $np > "${file%.*}.gz" ;} && rm -f "$file" && filen="${file%.*}.gz"
  	    ;;

      fa | fna | fasta)	# Appears likely not compressed; but will check
  	    if grep -- "-compressed" < "HTS.tmp" &> /dev/null ; then
	echo "$fbn: ***WARNING: File is compressed but with only a FASTA extension. Need a .gz/.zip extension. Skipping. "
	      continue
	    fi
  	    bgzip -if@ $np "$file"
	    filen="$file.gz"
  	    ;;

      *)		# Well that is embarrising
	    echo "$fbn: ***WARNING: Unknown file extension (internal error)"
	    ;;
    esac
    # Above does not handle .zip / .7z archive where multiple files (i.e. like tar+gzip combo from original pkzip)
    #  Singular, internal file name must be the same name as archive file without extension

    # Check again if matched file with extension is a valid, compressed FASTA file just to make sure
    if [[ ! -f "$filen" ]] ; then
      echo "$fbn: ***WARNING: Compression conversion failed (missing file). Skipping"
      continue
    fi
    #HTS=$(htsfile "$filen")
    htsfile "$filen" > "HTS.tmp"
    if grep -v "FASTA BGZF" < "HTS.tmp" &> /dev/null ; then
      echo "$fbn: ***WARNING: Compression conversion failed (not BGZF). Skipping"
      continue
    fi
  fi

  rm -f "HTS.tmp" || true
  # By this point, Ref Genome file ends in .gz, is FASTA and is BGZF compressed. Ready to create supporting indices, etc.
  #  Does allow explicitely named file that is FASTA and BGZF compressed but without the proper extension. Some programs
  #  will not handle that but who are we to judge.

  # Now to (re)create any missing or old indices (using a possibly new filename $filen if compression changed)
  # GATK that uses the .dict file wants it to be named without {fasta,fa,fna}.gz extension so we oblige
  # Reminder: BASH -nt is the "file newer than" operator
  filed=$(echo "$filen" | sed "s/.fasta.gz//;s/.fna.gz//;s/.fa.gz//")  # strip known extension
  fbnn=${filen##*/}         # We use the whole filename for status reporting; so if changed need a new basename
  [ "$filen" -nt "$filed.dict" ] && echo "$fbnn: Creating FA DICTionary"      && samtools dict "$filen" -o "$filed.dict"
  [ "$filen" -nt "$filen.fai" ]  && echo "$fbnn: Creating FA (FAI) Index"     && samtools faidx "$filen"    # Usually creates .gzi as well
  [ "$filen" -nt "$filen.gzi" ]  && echo "$fbnn: Creating BGZip (GZI) Index"  && bgzip -r "$filen"          # Also samtools index works

  # Those were the quick and easy ones; now the BWA Index which adds 5.5 GB of files and takes near an hour
  #   BWA Index creates .bwt (30 min, 3 GB), .pac (800MB), .ann, .amb, and .sa (10 min, 1.5GB); with 900M original, 6.4GB total!
  # [ "$filen" -nt "$filen.bwt" ]  && echo "$filen: Creating BWA Indices: 45 min, 5.5 GB" && bwa index "$filen"
  # Todo add parameter to optionally turn on BWA index capability

  # Finally, lets collect all the special WGS Extract desired stats to understand the FASTA content; all from DICT file
  [ "$filen" -nt "$filed.wgse" ] && echo "$fbnn: Creating WGS Extract Info"  && {
 
    # Throw out first line, take only columns 2 to 4, upcase and sort. Store temporarily. This is the key file.
    tail -n +2 "$filed.dict" | cut -f2-4 | awk '{print toupper($0)}' | sort -d > "$filed.tmp"
    # =sort -V would be optimal when manually viewing the result; but -V is not portable (GNU Util only)

    md5b=$(cut -f1-2 < "$filed.tmp" | md5sum | cut -d' ' -f1)     # MD5Sum BAM   (SN, LN)
    md5c=$(cut -f1-3 < "$filed.tmp" | md5sum | cut -d' ' -f1)     # MD5Sum CRAM  (SN, LN, MD)
    md5f=$(cut -f2-3 < "$filed.tmp" | md5sum | cut -d' ' -f1)     # MD5Sum FASTA (    LN, MD)
    # DO NOT sort in the above; original sort is based on the SN field. So md5f changes if you sort again

    # SN Count (whole DICT), Primary SN count (just chromosomes and mito), Y chromosome count, Mitochondria count
    declare -i scnt pcnt ycnt mcnt

    # Smallest Analysis model is 84 SN entries; if less than likely an error (25 is exception for hg_wgse and t2t)
    snct=$(wc -l < "$filed.tmp")                    # True SN count
	if (( $snct < 83 && $snct != 25 )); then
	  errp=" ***WARN: Too few SN entries (<84, !=25)***"      # First setting; no need to concat
	else
  	  errp=""		# Null out error report string that will become tail of WGSE entry
	fi

    # One long pattern line for GREP of Chromosome LN fields in major model builds.
    # Used to extract just the needed primary chromosomes. In lieu of using SNs which may be non-standard
    # Note: some alt contigs have the same length as the chromosomes so must use second filter based on name anyway
    # chrN will not match because the DICT file content was up-cased. Purely there for our documentation
    # Build         36,          37,          38   T2T chm13: v1.1     v1.0    (also T2T HG002 v2 and v2.7)
    chrs="LN:247249719\|LN:249250621\|LN:248956422\|LN:248387328\|LN:248387497"         # chr1
    chrs+="\|LN:242951149\|LN:243199373\|LN:242193529\|LN:242696752\|LN:242696747"      # chr2
    chrs+="\|LN:199501827\|LN:198022430\|LN:198295559\|LN:201105948\|LN:201106605"      # chr3
    chrs+="\|LN:191273063\|LN:191154276\|LN:190214555\|LN:193574945\|LN:193575430"      # chr4
    chrs+="\|LN:180857866\|LN:180915260\|LN:181538259\|LN:182045439\|LN:182045437"      # chr5
    chrs+="\|LN:170899992\|LN:171115067\|LN:170805979\|LN:172126628\|LN:172126870"      # chr6
    chrs+="\|LN:158821424\|LN:159138663\|LN:159345973\|LN:160567428\|LN:160567423"      # chr7
    chrs+="\|LN:146274826\|LN:146364022\|LN:145138636\|LN:146259331\|LN:146259322"      # chr8
    chrs+="\|LN:140273252\|LN:141213431\|LN:138394717\|LN:150617247\|LN:150617274"      # chr9
    chrs+="\|LN:135374737\|LN:135534747\|LN:133797422\|LN:134758134\|LN:134758122"     # chr10
    chrs+="\|LN:134452384\|LN:135006516\|LN:135086622\|LN:135127769\|LN:135127772"     # chr11
    chrs+="\|LN:132349534\|LN:133851895\|LN:133275309\|LN:133324548\|LN:133324781"     # chr12
    chrs+="\|LN:114142980\|LN:115169878\|LN:114364328\|LN:113566686\|LN:114240146"     # chr13
    chrs+="\|LN:106368585\|LN:107349540\|LN:107043718\|LN:101161492\|LN:101219177"     # chr14
    chrs+="\|LN:100338915\|LN:102531392\|LN:101991189\|LN:99753195\|LN:100338308"      # chr15
    chrs+="\|LN:88827254\|LN:90354753\|LN:90338345\|LN:96330374\|LN:96330493"          # chr16
    chrs+="\|LN:78774742\|LN:81195210\|LN:83257441\|LN:84276897\|LN:84277185"          # chr17
    chrs+="\|LN:76117153\|LN:78077248\|LN:80373285\|LN:80542538\|LN:80542536"          # chr18
    chrs+="\|LN:63811651\|LN:59128983\|LN:58617616\|LN:61707364\|LN:61707359"          # chr19
    chrs+="\|LN:62435964\|LN:63025520\|LN:64444167\|LN:66210255\|LN:66210247"          # chr20
    chrs+="\|LN:46944323\|LN:48129895\|LN:46709983\|LN:45090682\|LN:45827691"          # chr21
    chrs+="\|LN:49691432\|LN:51304566\|LN:50818468\|LN:51324926\|LN:51353906"          # chr22
    chrs+="\|LN:154913754\|LN:155270560\|LN:156040895\|LN:154259566\|LN:154259625\|LN:154343774\|LN:154349815" # chrX
    chrs+="\|LN:57772954\|LN:59373566\|LN:57227415\|LN:62456832\|LN:62460029" # chrY (no chm13 entry)
    # Some alt contigs have same length as chromosome; so further refine by known primary entry SNs
    chrsn="SN:CHR[1-9XY]\>\|SN:CHR1[0-9]\>\|SN:CHR2[0-2]\>"      # HGP naming
    chrsn+="\|SN:[1-9XY]\>\|SN:1[0-9]\>\|SN:2[012]\>"            # EBI naming
    chrsn+="\|CM0006\|SN:NC_0000[0-2][0-9]\|CP0682[567]\|CP08656[89]\|SN:CHR[XY]_HG002"    # Accession, begins with
    # Alt contig names in HGP and EBI start with full chromosome name; so make sure checking whole field
    # In accession naming, alt contigs have their own accession entries that start differently

    # Capture Primary (Chromosomes) key stats
    # <<< here-document is  not working in subsetted cygwin64 of Win10; so change to temp file use
    #primary=$(grep -w -e "$chrs" "$filed.tmp" | grep -e "$chrsn") 	# Length is not enough; need to then filter by SN
    grep -w -e "$chrs" "$filed.tmp" | grep -e "$chrsn" > "$filed.ptmp"  # Need to filter by LN then SN also
    # Do not sort again; keep same order as original sort of whole DICT SQ entries above.

    md5s=$(cut -f1   "$filed.ptmp" | md5sum | cut -d' ' -f1)  # MD5Sum Chrs-only Names (SN) (not really needed)
    md5p=$(cut -f2-3 "$filed.ptmp" | md5sum | cut -d' ' -f1)  # MD5Sum Chrs-only FASTA (   LN, MD) (IMPORTANT)

	# LNs are for Build 36-38 and T2T HG002 Y model v2 and v2.7; if cannot find 24 primary then an issue
    pcnt=$(wc -l < "$filed.ptmp")
    ycnt=$(grep -c -w -e "LN:57772954\|LN:59373566\|LN:57227415\|LN:62456832\|LN:62460029" "$filed.ptmp")
    # First check if missing Y or more than one Y entry; report that as issue if so before checking pcnt
    if (( $ycnt != 1 )) ; then
      # M and Y are only duplicates found so far; so give a special error if found to be Y here; M is later
      echo "$fbnn: ***ERROR: 1 expected, $ycnt Y chromosome entries found in ref model"
      errp="$errp *** ERROR:Y $ycnt!=1 ***"
    elif (( $pcnt != 24 )) ; then
      echo "$fbnn: ***ERROR: 24 expected, $pcnt chromosomes found in primary ref model"
      errp="$errp *** ERROR:P $pcnt!=24 ***"
    fi
    rm -f "$filed.ptmp" || true

	# Capture the MT line to get SN name, LN, and M5 hash -- not directly captured in primary above
	# <<< here-document no longer works in subsetted cygwin64 bootstrap release here; so change to temp file
	# Variable may be multi-line (contain \n). So dump to file and use head to just take first one.
    chrM=$(grep -w -e "LN:16569\|LN:16571" < "$filed.tmp")	# Keep SN and MD columns to help capture detail
    echo "$chrM" > "$filed.mtmp"
    mcnt=$(wc -l < "$filed.mtmp")
    if (( $mcnt != 1 )) ; then
      echo "$fbnn: ***ERROR: 1 expected, $mcnt mitrochondrial entries found in ref model"
      errp="$errp *** ERROR:M $mcnt!=1 ***"
	  chrM=$(head -n 1 "$filed.mtmp")       # Truncate to first line for further reporting
	  echo "$chrM" > "$filed.mtmp"
    fi

	# Setup trailing ID based on Sequence Naming convention (could use md5s from primary chromosomes also)
    chrMSN=$(cut -f1 "$filed.mtmp")
    case $chrMSN in
	  SN:CHRM)		    # UCSC chrN with single numeric / alphabetic naming
	    msn="h"
	    ;;
	  SN:CHRMT)		    # UCSC chrN naming with oddball chrMT
	    msn="ht"
	    errp="$errp ***WARN: chrMT name is non-standard ***"
	    ;;
	  SN:MT)          # EBI single mumeric / alphabetic naming except MT
	    msn="g"
	    ;;
	  SN:NC_012920.1)	# NCBI RefSeq Acquisition ID naming
	    msn="n"
	    ;;
	  SN:J01415.2|SN:CP068254.1|"SN:GI|113200490|GB|J01415.2|HUMMTCG")	# NCBI GenBank ID naming (both GRCh and T2T)
	    msn="c"
	    ;;
	  *)			    # Unrecognized naming convention; or no MT in Fasta model file?
	    msn="x"
	    echo "$fbnn: ***ERROR: Unregonized Sequence Naming (${chrMSN})"
	    errp="$errp ***ERROR: Unrecognized SN Name***"
	    ;;
	esac

	chrMLN=$(cut -f2 "$filed.mtmp")

	# Determine model Build variation based on mitochondrial model found
	#  Use special of 19 for Yoruba Build 19/37 and 37 for rCRS whether Build 37 or 38
    chrMM5=$(cut -f3 "$filed.mtmp")
	case $chrMM5 in
	  M5:C68F52674C9FB33AEF52DCF399755519)	# rCRS
	    mbuild="37"
	    ;;
	  M5:D2ED829B8A1628D16CBEEE88E88E39EB)	# Yoruba
	    mbuild="19"
	    ;;
	  M5:EC493A132AC4823AA696E37109F64972|M5:2AEA08C58600A30435E4302A82481DC0)    # chm13 T2T model
	    mbuild="99"
	    ;;
	  *)		# Unrecognized model; RSRS?
	    mbuild="xx"
	    echo "$fbnn: ***ERROR: Unregonized Mitochondrial Model (${chrMM5})"
	    errp="$errp ***ERROR: Unrecognized Mito Model***"
	    ;;
	esac

	rm -f "$filed.mtmp" || true
	
	# Based on Primary-Chromosome-only MD5sum of LN and M5 fields; determine major/minor Build type (ignoring SN names)
	# Bug introduced by removing sort caused fields to change mid-release,so some entries have old and new.
	case $md5p in
	  04ecc3be916bf559dd924885b94e52e3|b05113b52031beadfb6737bc1185960b)	# HG19/37
	    build="HG${mbuild}"   ;;
	  1710809e0f974613edbbf77bf161ff47|5a23f5a85bd78221010561466907bf7d)	# EBI37
	    build="EBI37"         ;;
	  1ed5687e6d6e4431e900e9c8487b7805|bee8aebc6243ff5963c30abbd738d1f6)	# NCBI 38 Genbank (all except p14)
	    build="NCB38"         ;;
	  22698c26453cea2733359b32e379bcf1|eec5eb2eeae44c48a31eb32647cd04f6)	# EBI38
	    build="EBI38" 	    ;;
	  31ee75abc53d3b9e8b612cf1c8c43d30|7a5eb72fb45c4567431651aa6f9edfef)	# 1K19/37
	    build="1K${mbuild}"   ;;
	  3edede0915ebe8c49b002cae10067e54)     # EBI37 with extra Y (error)
	    build="EBI37p"	    ;;
	  4136c29467b6757938849609bedd3996)     # NCBI38p14 GenBank
	    build="NCB38"         ;;
	  46cf0768c13ec7862c065e45f58155bf)     # EBI18
        build="EBI18"         ;;
	  47cf1169bb18701afada92e43f5daba7)     # NCB37
	    build="NCB37" 	    ;;
	  4bf6c704e4f8dd0d31a9bf305df63ed3)     # T2T CHM13 v1.1 with HG002 xy v2.7
	    build="THGv27"        ;;
	  4d0aa9b8472b69f175d279a9ba8778a1)     # HPP CHM13 v1.1 with GRCh38 Y
        build="HPPv11"        ;;
      5e16e3cbdcc7b69d21420c332deecd3b)     # CHM13 v1.0 from T2T Consortia (original)
        build="T2Tv10"        ;;
      591bb02c89ed438566ca68b077fee367)     # Errored EBI GRCh37 models (extra Y)
        build="1K37"          ;;
	  65a05319ad475cf51c929d3b55341bc2)     # T2T CHM13 v1.1 with HG002 xy v2
	    build="THGv20"        ;;
	  6b3d36e85f3df87d3a0c45327fe5bc41)     # CHM13 v1.0 from T2T Consortia (GENBANK)
	    build="T2Tv10"        ;;
	  76419821763d9ae311d30a2790f0ae5e|4bdbf8a3761d0cd03b53a398b6da026d)	# HG38
	    build="HG38"  	    ;;
	  7cee777f1939f4028926017158ed5512|13cbd449292df5bd282ff5a21d7d0b8f)  # Final T2T v2.0 (CHM13 v1.1 w/ HG002 v2.7 Y)
	    build="T2Tv20"      ;;              # Original Accession naming as well as standard HG naming we made here
	  84e78573982f3ea293bfeb54cd529309)     # Verily oddball GRCh38
	    build="1K38p"         ;;
	  8716c91fb4131e44549ae12196fb6b8a)     # EBI38p
	    build="EBI38p"	    ;;
	  a001f69d019f8a2c4f9591bb414fc56e)     # EBI18 (Yoruba)
	    build="EBI18"
	    mbuild="18"   	    ;;
	  a2fe6ab831d884104783f9be437ddbc0)     # Errored EBI GRCh38 models
	    build="EBI38p"        ;;
	  b7884451f3069579e5f2e885582b9434)     # 1K38
	    build="1K38"          ;;
	  bbd2cf1448ccc0eaa2472408fa9d514a)     # ySeq HG38 w/ HG002 v2 Y
        build="THGv20p"       ;;
	  ca2e97bc5ecff43a27420eee237dbcc3)     # Errored EBI GRCh37 models (extra Y)
	    build="EBI37"         ;;
	  c182b40ef3513ef9a1196881a4315392)     # HPP CHM13 v1 with GRCh38 Y
        build="HPPv1"         ;;
      d41d8cd98f00b204e9800998ecf8427e)     # CHM13 v0.9 from T2T Consortia (original) (ERROR: no LN's in array above)
        build="T2Tv09"        ;;
	  e9438f38ad1b9566c15c3c64a9419d9d)     # CHM13 v1.1 from T2T Consortia (original)
        build="T2Tv11"        ;;
	  e629a1ffd866c84db801688ecb178174)     # 1k38p (odd Verily one)
	    build="1K38p" 	    ;;
	  f65dd5cce2a25e509dd3c5d78e475c25)     # 1K19/37 with extra Y (error)
	    build="1K${mbuild}p"  ;;
      f7c76dbcf8cf8b41d2c1d05c1ed58a75)     # NCBI RefSeq GRCh37
        build="NCB37"        ;;
      fdf62efdd97b08a1b4a0e9bc183820e6)     # CHM13 v1.1 from T2T Consortia (GenBANK)
        build="T2Tv11"        ;;
	  *)
	    build="UNK"
	    echo "$fbnn: ***ERROR: Unknown Build (${md5p})"
	    errp="$errp ***ERROR: Unrecognized Build Model***"
	    ;;
	esac
	build="${build}${msn}"			# Append sequence naming convention used
	
    # Save all the values as single line, tab separated text file; append to project file for directory
	  #  Note: chrM is already a 3 value, tab separated variable that we simply use to create 3 columns
    printf $'\"%s\"\t%s\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\"%s\"\n' \
      "$filen" $build $snct $md5b $md5c $md5f $md5p $md5s $chrMSN $chrMLN $chrMM5 "$errp" > "$filed.wgse"

    # Final error check -- in case previous uncaught error left intermediate .wgse file or multi-line stats (mtDNA)
    wcnt=$(wc -l < "$filed.wgse")
    if (( $wcnt != 1 )) ; then
      echo "$fbnn: ***ERROR: Failed generating final WGSE stats file"
      rm -f "$filed.wgse" || true
    else
      cat "$filed.wgse" >> "${sumFPB}.csv"
    fi

    rm "$filed.tmp" || true
  }
done

# Some basic stats for the Reference Genome study document when run on a directory with all known models
if [ "$wholedir" = "True" ]; then
  echo "Entries: " $(wc -l "${sumFPB}.csv")
  cut -f2    "${sumFPB}.csv" | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $1,$2}'> "${sumFPB}_uniq_build.csv"
  echo "Entries: " $(wc -l "${sumFPB}_uniq_build.csv")
  cut -f8    "${sumFPB}.csv" | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $1,$2}'> "${sumFPB}_uniq_ChrLNM5.csv"
  echo "Entries: " $(wc -l "${sumFPB}_uniq_ChrLNM5.csv")
  cut -f3    "${sumFPB}.csv" | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $1,$2}'> "${sumFPB}_uniq_SNcnt.csv"
  echo "Entries: " $(wc -l "${sumFPB}_uniq_SNcnt.csv")
  # For convenience, create a table of all SNs in all FASTA's in the specified directory
  grep "SN:" "$1"/*dict | sort --key=3.4bn,3 -k4.4b,4 -k2.4b,2 > "${sumFPB}_dict.csv"
  echo "Entries: " $(wc -l "${sumFPB}_dict.csv")
  cut -f2-4 "${sumFPB}_dict.csv" | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}'> "${sumFPB}_dict_uniq_SNLNM5.csv"
  echo "Entries: " $(wc -l "${sumFPB}_dict_uniq_SNLNM5.csv")
  cut -f3-4 "${sumFPB}_dict.csv" | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $1,$2,$3}'> "${sumFPB}_dict_uniq_LNM5.csv"
  echo "Entries: " $(wc -l "${sumFPB}_dict_uniq_LNM5.csv")
fi
