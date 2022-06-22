#!/usr/bin/env bash
#
# Merger of function previously in end of Upgrade_UbuntuLinux.sh, get_reference_genomes.sh and the 14 individual
#  get_and_process_refgenname.sh files. Simplified Upgrade to simply call this script.  Any presented selection
#  menu comes from here.  Needed to provide one-click, command window control for users struggling with reference
#  genome downloads. Some, even in North America, are showing issues with accessing the NIH server. Most still seem
#  to operate fine with the typical "All" that uses the NIH server for the 1K genome models.
#
# Todo Change to a Python function and make a Reference Library manager tab in the main program.
#
# Part of the Reference Genome package in WGS Extract (https://wgsextract.github.io/)
# Copyright (c) 2021-22 Randy Harr
#

if [[ $# -ne 1 ]] ; then
  printf "Usage: %s <path to WGSE>" "$0"
  printf "  Requests input from user on Reference Genomes to download and process.\n"
  printf "  Needs the path to the WGS Extract installation as the Reference Library may have been relocated.\n"
  exit
fi
WGSEFIN=$1

# MacOS is overriding PATH and using builtin (ancient) BASH; so need absolute path override.
# Make sure we have the bioinformatic tools on the path
case $OSTYPE in
  linux*)
    bashx="/usr/bin/env bash"
    home=~  ;;
  darwin*)
    bashx="/opt/local/bin/bash"
    home=~
    if [ -z "$(echo $PATH | grep "/opt/local/bin")" ]; then
      PATH="/opt/local/bin:${PATH}"
    fi   ;;
  msys*|cygwin*)
    bashx="/bin/bash.exe"       ;
    homew=$(env | grep USERPROFILE | cut -d '=' -f2)
    home=$(cygpath -u "${homew}")
    if [ -z "$(echo $PATH | grep "/usr/local/bin")" ]; then
      PATH="/usr/local/bin:/bin:${PATH}"
    fi   ;;
  *)  printf "*** Error: unknown OSTYPE of %s\n" "$OSTYPE" &&  exit 1  ;;
esac

curlx="curl -kLZC - --retry 5"

newreflib="null"                        # If command below fails; then returns 0/"null"
if [[ -e "${home}/.wgsextract" ]]; then
  newreflib=$(cat ${home}/.wgsextract | jq -r '."reflib.FP"')   # Return string from settings (else 0/"null")
fi
reflib="${WGSEFIN}/reference/"          # Default location in installation directory
if [[ "$newreflib" != "null" ]]; then
  reflib="${newreflib}"                 # From settings file -- reference library was moved
fi

# Moved out of separate get_reference_genomes.sh script to here to simply keep duplicate curl code together
# Parameterized to allow NIH or EBI source selection as an option
get_all_refgenomes() {
  # Start with 1K Genome models so users likely to see failures right up front. Are the recommended anyway.
  case $1 in
    "NIH")
      [ ! -f hs38.fa.gz ] && echo Grabbing hs38 NIH && \
        ${curlx} -o hs38.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
      [ ! -f hs37d5.fa.gz ] && echo Grabbing hs37d5 NIH && \
        ${curlx} -o hs37d5.fa.gz ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
      [[ ! -f hs38DH.fa.gz && ! -f hs38DH.fa ]] && echo Grabbing hs38DH NIH && \
        ${curlx} -o hs38DH.fa ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
      ;;
    "EBI")    # Actually hs38 is not available on the EBI server so we supply a WGSE Google Drive version
      [ ! -f hs38.fa.gz ] && echo Grabbing hs38 EBI && \
        ${curlx} -o hs38.fa.gz 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjgQbu7wjUVzD86txD/root/content'
      [ ! -f hs37d5.fa.gz ] && echo Grabbing hs37d5 EBI && \
        ${curlx} -o hs37d5.fa.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
      [[ ! -f hs38DH.fa.gz && ! -f hs38DH.fa ]] && echo Grabbing hs38DH EBI && \
        ${curlx} -o hs38DH.fa ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
      ;;
    *)
      echo "*** Unknown server selection: $1"
      exit
      ;;
  esac
  [ ! -f chm13v2.0.fa.gz ] && echo "Grabbing T2T v2 chm13 with Y - chrN naming" && \
    ${curlx} -o chm13v2.0.fa.gz 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz'
  # UCSC server for HG models tends to be slow but is still reliable
  [ ! -f hg38.fa.gz ] && echo Grabbing hg38 && \
    ${curlx} -o hg38.fa.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  [ ! -f hg19.fa.gz ] && echo Grabbing hg19 yoruba && \
    ${curlx} -o hg19.fa.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
  # ySeq is the only source of an HG19 model with an rCRS mitochondrial model (what we call hg37 here)
  [[ ! -f hg19_yseq.fa.gz && ! -f hg19_yseq.fa.zip ]] && echo "Grabbing hg37 rCRS yseq" && \
    ${curlx} -o hg19_yseq.fa.zip http://genomes.yseq.net/WGS/ref/hg19/hg19.zip
  # Unique EBI models we likely will not encounter but ....
  [ ! -f Homo_sapiens.GRCh37.dna.toplevel.fa.gz ] && echo Grabbing Homo_sapiens.GRCh37.dna.toplevel && \
    ${curlx} -o Homo_sapiens.GRCh37.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
  [ ! -f Homo_sapiens.GRCh38.dna.toplevel.fa.gz ] && echo Grabbing Homo_sapiens.GRCh38.dna.toplevel && \
    ${curlx} -o Homo_sapiens.GRCh38.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
}


# This selection and the options All (NIH) or None is what was in the original v3 release) We have added the other
# options and EBI server processing to help those having trouble downloading from the NIH servers. More in EU have the
# issue with the NIH server than elsewhere.  But even US based people have issues with access to NIH (rarely)).
# Suspect it is jitter / latency; especially on a wireless connections.  But not enough data to characterize it yet.
echo
echo --------------------------------------------------------------------------------
echo WGS Extract Reference Genome Library Installation and Update
echo --------------------------------------------------------------------------------
PS3='Choose which Reference Genome(s) to process now (1 to Exit): '
# Main group, controlling options. All terminate program after being selected and run.
option=("Exit" "Recommended (@US NIH)" "Recommended (@EU EBI)" "First 9 (@US NIH)" "First 9 (@EU EBI)" )
# Recommended 3 individually and first of First 9
option+=("hs38 (Nebula) (R)" "hs37d5 (Dante) (R)" "T2T_v2 (PGP/HPP chrN) (R)" "hs38DH (aDNA)" )
# Rest of First 9 (was traditionally ALL in v3 release; v2 had 5 fixed and delivered with the tool)
option+=("hg38 (ySeq)" "hg37 (rCRS, ySeq)" "hg19 (yoruba, ucsc)" "GRCh38 (@Ensembl)" "GRCh37 (@Ensembl)" )
# EBI versions from Recommended and First 9; T2T v2 Genbank original with accession naming
option+=("hs38 (@EU EBI)" "hs37d5 (@EU EBI)" "hs38DH (@EU EBI)" "T2T_v2 (PGP/HPP Genbank)" )
# Additional oddball T2T / HPP / PGP models that should be re-aligned to T2Tv2
option+=("hg002xy_v2.7 (T2T)" "hg002xy_v2 (T2T)" "chm13y_v1.1 (HPP)" "chm13y_v1 (HPP)")
# Human_g1k model from original release and 1K Genome Project; was in original v2 release (use hs37d5 instead)
option+=("human_g1k_v37 (@EU EBI)" "human_g1k_v37 (@US NIH)")
# Unique, custom, not used anywhere else models (use to load but immediately re-align to Recommended)
option+=("hs38s (by Sequencing)" "hg38+CP086569 (by ySeq)" "hg19_wgse (in v1)")

#  We only redisplay the menu after every 3 downloads; when menu will likely scroll off the screen. So keep track.
declare -i menu_cnt
menu_cnt=0

select rg in "${option[@]}"; do
  case $rg in
    "Exit")
      echo "Exiting. You can run the Library script later to add or update the reference genomes."
      echo ""
      break
      ;;

    "Recommended (@US NIH)")
      echo "You can start the WGS Extract program while the Reference Library is downloading."
      echo "When selecting Recommended, we always download and replace any existing files."
      [ ! -f hs37d5.fa.gz ] && echo "Downloading hs37d5 (default NIH server)" && \
        ${curlx} -o hs37d5.fa.gz ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
      [ ! -f hs38.fa.gz ] && echo "Downloading hs38 (default NIH server)" && \
        ${curlx} -o hs38.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
      [ ! -f chm13v2.0.fa.gz ] && echo "Downloading T2T v2.0: chm13 v1.1 with HG002 Y v2.7 (T2T AWS)" && \
        ${curlx} -o chm13v2.0.fa.gz 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz'
      ${bashx} process_reference_genomes.sh hs37d5.fa.gz hs38.fa.gz chm13v2.0.fa.gz
      echo "Finished with Recommended (@US NIH)."
      menu_cnt+=2
      ;;
    "Recommended (@EU EBI)")
      echo "You can start the WGS Extract program while the Reference Library is downloading."
      echo "When selecting Recommended, we always download and replace any existing files."
      [ ! -f hs37d5.fa.gz ] && echo "Downloading hs37d5 (EBI server)" && \
        ${curlx} -o hs37d5.fa.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
      [ ! -f hs38.fa.gz ] && echo "Downloading hs38 (from WGS Extract Google Drive; EBI copy does not exist)" && \
        ${curlx} -o hs38.fa.gz 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjgQbu7wjUVzD86txD/root/content'
      [ ! -f chm13v2.0.fa.gz ]  && echo "Downloading T2T v2.0: chm13 v1.1 with HG002 Y v2.7 (T2T AWS)" && \
        ${curlx} -o chm13v2.0.fa.gz 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz'
      ${bashx} process_reference_genomes.sh hs37d5.fa.gz hs38.fa.gz chm13v2.0.fa.gz
      echo "Finished with Recommended (@EU EBI)."
      menu_cnt+=2
      ;;

    "First 9 (@US NIH)")
      echo "You can start the WGS Extract program while the Reference Library is downloading."
      echo "When selecting First 9, we only download any missing files; even if the existing ones are bad or partial"
      # Compared to v3 ALL, removed human_g1k_v37 and hg19_wgsev1; added T2Tv2)
      get_all_refgenomes NIH              # Download any missing reference genomes from their (NIH) source
      ${bashx} process_reference_genomes.sh . # Fix compression and create indices and support files for each genome
      echo "Finished with First 9 (@US NIH)."
      menu_cnt+=8
      ;;
    "First 9 (@EU EBI)")
      echo "You can start the WGS Extract program while the Reference Library is downloading."
      echo "When selecting First 9, we only download any missing files; even if the existing ones are bad or partial"
      # Compared to v3 ALL, removed human_g1k_v37 and hg19_wgsev1; added T2Tv2)
      get_all_refgenomes EBI              # Download any missing reference genomes from their (EBI) source
      ${bashx} process_reference_genomes.sh . # Fix compression and create indices and support files for each genome
      echo "Finished with First 9 (@EU EBI)."
      menu_cnt+=8
      ;;

    # Individual models covered by the Recommended (3) and First 9 (6 additional) commands (NIH and EBI source)
    # Compared to v3 ALL, removed human_g1k_v37 and hg19_wgse; added T2Tv2)
    "hs38 (Nebula) (R)")
      echo "Downloading and Processing hs38 (default NIH server)"
      [ -f hs38.fa.gz ] && rm -f hs38.fa.gz
      ${curlx} -o hs38.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
      [ -f hs38.fa.gz ] && ${bashx} process_reference_genomes.sh hs38.fa.gz
      ;;
    "hs37d5 (Dante) (R)")
      echo "Downloading and Processing hs37d5 (default NIH server)"
      [ -f hs37d5.fa.gz ] && rm -f hs37d5.fa.gz
      ${curlx} -o hs37d5.fa.gz ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
      [ -f hs37d5.fa.gz ] && ${bashx} process_reference_genomes.sh hs37d5.fa.gz
      ;;
    "T2T_v2 (PGP/HPP chrN) (R)")
      echo "Downloading and Processing T2T v2.0: chm13 v1.1 with HG002 Y v2.7 using UCSC naming (T2T AWS)"
      [ -f chm13v2.0.fa.gz ] && rm -f chm13v2.0.fa.gz
      ${curlx} -o chm13v2.0.fa.gz 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz'
      [ -f chm13v2.0.fa.gz ] && ${bashx} process_reference_genomes.sh chm13v2.0.fa.gz
      ;;

    "hs38DH (aDNA)")
      echo "Downloading and Processing hs38DH (default NIH serverl uncompressed original so 3x the size)"
      [ -f hs38DH.fa.gz ] && rm -f hs38DH.fa.gz
      ${curlx} -o hs38DH.fa ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
      [ -f hs38DH.fa ] && ${bashx} process_reference_genomes.sh hs38DH.fa
      ;;
    "hg38 (yseq)")
      echo "Downloading and Processing hg38 (UCSC) (used by ySeq)"
      [ -f hg38.fa.gz ] && rm -f hg38.fa.gz
      ${curlx} -o hg38.fa.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
      [ -f hg38.fa.gz ] && ${bashx} process_reference_genomes.sh hg38.fa.gz
      ;;
    "hg37 (rCRS, ySeq)")
      echo "Downloading and Processing hg19_yseq (rCRS model) (used by ySeq)"
      [ -f hg19_yseq.fa.gz ] && rm -f hg19_yseq.fa.gz     # Post-processed version
      [ -f hg19_yseq.fa.zip ] && rm -f hg19_yseq.fa.zip   # Downloaded version, before processing
      ${curlx} -o hg19_yseq.fa.zip http://genomes.yseq.net/WGS/ref/hg19/hg19.zip
      [ -f hg19_yseq.fa.zip ] && ${bashx} process_reference_genomes.sh hg19_yseq.fa.zip
      ;;
    "hg19 (yoruba, ucsc)")
      echo "Downloading and Processing hg19 (Yoruba, UCSC) (used by Dante in 2018)"
      [ -f hg19.fa.gz ] && rm -f hg19.fa.gz
      ${curlx} -o hg19.fa.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
      [ -f hg19.fa.gz ] && ${bashx} process_reference_genomes.sh hg19.fa.gz
      ;;

    "GRCh38 (@Ensembl)")
      echo "Downloading and Processing Homo_sapiens.GRCh38 (EBI Ensembl)"
      [ -f Homo_sapiens.GRCh38.dna.toplevel.fa.gz ] && rm -f Homo_sapiens.GRCh38.dna.toplevel.fa.gz
      ${curlx} -o Homo_sapiens.GRCh38.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
      [ -f Homo_sapiens.GRCh38.dna.toplevel.fa.gz ] && ${bashx} process_reference_genomes.sh Homo_sapiens.GRCh38.dna.toplevel.fa.gz
      ;;
    "GRCh37 (@Ensembl)")
      echo "Downloading and Processing Homo_sapiens.GRCh37 (EBI Ensembl)"
      [ -f Homo_sapiens.GRCh37.dna.toplevel.fa.gz ] && rm -f Homo_sapiens.GRCh37.dna.toplevel.fa.gz
      ${curlx} -o Homo_sapiens.GRCh37.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
      [ -f Homo_sapiens.GRCh37.dna.toplevel.fa.gz ] && ${bashx} process_reference_genomes.sh Homo_sapiens.GRCh37.dna.toplevel.fa.gz
      ;;

    "hs38 (@EU EBI)")
      echo "Downloading and Processing hs38 (WGSE MS Onedrive; EBI does not exist)"
      [ -f hs38.fa.gz ] && rm -f hs38.fa.gz
      ${curlx} -o hs38.fa.gz 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjgQbu7wjUVzD86txD/root/content'
      [ -f hs38.fa.gz ] && ${bashx} process_reference_genomes.sh hs38.fa.gz
      ;;
    "hs37d5 (@EU EBI)")
      echo "Downloading and Processing hs37d5 (EBI server)"
      [ -f hs37d5.fa.gz ] && rm -f hs37d5.fa.gz
      ${curlx} -o hs37d5.fa.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
      [ -f hs37d5.fa.gz ] && ${bashx} process_reference_genomes.sh hs37d5.fa.gz
      ;;
    "hs38DH (@EU EBI)")
      echo "Downloading and Processing hs38DH (EBI server, uncompressed original so 3x the size)"
      [ -f hs38DH.fa.gz ] && rm -f hs38DH.fa.gz
      [ -f hs38DH.fa ] && rm -f hs38DH.fa
      ${curlx} -o hs38DH.fa ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
      [ -f hs38DH.fa ] && ${bashx} process_reference_genomes.sh hs38DH.fa
      ;;

    # 10 Secondary models only downloaded by individual request here (not part of Recommended or First 9)
    # Should only be used if required to process an existing BAM / CRAM and transform into a recommended.
    "T2T_v2 (PGP/HPP Genbank)")
      echo "Downloading and Processing the T2T v2.0 (Genbank accession named; UCSC file name)"
      [ -f GCA_009914755.4.fa.gz ] && rm -f GCA_009914755.4.fa.gz
      ${curlx} -o GCA_009914755.4.fa.gz 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz'
      [ -f GCA_009914755.4.fa.gz ] && ${bashx} process_reference_genomes.sh GCA_009914755.4.fa.gz
      ;;
    "hg002xy_v2.7 (T2T)")
      echo "Downloading and Processing T2T chm13 v1.1 with HG002 XY v2.7 (WGSE MS OneDrive)"
      [ -f hg002xy_v2.7.fasta.gz ] && rm -f hg002xy_v2.7.fasta.gz
      ${curlx} -o hg002xy_v2.7.fasta.gz 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjgRHW0-TxY1_mhGV-/root/content'
      [ -f hg002xy_v2.7.fasta.gz ] && ${bashx} process_reference_genomes.sh hg002xy_v2.7.fasta.gz
      ;;
    "hg002xy_v2 (T2T)")
      echo "Downloading and Processing T2T chm13 v1.1 with HG002 XY v2 (WGSE MS OneDrive)"
      [ -f hg002xy_v2.fasta.gz ] && rm -f hg002xy_v2.fasta.gz
      ${curlx} -o hg002xy_v2.fasta.gz 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjgRB_ilxljbLqzCX6/root/content'
      [ -f hg002xy_v2.fasta.gz ] && ${bashx} process_reference_genomes.sh hg002xy_v2.fasta.gz
      ;;
    "chm13y_v1.1 (HPP)")
      echo "Downloading and Processing HPP chm13 T2T v1.1 with GRCh38 Y"
      [ -f CHM13v11Y.fa.gz ] && rm -f CHM13v11Y.fa.gz
      ${curlx} -o CHM13v11Y.fa.gz 'https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph/CHM13v11Y.fa.gz'
      [ -f CHM13v11Y.fa.gz ] && ${bashx} process_reference_genomes.sh CHM13v11Y.fa.gz
      ;;
    "chm13y_v1 (HPP)")
      echo "Downloading and Processing HPP chm13 T2T v1 with GRCh38 Y"
      [ -f CHM13v1Y.fa.gz ] && rm -f CHM13v1Y.fa.gz
      ${curlx} -o CHM13v1Y.fa.gz 'https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph/CHM13v11Y.fa.gz'
      [ -f CHM13v1Y.fa.gz ] && ${bashx} process_reference_genomes.sh CHM13v1Y.fa.gz
      ;;

    "human_g1k_v37 (@US NIH)")
      echo "Downloading and Processing human_g1k_v37"
      [ -f human_g1k_v37.fasta.gz ] && rm -f human_g1k_v37.fasta.gz
      ${curlx} -o human_g1k_v37.fasta.gz ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
      [ -f human_g1k_v37.fasta.gz ] && ${bashx} process_reference_genomes.sh human_g1k_v37.fasta.gz
      ;;
    "human_g1k_v37 (@EU EBI)")
      echo "Downloading and Processing human_g1k_v37 (EBI server)"
      [ -f human_g1k_v37.fasta.gz ] && rm -f human_g1k_v37.fasta.gz
      ${curlx} -o human_g1k_v37.fasta.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
      [ -f human_g1k_v37.fasta.gz ] && ${bashx} process_reference_genomes.sh human_g1k_v37.fasta.gz
      ;;

    "hs38s (by Sequencing)")
      # Sequencing.com unique hs38 model (GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz with 22_KI270879v1_alt from hs38DH added in)
      echo "Downloading and Processing hs38s (Sequencing.com) (WGSE MS OneDrive)"
      [ -f hs38s.fa.gz ] && rm -f hs38s.fa.gz
      ${curlx} -o hs38s.fa.gz 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjgR0QualUlHx53-0U/root/content'
      [ -f hs38s.fa.gz ] && ${bashx} process_reference_genomes.sh hs38s.fa.gz
      ;;
    "hg38+CP086569 (by ySeq)")
      echo "Downloading and Processing T2T yseq (hg38 with hg002y v2.0)"
      [ -f hg38_CP086569.fasta.gz ] && rm -f hg38_CP086569.fasta.gz     # Post-processed version
      [ -f hg38_CP086569.fasta ] && rm -f hg38_CP086569.fasta   # Downloaded version, before processing
      ${curlx} -o hg38_CP086569.fasta http://genomes.yseq.net/WGS/ref/CP086569.1/hg38_CP086569.fasta
      [ -f hg38_CP086569.fasta.gz ] && ${bashx} process_reference_genomes.sh hg38_CP086569.fasta
      ;;
    "hg19_wgse (in v1)")
      # Was in original WGS Extract v1 and v2 release; unique. Include simply for historical reasons. Should not be used.
      echo "Downloading and Processing hg19_wgse (custom) (WGSE MS OneDrive)"
      [ -f hg19_wgse.fa.gz ] && rm -f hg19_wgse.fa.gz
      ${curlx} -o hg19_wgse.fa.gz 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjde1tXQ1bDvCrlFo/root/content'
      [ -f hg19_wgse.fa.gz ] && ${bashx} process_reference_genomes.sh hg19_wgse.fa.gz
      ;;
    *)
      echo "Please enter a valid option (emter 1 to exit without further processing)"
      ;;
  esac

  menu_cnt+=1
  if (( menu_cnt > 3 )); then
    REPLY=""      # Causes the menu to be redisplayed again before the prompt
    menu_cnt=0
  fi

done
echo "Finished installing and processing Reference Genomes."
