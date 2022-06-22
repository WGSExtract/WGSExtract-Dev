Content of the Reference directory provided with WGS Extract v4 Installation:

UCSC Liftover Chain file (ucsc.edu):
-- used by pyLoftover in the Microarray generator when starting with Build 38 model
* hg38ToHg19.over.chain.gz
  (and corresponding .gzi index file)

Y Chromosome SNPs file (yBrowse.org):
-- used by Y DNA VCF file generator
* snps_hg38.vcf.gz, snps_grch38.vcf.gz -- Build 38 files in chrN and N nomenclature; respectively
* snps_hg19.vcf.gz, snps_grch38.vcf.gz -- Build 37 files in chrN and N nomenclature; respectively
   (and corresponding .gzi and .tbi index files)
   (note, original hg19 had some errors and have not been updated since 2017. Corrected. .orig original there)

Exome BED Files for determing WES regions ():
-- Original TruSeq from Illumina for Build 37; converted form for Build38 from BioBank
        https://support.illumina.com/downloads/truseq-exome-product-files.html
        https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=3803
* TruSeq_Exome_TargetedRegions_v1.2.bed -- Build 37 chrN nomenclature
* TruSeq_Exome_TargetedRegions_v1.2_GRCh.bed -- Build 37 N nomenclature
* xgen_plus_spikein.GRCh38.bed, xgen_plus_spikein.GRCh38.GRCh.bed -- Build 38 chrN and N nomenclature

Y Chromosome BED files (with added mito) for use in the WES buttons on Y (or Y and MT) only BAM files.
Instead of using the Y Exome BED in those cases, will use the Phylogenetic Tree of Haplogroup
communities regions of interest.  Original from a spreadsheet of HG38 regions from David Vance. See
manual for more details.

Genomes:
-- See the shell scripts for the 10 reference genomes that are downloaded and indexed at install time
* process_reference_genomes.sh -- sript to index any available reference genomes (called by Upgrade)
* compare_reference_genomes.sh -- script to compare two reference genomes based on extracted stats
* get_and_process_*.sh -- separate script for each reference genome to individually download and process
* bwa-kit-gen-ref.sh   -- modified (updated, corrected) script from bwa-kit to make 1K Genome project reference models

Microarray:
-- All Variants file to create CombinedKit from BAM; one for each build and naming style
-- Microarray RAW file templates for headers and bodies to target generating microazray file types
-- Ploidy definition file
