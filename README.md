# ebi-parasite


#### General introduction of how to use the Crypto analysis pipeline

In crypto genome analysis, the first command "download" prepares the reference file, so need to be run firstly. After that, You can start to the basis analsis (command 2-5) for each individual isolate genome. We recommend you always run command 2 "quality control" as the first step. Command 3 "assembly" is a relatively independent step, while command 2,4,5 need to be run sequentially. After all the isolate genomes go through basis analysis individually, they will be ready for the independent advanced analysis (command 6-10) as a whole. 

#### File examples

The input fastq files can be found in https://www.ebi.ac.uk/ena/data/view/PRJEB15112.
In github, you can find the example of the following files: 
 - properties.txt
 - map_hominis_genotype_A10G2.txt
 - 28Iso_map_for_visua.txt
 - hominis_go_26Jul2019.csv

#### Commands
1. [download](#download)
2. [quality_control](#quality control)
3. [assembly](#assembly)
4. [reference_mapping](#reference mapping)
5. [variation_(snp_and_indel)_call_by_using_GATK](#variation (snp and indel) call by using GATK)
6. [dNdS_analysis](#dNdS analysis)
7. [short_repeats_variation_analysis](#short repeats variation analysis)
8. [multiple_alignment_for_all_individual_chromosomes](#multiple alignment for all individual chromosomes) 
9. [variation_visualization](#variation visualization)
10. [multiQC](#multiQC)               

### download 

The scripts downloads the latest genome dna fasta and gff3 file and index fasta file:
```
download_and_index_genome_files.py -p $full_dir/properties.txt -g cryptosporidium_hominis
```
output files:
 - $workdir/ref/$genome_name+".fasta"
 - $workdir/ref/$genome_name+".gff3"
 - $workdir/ref/$genome_name+".fasta.fai"
 - $workdir/ref/$genome_name+".fasta.pac"
 - $workdir/ref/$genome_name+".fasta.bwt"
 - $workdir/ref/$genome_name+".fasta.ann"
 - $workdir/ref/$genome_name+".fasta.amb"
 - $workdir/ref/$genome_name+".fasta.sa"
 - $workdir/ref/$genome_name+".*.bt2"

### quality control

The script provides quality control on fastq files using trim_galore and remove
duplicated reads using clumpify. The quality of the original and filtered
reads are monitored by fastQC and visualized by multiQC.
```
quality_control.py -p $full_dir/properties.txt -fq1 $full_dir/ERR2889329_1.fastq -fq2 $full_dir/ERR2889329_2.fastq -pre test -de -m $full_dir/map_hominis_genotype_A10G2.txt
```
output files
 - $workdir/quality/qc/$prefix/($sample_name_)$runID+"*_fastqc.html"
 - $workdir/quality/qc/$prefix/($sample_name_)$runID+"*_fastqc.zip"
 - $workdir/quality/qc/$prefix/($sample_name_)$runID+".multiQC.html"

single_end:
 - $workdir/quality/out/$prefix_$runID.fastq_trimming_report.txt
 - $workdir/quality/out/$prefix_$runID_trimmed.fq

paired_end:
 - $workdir/quality/out/$prefix_$runID_1or2.fastq_trimming_report.txt
 - $workdir/quality/out/$prefix_$runID_1or2_val_1or2.fq

### assembly

The script assembles short reads by using spades and provides statistics summary
based on the result assemblies by using QUAST, and visualized by multiQC.
```
assembly.py -p $full_dir/properties.txt -g cryptosporidium_hominis -fq1 $full_dir/ERR2889329_1.fastq -fq2 $full_dir/ERR2889329_2.fastq -pre test -m $full_dir/map_hominis_genotype_A10G2.txt
```
output files:
 - $workdir/assembly/out/$prefix_($smple_name_)$runID+"_scaffolds.fasta"
 - $workdir/assembly/out/$prefix_($smple_name_)$runID+"_spades.log"
 - $workdir/assembly/qc/$prefix/report.txt
 - $workdir/assembly/qc/$prefix/report.tsv
 - $workdir/assembly/qc/$prefix/report.pdf
 - $workdir/assembly/qc/$prefix/report.html
 - $workdir/assembly/qc/$prefix/($sample_name_)$runID.multiQC.html

### reference mapping

The script mapping short reads to reference genomes using BWA or bowtie2, and creates 
statistics summary files using fastQC, qualiMap, and multiQC

reference mapping by using bwa:
```
reference-mapping_27Aug19.py -p $full_dir/properties.txt -t bwa -g cryptosporidium_hominis -fq1 $full_dir/ERR2889329_1_val_1.fq -fq2 $full_dir/ERR2889329_2_val_2.fq -pre xin_test_bwa_pair_ERR970586_trimmed -de -f illumina -l hominis
```
reference mapping by using bowtie2:
```
reference-mapping.py -p $full_dir/properties.txt -t bowtie2 -g cryptosporidium_hominis -fq1 $full_dir/ERR2889329_1_val_1.fq -fq2 $full_dir/ERR2889329_2_val_2.fq -pre test24Sep19 -f illumina -l hominis -m $full_dir/map_hominis_genotype_A10G2.txt
```
reference mapping by using bowtie2 for recombination
```
reference-mapping.py -p $full_dir/properties.txt -t bowtie2 -g cryptosporidium_hominis -fq1 $full_dir/ERR2889329_1_val_1.fq -fq2 $full_dir/ERR2889329_2_val_2.fq -pre test24Sep19 -f illumina -l hominis -m $full_dir/map_hominis_genotype_A10G2.txt -recom
```
output files:
 - $workdir/reference_mapping/out/$prefix_($sample_name_)$runID_*.bam
 - $workdir/reference_mapping/qc/$prefix_($sample_name_)$runID_fastqc.zip
 - $workdir/reference_mapping/qc/$prefix_($sample_name_)$runID_fastqc.html
 - $workdir/reference_mapping/qc/$prefix_($sample_name_)$runID.log
 - $workdir/reference_mapping/qc/$prefix_($sample_name_)$runID.multiQC.html

### variation (snp and indel) call by using GATK:

The script creates filtered or unfiltered SNP and INDEL vcf and gvcf files
from bam files using gatk, and then creates statistics summary by using
bcf_tools and multiQC
```
variation_gatk.py -p $full_dir/properties.txt -g cryptosporidium_hominis -bam $full_dir/ERR2889329.dedup.bam -f -pre test -m $full_dir/map_hominis_genotype_A10G2.txt
```
output files
 - $workdir/snp/out/$prefix_($sample_name_)$runID*.vcf
 - $workdir/snp/qc/$prefix/($sample_name_)$runID+".bcf_stats"
 - $workdir/snp/qc/$prefix/($sample_name_)$runID+".multiQC.html"

After all above commands are run sequencially for each isolate, the following advanced analysis can be run on all isolates independently.

### dNdS analysis

The script invests genes under selection pressure within species through dNdS. It
creates variation annotation file for each vcf file, and gene variation
annotation summary file based on all vcf files by using snpEff.
```
dNdS_within_species.py -p /$full_dir/properties.txt -g 'cryptosporidium_hominis' -vp '$full_dir/*.snp.vcf' -go $full_dir/hominis_go_26Jul2019.csv -m $full_dir/map_hominis_genotype_A10G2.txt -pre 28Iso
```
output files:
 - $workdir/dNdS/out/$prefix + "_gene_matrix.csv"
 - $workdir/dNdS/out/$prefix + "_gene_summary.csv"
 - $workdir/dNdS/out/$prefix + "_posi.csv"
 - $workdir/dNdS/qc/$prefix/$runID or $sample_$runID

### short repeats variation analysis

The script creates the genome Short Tandem Repeat (STR) variation summary file
based on all vcf files and creates multiple alignment files for all repeat regions.
```
get_repeat_variations.py -p $full_dir/properties.txt -g cryptosporidium_hominis -vp '$full_dir/*.vcf' -pre 28Iso -m $full_dir/map_hominis_genotype_A10G2.txt
```
output files
 - $workdir/repeats/out/$prefix+".summary"
 - $workdir/repeats/out/multi_align/*clustalo_num

### multiple alignment for all individual chromosomes

The script creates all individual chromosome multiple alignment for recombination.
```
build_chr_multiAlign_for_recombi.py -p $full_dir/properties.txt -g cryptosporidium_hominis -bp '$full_dir/*.bam' -pre forSimone -m $full_dir/map_hominis_genotype_A10G2.txt
```
output files:
 - $workdir/recombination/out/$prefix_$chromosome+".xmfa"

### variation visualization

a) prepare gvcf files for creating image:

The script merges gvcf files, seperating the variations into SNP and INDEL,and then do filtering if requested.
```
variation_gatk_gvcf.py -p $full_dir/properties.txt -g cryptosporidium_hominis -gv '$full_dir/*.snps.indels.g.vcf' -f -pre 28Iso
```
b) Create images

The script creates the following images: phylogeny tree, PCA, heatmap, upset, and SNP
distribution on chromosome for vcf files
```
create_all_images.py -p /$full_dir/properties.txt -g cryptosporidium_hominis -v '$full_dir/*.snp.vcf' -gv '$full_dir/*.snp.g.vcf' -gt $full_dir/28Iso_map_for_visua.txt -pre 28Iso
```
output files:
 - $workdir/visualization/out/$prefix+"_phylo_PCA_upset_heatmap.pdf"
 - $workdir/visualization/out/$prefix+"snp_dist_on"$chromosome".pdf"

### multiQC 

The script creates multiQC html file based on fastqc, bcftools, snpEff, QUAST, and
QualiMap output files. This script can only be run after the all isolate genomes go through commands 1-6.
```
multiQC.py -p $full_dir/properties.txt -pre test
```
output files:
 - $workdir/qc_report/out/$prefix+".multiQC.html"


