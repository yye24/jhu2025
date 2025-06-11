#!/bin/bash

### this script is modified from https://github.com/mingming-cgz/ARTR-seq/
## trimming
# {}_R1_001.fastq.gz & {}_R2_001.fastq.gz as raw PE reads
ls *_R1_001.fastq.gz \
| awk -F "_R1_001.fastq.gz" '{print $1}' \
| sort -u \
| parallel "cutadapt -j $n_threads --nextseq-trim=20 --action=trim \
		-a AGATCGGAAGAGCACACGTCTGAACTCCAG \
		-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAG \
		-o {}_R1_001.trim.fastq.gz -p {}_R2_001.trim.fastq.gz {}_R1_001.fastq.gz {}_R2_001.fastq.gz"

## extracting umi
ls *_R1_001.trim.fastq.gz \
| awk -F "_R1_001.trim.fastq.gz" '{print $1}' \
| sort -u \
| parallel "cutadapt -j $n_threads -q 20 -m 20 --action=trim \
		-u 8 -u -4 -U -8 -U 4 \
		--rename='{id}_{r1.cut_prefix} {comment}' \
		-o {}_R1_001.trimMI.fastq.gz -p {}_R2_001.trimMI.fastq.gz {}_R1_001.trim.fastq.gz {}_R2_001.trim.fastq.gz"

## mapping to rRNA using bowtie2
total_files=`find -name '*.gz' | wc -l`
arr=( $(ls *.gz) )
for ((i=0; i<$total_files; i+=2))
{
    sample_name=$(echo ${arr[$i]} | sed 's/_R[12].*//')
    echo "[bowtie2 mapping running for sample] $sample_name"

    date && time bowtie2 --threads $n_threads --seedlen=15 \
        -x hg38_rRNA \
        -1 ${arr[$i]} -2 ${arr[$i+1]} \
        --un-conc-gz ${sample_name}_norRNA

    printf "\n\n"
} > /dev/null

## re-mapping norRNA reads using STAR
for r1_file in *_R1.fastq.gz; do
r2_file="${r1_file/_R1.fastq.gz/_R2.fastq.gz}"
sample_name=$(basename "$r1_file" | sed 's/_R1.fastq.gz//')
output_dir="STAR/${sample_name}"
mkdir -p "$output_dir"
STAR --runMode alignReads --runThreadN $n_threads \
	--readFilesCommand zcat \
	--genomeDir refgenome/GRCh38 \
	--alignEndsType EndToEnd \
	--genomeLoad NoSharedMemory \
	--quantMode TranscriptomeSAM \
	--alignMatesGapMax 15000 \
	--readFilesIn "$r1_file" "$r2_file" \
	--outFileNamePrefix "${output_dir}/" \
	--outFilterMultimapNmax 1 \
	--outSAMattributes All \
	--outSAMtype BAM SortedByCoordinate \
	--outFilterType BySJout \
	--outReadsUnmapped Fastx \
	--outFilterScoreMin 10 \
	--outFilterMatchNmin 24 \
	> "${output_dir}/STAR_log.txt" 2>&1
done

## dedup using UMItools
ls *.bam | parallel samtools index -@ $n_threads '{}'
find . -name "*.sortedByCoord.out.bam" | parallel -j $n_threads ' \
    filename=$(basename {} .sortedByCoord.out.bam); \
    umi_tools dedup --method unique \
    -I {} \
    --output-stats=${filename}_dedup \
    -L ${filename}_umitools.log \
    --temp-dir=dedup_temp \
    -S ${filename}.sort.dedup.bam'

## splitting strands
for file in *_L001_norRNA_Aligned.sort.dedup.bam
do filename="${file%%.*}"
samtools view -@ $n_threads -f 16 $file -b -o ${filename}.sort.dedup.reverse.bam
samtools view -@ $n_threads -F 16 $file -b -o ${filename}.sort.dedup.forward.bam
done
ls *.bam | parallel samtools index -@ $n_threads '{}'

## call peaks using macs3: do fwd and rev seperately
macs3 callpeak --treatment IP_1.bam IP_2.bam IP_3.bam\
       	--control Input_1.bam Input_2.bam Input_3.bam\
	-f BAM -n $sample_name -g hs -B \
    --keep-dup all \
	--nomodel --extsize 30

## combine peaks
for fwd_file in *_fwd_peaks.narrowPeak; do
    # Derive the sample name by removing the "_fwd_peaks.narrowPeak" suffix
    sample_name="${fwd_file%_fwd_peaks.narrowPeak}"

    # Find the corresponding reverse file
    rev_file="${sample_name}_rev_peaks.narrowPeak"

    # Check if the reverse file exists
    if [[ -f "$rev_file" ]]; then
        # Combine forward and reverse files
        combined_file="${sample_name}_combined_peaks.narrowPeak"

        echo "Combining $fwd_file and $rev_file into $combined_file"

        # Annotate strands and combine
        awk 'BEGIN{FS="\t"; OFS="\t"} {$6 = "+"; print $0}' "$fwd_file" > "$combined_file"
        awk 'BEGIN{FS="\t"; OFS="\t"} {$6 = "-"; print $0}' "$rev_file" >> "$combined_file"
    else
        echo "Warning: Reverse file for $sample_name not found"
    fi
done

## making a ref peak saf file
cat *.narrowPeak > all_peaks.narrowPeak
sort -k1,1 -k2,2n -k6,6 all_peaks.narrowPeak > all_peaks_sorted.narrowPeak
bedtools merge -i all_peaks_sorted.narrowPeak -s -c 6 -o distinct > all_peaks_merged.bed
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
{ print "Peak" NR, $1, $2, $3, $4 }' all_peaks_merged.bed > all_peaks_merged.saf

## count peaks using subread
featureCounts -F SAF -a all_peaks_merged.saf -o peak_counts_s1.txt -s 1 -T $n_threads -p -B -C *.bam
featureCounts -F SAF -a all_peaks_merged.saf -o peak_counts_s2.txt -s 2 -T $n_threads -p -B -C *.bam
