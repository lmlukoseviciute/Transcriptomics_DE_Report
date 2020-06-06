# Quality control

rule fastqc:
    input:
        "inputs/{sample}.fastq.gz"
    output:
        "outputs/fastqc_outputs/raw/{sample}_fastqc.html"
    shell:
        "fastqc {input} -o ./outputs/fastqc_outputs/raw/"

rule multiqc:
    input:
        "outputs/fastqc_outputs/raw/"
    output:
        "outputs/multiqc_outputs/raw/multiqc_report.html"
    shell:
        "multiqc {input} -o outputs/multiqc_outputs/raw"

rule bbduk:
    input:
        "inputs/{sample}.fastq.gz"
    output:
        "outputs/bbduk_outputs/{sample}_trimmed.fastq.gz"
    shell:
        "bbduk.sh in={input} out={output} ref=inputs/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"

rule fastqc_trimmed:
    input:
        "outputs/bbduk_outputs/{sample}.fastq.gz"
    output:
        "outputs/fastqc_outputs/trimmed/{sample}_fastqc.html"
    shell:
        "fastqc {input} -o ./outputs/fastqc_outputs/trimmed/"

rule multiqc_trim:
    input:
        "outputs/fastqc_outputs/trimmed/"
    output:
        "outputs/multiqc_outputs/trimmed/multiqc_report.html"
    shell:
        "multiqc {input} -o outputs/multiqc_outputs/trimmed/"

# Reads alligment

rule genome_indexing:
    input:
        ref_genome="inputs/chr19_20Mb.fa",
        gtf_file="inputs/chr19_20Mb.gtf"
    output:
        "outputs/indexing/SAindex"
    shell:
        "./bins/check_dir.py 'outputs/indexing';  STAR --runThreadN 4 "
        "--runMode genomeGenerate --genomeDir outputs/indexing "
        "--genomeFastaFiles {input.ref_genome} --sjdbGTFfile {input.gtf_file} "
        "--sjdbOverhang 150 --genomeSAindexNbases 11"

rule read_alligment:
    input:
        "outputs/bbduk_outputs/{sample}_R1_001_trimmed.fastq.gz",
        "outputs/bbduk_outputs/{sample}_R2_001_trimmed.fastq.gz",
        "outputs/indexing"
    output:
        "outputs/alligments/{sample}/Aligned.sortedByCoord.out.bam"
    shell:
        "./bins/check_dir.py outputs/alligments/{wildcards.sample}/; STAR --runThreadN 4 "
        "--outFileNamePrefix outputs/alligments/{wildcards.sample}/ "
        "--genomeDir {input[2]} --readFilesCommand zcat --runMode alignReads "
        "--readFilesIn {input[0]} {input[1]} "
        "--outSAMtype BAM SortedByCoordinate"

rule bam_indexing:
    input:
        "outputs/alligments/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "outputs/alligments/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"
rule feature_counts:
    input:
        gtf_file="inputs/chr19_20Mb.gtf",
        bam_file="outputs/alligments/{sample}/Aligned.sortedByCoord.out.bam"
        
    output:
        "outputs/counts/{sample}_counts.txt"
    shell:
        "featureCounts -p -t exon -g gene_id -a {input.gtf_file} -o {output} {input.bam_file}"

rule extract_counts:
    input:
        "outputs/counts/{sample}_counts.txt"
    output:
        "outputs/counts/{sample}_counts_extracted.txt"
    shell:
        """awk '{{print $1 " " $7}}' {input} | sed '1d'| sed '1d' > {output}"""

rule de_analyses:
    input:
        "inputs/colnames.csv",
        "inputs/counts_colli.csv",
        "inputs/counts_kapa.csv"
    output:
        "outputs/plots/venn_diag.png"
    shell:
        "./bins/deseq_rscript.R"
