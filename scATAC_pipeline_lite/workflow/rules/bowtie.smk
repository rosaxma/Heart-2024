rule match_barcodes: 
    """
    Barcode correction and filtering
    """
    input: 
        fq_R1 = "fastqs/{sample}/r1.fastq.gz",
        fq_R2 = "fastqs/{sample}/r2.fastq.gz",
        fq_BC = "fastqs/{sample}/bc.fastq.adapter.removed.gz",
        whitelist ="bc_whitelist.txt",
    output: 
        fastq1_bc = "results/{sample}/bwt2/R1_bc_full.fastq.gz",
        fastq2_bc = "results/{sample}/bwt2/R2_bc_full.fastq.gz",
        qc_matching = "results/{sample}/bwt2/barcode_matching_full.tsv"
    params:
        barcode_dist = lambda w: config["max_barcode_dist"]
    threads:
        max_threads//4
    resources:
        mem_gb = 1,
	runtime_hr=24
    conda:
        "../envs/bwt2.yaml"
    script:
        "../scripts/match_barcodes.py"

rule trim_adapter:
    """
    Read adapter trimming
    """
    input:
        fastq1_bc = "results/{sample}/bwt2/R1_bc_full.fastq.gz",
        fastq2_bc = "results/{sample}/bwt2/R2_bc_full.fastq.gz",
    output:
        fastq1_trim = "results/{sample}/bwt2/R1_trim.fastq.gz",
        fastq2_trim = "results/{sample}/bwt2/R2_trim.fastq.gz",
        stats = "results/{sample}/bwt2/trim_adapters.txt"
    log:
        html = "logs/{sample}/bwt2/fastp.html",
        json = "logs/{sample}/bwt2/fastp.json"
    threads:
        max_threads//4
    resources:
        mem_gb = 256,
        runtime_hr=24
    conda:
        "../envs/bwt2.yaml"
    shell:
        "fastp -i {input.fastq1_bc} -I {input.fastq2_bc} -o {output.fastq1_trim} -O {output.fastq2_trim}"
        " -h {log.html} -j {log.json} -G -Q -L -w $(({threads} * 2)) 2> {output.stats}"

rule bowtie_index:
    """
    Build bowtie_index 
    """
    input:
        "genomes/genome.fa"
    output:
        directory("genomes/bwt2_idx")
    log:
        index = "logs/bwt_idx.log",
        index_err = "logs/bwt_idx_err.log",
    resources:
        mem_gb = 10,
        runtime_hr = 24
    conda:
        "../envs/bwt2.yaml"
    shell:
        "mkdir -p {output}; "
        "bowtie2-build {input} {output}/index "
        "> {log.index} 2> {log.index_err}"

rule bowtie2:
    """
    Read mapping (Bowtie2 aligner)
    """
    input:
        fastq1 = "results/{sample}/bwt2/R1_trim.fastq.gz",
        fastq2 = "results/{sample}/bwt2/R2_trim.fastq.gz",
        index = "genomes/bwt2_idx",
    output:
        touch("results/{sample}/bwt2/alignment.done"),
        bam_raw = "results/{sample}/bwt2/raw_collated.bam",
        qc = "results/{sample}/bwt2/bwt2_stats.txt"
    params:
        k = 1 + config["multimapping"]
    threads: 32
    resources:
        mem_gb = 20,
        runtime_hr = 24
    conda:
        "../envs/bwt2.yaml"
    shell:
        "bowtie2 -X 2000 --threads {threads} -x {input.index}/index "
        "-1 {input.fastq1} -2 {input.fastq2} --sam-append-comment -k {params.k} 2> {output.qc} | "
        "samtools view -b -S -o {output.bam_raw} -"

rule filter_multimappers:
    """
    Remove or assign multimapping reads
    """
    input:
        "results/{sample}/bwt2/raw_collated.bam"
    output:
        "results/{sample}/bwt2/primary_align.bam"
    params:
        multimapping = config["multimapping"],
        mmp_path = script_path("scripts/assign_multimappers.py")
    threads: 4
    resources:
        mem_gb = 8,
        runtime_hr=24
    conda:
        "../envs/bwt2.yaml"
    shell:
        "samtools view -F 524 -f 2 -h {input} | "
        "python {params.mmp_path} --paired-end -k {params.multimapping} | "
        "samtools view -u - | "
        "samtools fixmate -r - {output}"

rule remove_duplicates:
    """
    Mark and remove PCR duplicates
    """
    input:
        "results/{sample}/bwt2/primary_align.bam"
    output:
        bam_nodup = "results/{sample}/bwt2/alignments_map_out.bam",
        markdup_stats = "results/{sample}/bwt2/markdup.txt"
    resources:
        mem_gb = 8,
        runtime_hr=24
    log:
        "logs/{sample}/bwt2/picard.log"
    conda:
        "../envs/bwt2.yaml"
    shell:
        "picard MarkDuplicates --INPUT {input} --OUTPUT {output.bam_nodup} --METRICS_FILE {output.markdup_stats} "
        "--VALIDATION_STRINGENCY LENIENT --ASSUME_SORT_ORDER queryname --REMOVE_DUPLICATES true --BARCODE_TAG CB 2> {log}"






