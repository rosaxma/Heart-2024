rule filter_mito:
    """
    Filter and count mitochondrial reads (and also fiter out secondary alignments)
    """
    input: 
        lambda w: f"results/{w.sample}/{config['aligner']}/alignments_map_out.bam"
    output: 
        bam = "results/{sample}/alignments_no_mito.bam",
        qc = "results/{sample}/frac_mito.tsv"
    params:
        mitochr = lambda w: config["mito_chr"]
    resources:
        mem_gb = 1
    conda:
        "../envs/fragments.yaml"
    script:
        "../scripts/filter_mito.py"

rule sort_alignments:
    """
    Sort alignments
    """
    input:
        "results/{sample}/alignments_no_mito.bam"
    output:
        "results/{sample}/alignments_sorted.bam"
    threads:
        max_threads//8
    resources:
        runtime_min = 8,
        mem_gb = 16
    conda:
        "../envs/fragments.yaml"
    shell:
        "samtools sort -@ $(({threads} * 2)) -o {output} {input}"

rule bam_to_fragments: 
    """
    Convert BAM to fragment file
    """
    input:
        "results/{sample}/alignments_sorted.bam"
    output:
        "results/{sample}/fragments.tsv"
    params:
        shift_plus = config["tn5_shift_plus"],
        shift_minus = config["tn5_shift_minus"]
    resources:
        mem_gb = 1
    conda:
        "../envs/fragments.yaml"
    script:
        "../scripts/bam_to_fragments.py"

rule compress_fragments:
    """
    Compress fragment file
    """
    input:
        "results/{sample}/fragments.tsv"
    output: 
        "results/{sample}/fragments.tsv.gz"
    resources:
        mem_gb = 1
    conda:
        "../envs/fragments.yaml"
    shell: 
        "bgzip -c {input} > {output}"

rule index_fragments:
    """
    Index fragment file
    """
    input: 
        "results/{sample}/fragments.tsv.gz"
    output: 
        "results/{sample}/fragments.tsv.gz.tbi"
    resources:
        mem_gb = 1
    conda:
        "../envs/fragments.yaml"
    shell: 
        "tabix --zero-based --preset bed {input}"