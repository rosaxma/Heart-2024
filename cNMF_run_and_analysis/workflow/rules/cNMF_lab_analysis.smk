rule analysis:
    input:
        spectra_consensus=os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "{sample}_allK.spectra.k_{k}.dt_" + density + ".consensus.txt"),
        spectra_tpm=os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "{sample}_allK.gene_spectra_tpm.k_{k}.dt_" + density + ".txt"),
        spectra_score=os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "{sample}_allK.gene_spectra_score.k_{k}.dt_" + density + ".txt"),
    output:
        cnmf_results=os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis", "cNMF_results.k_{k}.dt_"+density+".RData"),
    params:
        outdir = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis"), 
        organism = config["lab_analysis"]["organism"],
        threshold = density,
        scriptdir=config["ScriptDir"]
    resources:
        mem_gb=16,
        runtime_hr=3
     conda:
        config["R_env"],
    shell:
        """
        Rscript {params.scripts}/cNMF_analysis.R \
        --spectra_consensus {input.spectra_consensus} \
        --spectra_tpm {input.spectra_tpm} \
        --spectra_score {input.spectra_score} \
        --sampleName {wildcards.sample} \
        --outdir {params.outdir} \
        --K {wildcards.k} \
        --density.thr {params.threshold} \
        --organism {params.organism}
        """


rule clusterProfiler_GSEA:
    input:
        cnmf_results=os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis", "cNMF_results.k_{k}.dt_"+density+".RData"),
    output:
        clusterProfiler_result = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis", "GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.txt")
    params:
        outdir = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis", "cluster_profiler"), 
        organism = config["lab_analysis"]["organism"],
        threshold = density,
        topGenes2look = config["lab_analysis"]["topGenes2look"],
        scriptdir=config["ScriptDir"],
    resources:
        mem_gb=16,
        runtime_hr=3
    conda:
        config["R_env"],
    shell:
        """
        Rscript {params.scripts}/cNMF_analysis_gsea_clusterProfiler.R \
        --cnmf_results {input.cnmf_results} \
        --sampleName {wildcards.sample} \
        --outdir {params.analysisdir} \
        --K {wildcards.k} \
        --density.thr {params.threshold} \
        --ranking.type {wildcards.ranking_type} \
        --GSEA.type {wildcards.GSEA_type} \
        --organism {params.organism} \
        --topGenes2look {params.topGenes2look} \
        """

rule clusterProfiler_GSEA_plot:
    input:
        clusterProfiler_result = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis", "GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.txt")
    output:
        clusterProfiler_plot = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis", "top10EnrichedPathways_GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.pdf")
    params:
        outdir = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis"),
        threshold = density,
        scriptdir=config["ScriptDir"],
    resources:
        mem_gb=16,
        runtime_hr=1
    conda:
        config["R_env"],
    shell:
        """
        Rscript {params.scriptdir}/plot_gsea_clusterProfiler.R \
        --clusterProfiler_result {input.clusterProfiler_result} \
        --sampleName {wildcards.sample} \
        --outdir {params.outdir} \
        --K {wildcards.k} \
        --density.thr {params.threshold} \
        --ranking.type {wildcards.ranking_type} \
        --GSEA.type {wildcards.GSEA_type}
        """
        
rule variance_explained:
    input:
        tpm_mtx=os.path.join(config["OutDir"], "{sample}", "{folder}", "{sample}_allK.tpm.h5ad"),
        spectra_consensus=os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "{sample}_allK.spectra.k_{k}.dt_" + density + ".consensus.df.npz"),
        usage_consensus=os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "{sample}_allK.usages.k_{k}.dt_" + density + ".consensus.txt"),
    output:
        Var_k_txt = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis", "metrics.varianceExplained.df.txt")
        Var_k_summary_txt = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis", "summary.varianceExplained.df.txt")
    params:
        cnmf_dir = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density),
        X_normalized_path = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
        outdir = os.path.join(config["OutDir"], "{sample}", "{folder}", "K{k}", "consensus"+density, "lab_analysis"),
        threshold = density,
        scriptdir=config["ScriptDir"],
    resources:
        mem_gb=16,
        runtime_hr=1
    conda:
        config["cnmf_env"],
    shell:
        """
        python {params.scriptdir}/variance_explained_v2.py \
            --consensus_spectra {input.spectra_consensus} \
            --consensus_usage {input.usage_consensus} \
            --topic_sampleName {wildcards.sample} \
            --X_normalized {params.X_normalized_path} \
            --outdir {params.outdir} \
            --k {wildcards.k} \
            --density_threshold {params.threshold} 
    	"""
    
rule aggregate_over_K:
    input:
        cNMF_Results = expand(os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "K{k}", "consensus"+density, "lab_analysis", "cNMF_results.k_{k}.dt_"+density+".RData"),k=cnmf_k),
        clusterProfiler_result = expand(os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "K{k}", "consensus"+density, "lab_analysis", "GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.txt"), k=cnmf_k,ranking_type = ["zscore"]GSEA_type = ["GOEnrichment", "ByWeightGSEA", "GSEA"]),
        Var_k_txt = expand(os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "K{k}", "consensus"+density, "lab_analysis", "metrics.varianceExplained.df.txt"),k=cnmf_k),
        Var_k_summary_txt = expand(os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "K{k}", "consensus"+density, "lab_analysis", "summary.varianceExplained.df.txt"),k=cnmf_k),
    output:
        aggregated_output = os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "lab_analysis", "aggregated.outputs.findK.RData")
    params:
        prev_analysis_dir=os.path.join(config["OutDir"], "{{sample}}", "{{folder}}")
        klist = ",".join([str(i) for i in cnmf_k]),
        scriptdir=config["ScriptDir"],
        outdir=os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "lab_analysis")
    resources:
        mem_gb=16,
        runtime_hr=1
    conda:
        config["R_env"],
    shell:
        """
        Rscript {params.scriptdir}/aggregate_across_K.R \
        --cnmf_results_dir {params.cnmf_results_dir} \
        --clusterProfiler_dir {params.clusterProfiler_dir} \
        --variance_dir {params.clusterProfiler_dir} \
        --var_dir {params.var_dir} \
        --density.thr {params.threshold} \
        --outdir {params.outdir} \
        --sampleName {wildcards.sample} \
        --K.list {params.klist} 
        """

rule findK_plot:
    input:
        toplot = os.path.join(config["analysisDir"], "{folder}/{sample}/acrossK/aggregated.outputs.findK.RData")
    output:
        GSEA_plots = os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "lab_analysis","FindKPlots", "cluster.topic.zscore.by.Pearson.corr.pdf"),
        variance_explained_plot = os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "lab_analysis","FindKPlots", "variance.explained.by.model.pdf"),
    params:
        outdir=os.path.join(config["OutDir"], "{{sample}}", "{{folder}}", "lab_analysis","FindKPlots")
    resources:
        mem_gb=16,
        runtime_hr=1
    conda:
        config["R_env"],
    shell:
        """
        Rscript workflow/scripts/cNMF_findK_plots.R \
        --outdir {params.outdir \
        --sampleName {wildcards.sample} \
        --aggregated.data {input.toplot}
        """