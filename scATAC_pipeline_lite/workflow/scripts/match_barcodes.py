import numpy as np

import matcha


def read_barcodes(path):
    with open(path, "rt") as file:
        bc = [b.rstrip("\n") for b in file]

    return bc


def match_bc(
    fastqs,
    whitelist,
    max_barcode_dist,
    fastq1_out_path,
    fastq2_out_path,
    qc_path,
    threads,
):
    f = matcha.FastqReader(threads=threads)
    f.add_sequence("R1", fastqs["R1"], output_path=fastq1_out_path)
    f.add_sequence("R2", fastqs["R2"])
    f.add_sequence("R3", fastqs["R3"], output_path=fastq2_out_path)

    barcode_sequences = read_barcodes(whitelist)
    cell_barcode = matcha.HashMatcher(
        sequences=barcode_sequences,
        labels=barcode_sequences,
        max_mismatches=max_barcode_dist,
        subsequence_count=2,
    )
    f.add_barcode("cell", cell_barcode, "R2")
    f.set_output_names("{read_name} CB:Z:{cell}")

    barcode_counts = np.zeros(max_barcode_dist + 2, int)

    total_reads = 0
    total_pass = 0

    chunk_size = 10000
    while f.read_chunk(chunk_size):
        pass_filter = (f.get_match_result("cell", "dist") <= max_barcode_dist) & (
            f.get_match_result("cell", "second_best_dist")
            > f.get_match_result("cell", "dist")
        )

        total_reads += len(pass_filter)
        total_pass += pass_filter.sum()
        values, counts = np.unique(
            f.get_match_result("cell", "dist"), return_counts=True
        )
        barcode_counts[np.minimum(values, max_barcode_dist + 1)] += counts

        f.write_chunk(pass_filter)

    with open(qc_path, "w") as stats_output:
        print(
            f"{total_pass}/{total_reads} reads passing, ({total_pass/total_reads*100:.2f}%)\n",
            file=stats_output,
        )
        print("mismatches\treads", file=stats_output)
        for dist in range(max_barcode_dist + 2):
            print(
                dist if dist <= max_barcode_dist else f">{max_barcode_dist}",
                barcode_counts[dist],
                sep="\t",
                file=stats_output,
            )


try:
    max_barcode_dist = snakemake.params["barcode_dist"]

    fastq1_out_path = snakemake.output["fastq1_bc"]
    fastq2_out_path = snakemake.output["fastq2_bc"]

    qc_path = snakemake.output["qc_matching"]

    threads = snakemake.threads*4

    fastqs = {
        "R1": snakemake.input["fq_R1"],
        "R2": snakemake.input["fq_BC"],
        "R3": snakemake.input["fq_R2"],
    }

    whitelist = snakemake.input["whitelist"]

    match_bc(
        fastqs,
        whitelist,
        max_barcode_dist,
        fastq1_out_path,
        fastq2_out_path,
        qc_path,
        threads,
    )

except NameError:
    pass
