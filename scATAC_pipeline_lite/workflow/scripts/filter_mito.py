import pysam


def filter_mito(in_path, out_path, qc_path, mito_chr):
    infile = pysam.AlignmentFile(in_path, "rb")
    outfile = pysam.AlignmentFile(out_path, "wb", template=infile)

    num_mito = 0
    num_non_mito = 0
    for a in infile:
        if a.flag & 260 == 0:  # Alignment is mapped and is primary
            if a.reference_name == mito_chr:
                num_mito += 1
            else:
                num_non_mito += 1
                outfile.write(a)

    with open(qc_path, "w") as qc_file:
        print("Non-Mitochondrial\tMitochondrial", file=qc_file)
        print(f"{num_non_mito}\t{num_mito}", file=qc_file)


(in_path,) = snakemake.input
out_path = snakemake.output["bam"]
qc_path = snakemake.output["qc"]
mito_chr = snakemake.params["mitochr"]

filter_mito(in_path, out_path, qc_path, mito_chr)
