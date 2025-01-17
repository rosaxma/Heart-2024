with open("/oak/stanford/groups/engreitz/Users/rosaxma/230828_heart_map/atac_alignment/pool_0/FT005/processing_QC/barcode_matching_full.tsv") as f:
	first_line = f.readline().strip('\n')
print(first_line)
total_reads=first_line.split(" ")[0].split("/")[1]
print(total_reads)
