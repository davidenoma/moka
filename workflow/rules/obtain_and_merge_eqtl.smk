from snakemake.io import expand

# # Rule: Obtain eqtl weights from bim file
rule obtain_eqtl_weights:
    input:
        bim = expand("{genotype_file_path}{genotype_prefix}.bim",genotype_file_path=config["genotype_file_path"],
            genotype_prefix=config["genotype_prefix"]),
        # bim = config["genotype_file_path"] + config["genotype_prefix"] + ".bim",
        weights = config["weights_file"]
    params:
        weights_source = config["weights_type"]
    dependencies:
        ""
    shell:
        """
        mkdir -p output_weights
        bim="{input.bim}"
        weights="{input.weights}"
        weights_source={params.weights_source}
        for chr in {{1..22}}; do
            echo "Processing chromosome $chr"
            weights_file="output_weights/{params.weights_source}_weight_file_$chr.csv"
            python scripts/obtain_weights_from_bim.py $bim $chr $weights_file $weights
        done
        """
# #Rule 2: Merge weight files into one
# #
# #
rule merge_eqtl_weights:
    input:
        eqtl_weights=expand("output_weights/{weight_source}_weight_file_{i}.csv",i=range(1,23),weight_source=config["weights_type"]),
        # merged_eqtl_weights=expand("{file_path}{weight_source}_merged_weights.csv",file_path=config["genotype_file_path"],weight_source=config["weights_type"])
    params:
        merged_eqtl_weights = expand("{file_path}{weight_source}_merged_weights.csv",file_path=config["genotype_file_path"],weight_source=config["weights_type"])

    shell:
        """
        # Check if the merged file already exists
        if [ ! -f {params.merged_eqtl_weights} ]; then
            echo 'SNP_ID,Chromosome,Position,Score' > {params.merged_eqtl_weights}
        fi
        # Concatenate contents of all input files to the merged file, skipping header from subsequent files
        for file in {input.eqtl_weights}; do
            tail -n +2 "$file" >> {params.merged_eqtl_weights}
        done
        """