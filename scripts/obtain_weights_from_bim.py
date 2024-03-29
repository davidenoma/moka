import subprocess
import argparse




def read_gene_score_file(score_file_path):
    # Create a dictionary to store SNP_ID and Chromosome as keys and Score as values
    score_dict = {}
    with open(score_file_path, 'r') as score_file:
        for line in score_file:
            if line.startswith("Gene_ID"):
                continue  # Skip the header line
            parts = line.strip().split(',')
            if len(parts) != 5:
                # print("Skipping line:", line)
                continue  # Skip lines with unexpected number of fields
            gene_id, snp_id, chromosome,position, score = parts
            chromosome = chromosome.replace("chr", "")
            key = (snp_id, chromosome)
            score_dict[key] = float(score)
    return score_dict

def search_bim_file(bim_file_path, chromosome, output_file_path, score_dict):
    # Open the output file in write mode
    with open(output_file_path, "w") as out_file:
        # Write the header to the output file
        out_file.write("SNP_ID,Chromosome,Position,Score\n")

        # Use subprocess to run the grep command to search the BIM file
        grep_command = f"grep -w -F '{chromosome}' '{bim_file_path}'"
        result = subprocess.run(grep_command, shell=True, stdout=subprocess.PIPE, text=True)

        # Process the matching lines from grep and write them to the output file
        for line in result.stdout.strip().split('\n'):
            parts = line.split()
            snp_id = parts[1]
            snp_chromosome = parts[0]
            snp_position = parts[3]
            key = (snp_id, snp_chromosome)
            if key in score_dict:
                score = score_dict[key]
            else:
                score = 0.0  # Set the score to 0.0 if there's no match in the score file
            output_line = f"{snp_id},{snp_chromosome},{snp_position},{score}\n"
            out_file.write(output_line)

def main():
    parser = argparse.ArgumentParser(description="Search for SNPs in a BIM file by chromosome and create an output file.")

    parser.add_argument("bim_file", help="Path to the BIM file.")
    parser.add_argument("chromosome", help="Chromosome to search for (e.g., 10, 20, etc.).")
    parser.add_argument("output_file", help="Path to the output CSV file for that chromosome.")
    parser.add_argument("score_file", help="Path to the weights file that has the weights for each SNV position")
    args = parser.parse_args()

    bim_file_path = args.bim_file
    chromosome = args.chromosome
    output_file_path = args.output_file
    score_file_path = args.score_file

    # Read the score file into a dictionary
    score_dict = read_gene_score_file(score_file_path)

    # Call the search function to search for SNPs in the BIM file and create the output file
    search_bim_file(bim_file_path, chromosome, output_file_path, score_dict)

    print("Done processing:", output_file_path)

if __name__ == "__main__":
    main()
