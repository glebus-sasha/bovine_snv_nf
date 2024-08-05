import sys
import pandas as pd

def process_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Initialize lists to store columns
    chromosomes = []
    start_positions = []
    end_positions = []

    # Process lines after the first [Probes] section
    in_probes_section = False

    for line in lines:
        if "[Probes]" in line:
            in_probes_section = True
            continue  # Skip the header line after [Probes]

        if in_probes_section:
            columns = line.split('\t')
            if len(columns) > 3:  # Ensure there are enough columns
                chromosome = columns[5]
                snp_id = columns[2]
                
                # Extract the SNP positions
                snp_info = snp_id.split('_')[-1].strip('()')
                start_position = int(snp_info) - 1  # Convert to 0-based
                end_position = int(snp_info)  # End position is start + 1
                
                chromosomes.append(chromosome)
                start_positions.append(start_position)
                end_positions.append(end_position)

    # Create DataFrame
    bed_df = pd.DataFrame({
        'chromosome': chromosomes,
        'start': start_positions,
        'end': end_positions
    })

    # Save DataFrame to BED file
    bed_file_path = output_file
    bed_df.to_csv(bed_file_path, sep='\t', header=False, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python manifest2bed.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_file(input_file, output_file)
    print(f"Converted {input_file} to {output_file} successfully.")

