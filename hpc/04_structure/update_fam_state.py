import csv
import os

# Files
metadata_file = 'clone_corrected_metadata.csv'
fam_file = 'karan_snps_filtered_plink.fam'
output_fam = 'karan_snps_filtered_plink_updated.fam'

# 1. Load the state metadata into a dictionary
state_map = {}
# utf-8-sig handles any hidden BOM characters at the start of CSVs
with open(metadata_file, 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = row['Sakalidis_lab_code'].strip()
        state = row['State'].strip()
        state_map[sample] = state

# 2. Rewrite the .fam file replacing the 1st column (Family ID) with State
updated_count = 0
with open(fam_file, 'r') as f_in, open(output_fam, 'w') as f_out:
    for line in f_in:
        parts = line.strip().split()
        if not parts:
            continue
        
        sample_id = parts[1] # Use the Individual ID to find a match
        
        # If the sample is in our metadata, replace the first column
        if sample_id in state_map:
            parts[0] = state_map[sample_id]
            updated_count += 1
            
        f_out.write(" ".join(parts) + "\n")

print(f"Successfully updated {updated_count} samples!")
print(f"New file saved as: {output_fam}")
