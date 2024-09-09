import sys

# Make new ID column 
# Clean ALT column
# Filter out rows where FILTER column is not PASS

def process_vcf(vcf_file, sample):
    with open(vcf_file, 'r') as infile:
        for line in infile:
            # Keep header lines intact
            if line.startswith('#'):
                print(line.strip())
                continue

            # Split the line into fields
            fields = line.strip().split('\t')

            # Only process lines where the FILTER column is "PASS"
            if fields[6] != "PASS":
                continue

            # Extract necessary fields for processing
            chrom = fields[0]
            pos = fields[1]
            svlen = next((x.split('=')[1] for x in fields[7].split(';') if x.startswith('SVLEN')), None)
            alt = fields[4]
            svtype = fields[4].strip('<>').split(':')[-1]

            # Append '_alt' to the new_id
            new_id = f"{sample}_{chrom}_{pos}_{svlen}_{svtype}"

            # Modify ALT and ID fields
            fields[2] = new_id

            # Print the modified line
            print('\t'.join(fields))

# Run the function
process_vcf(sys.argv[1], sys.argv[2])

