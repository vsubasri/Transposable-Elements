import sys

# Make new ID column
# Clean ALT column
# Remove lines containing "INCOMPLETE_ASSEMBLY" because length can't be estimated and need that for merging

def process_vcf(vcf_file, sample):
    with open(vcf_file, 'r') as infile:
        for line in infile:
            # Keep header lines intact
            if line.startswith('#'):
                print(line.strip())
                continue

            # Skip lines containing "INCOMPLETE_ASSEMBLY"
            if "INCOMPLETE_ASSEMBLY" in line:
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            alt = "INS" 
            svlen_info = next((x.split('=')[1] for x in fields[7].split(';') if x.startswith('SVLEN')), None)
            new_id = f"{sample}_{chrom}_{pos}_{svlen_info}_{alt}"  # Use INS in the ID as well

            # Modify ID fields
            fields[2] = new_id

            # Print the modified VCF line
            print('\t'.join(fields))

# Simply pass the arguments directly to the function
process_vcf(sys.argv[1], sys.argv[2])

