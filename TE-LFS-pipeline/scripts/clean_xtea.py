import sys

# Make new ID column 
# Clean ALT column
# Clean SVTYPE info
# Filter out lines containing "orphan"

def process_vcf(vcf_file, sample):
    with open(vcf_file, 'r') as infile:
        for line in infile:
            # Keep header lines intact
            if line.startswith('#'):
                print(line.strip())
                continue

            # Filter out lines containing the word "orphan"
            if "orphan" in line:
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]

            # SVTYPE 
            alt = fields[4]
            svtype = fields[4].strip('<>').split(':')[-1]

            # Create new ID column using sample, chrom, pos, svlen, and alt
            svlen_info = next((x.split('=')[1] for x in fields[7].split(';') if x.startswith('SVLEN')), None)
            new_id = f"{sample}_{chrom}_{pos}_{svlen_info}_{svtype}"

            # Modify ALT, SVTYPE, and ID fields
            fields[2] = new_id  

            # Rebuild INFO with modified SVTYPE
            info_fields = []
            for item in fields[7].split(';'):
                if item.startswith('SVTYPE'):
                    info_fields.append(f"SVTYPE={svtype}")  # Replace SVTYPE 
                else:
                    info_fields.append(item)
            fields[7] = ';'.join(info_fields)

            # Print the modified VCF line
            print('\t'.join(fields))

# Simply pass the arguments directly to the function
process_vcf(sys.argv[1], sys.argv[2])

