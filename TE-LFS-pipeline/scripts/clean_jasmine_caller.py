import sys
import math

# change start, end, len to average
# update id to new coordinates and length
# update SVtype because was INS when insurveyor involved

def process_vcf(vcf_file, sample):
    with open(vcf_file, 'r') as infile:
        for line in infile:
            # Keep header lines intact
            if line.startswith('#'):
                print(line.strip())
                continue

            # Split the line into fields
            fields = line.strip().split('\t')
            chrom = fields[0]
            info = fields[7]

            # Parse the INFO field into a dictionary, split only on the first '='
            info_dict = {}
            for kv in info.split(';'):
                if '=' in kv:
                    key, value = kv.split('=', 1)  # Split on the first '=' only
                    info_dict[key] = value
                else:
                    # Handle flag (e.g., PRECISE) as key with a placeholder or True value
                    info_dict[kv] = True

            # Replace POS with rounded AVG_START using math.ceil()
            avg_start = math.ceil(float(info_dict['AVG_START']))
            fields[1] = str(avg_start)

            # Replace END with rounded AVG_END using math.ceil()
            avg_end = math.ceil(float(info_dict['AVG_END']))
            info_dict['END'] = str(avg_end)

            # Replace SVLEN with rounded AVG_LEN using math.ceil()
            avg_len = math.ceil(float(info_dict['AVG_LEN']))
            info_dict['SVLEN'] = str(avg_len)

            # Remove SVINSLEN from the INFO field if present
            if 'SVINSLEN' in info_dict:
                del info_dict['SVINSLEN']

            # Search for ALU, LINE, or SVA in IDLIST and set SVTYPE and ALT accordingly
            idlist = info_dict['IDLIST']
            if 'ALU' in idlist:
                svtype_found = 'ALU'
            elif 'LINE' in idlist:
                svtype_found = 'LINE'
            elif 'SVA' in idlist:
                svtype_found = 'SVA'
            else:
                svtype_found = fields[4]  # Keep original ALT if none found

            # Set SVTYPE and ALT to the found value
            info_dict['SVTYPE'] = svtype_found
            fields[4] = f"<{svtype_found}>" # Modify the ALT column

            # Update the ID field with the new format
            new_id = f"{sample}_{chrom}_{avg_start}_{avg_len}_{svtype_found}"
            fields[2] = new_id  # Modify the ID field

            # Rebuild the INFO field
            info_rebuilt = ';'.join([f"{k}={v}" if v is not True else k for k, v in info_dict.items()])
            fields[7] = info_rebuilt

            # Print the modified VCF line
            print('\t'.join(fields))

# Run the function with the provided VCF file and sample name
process_vcf(sys.argv[1], sys.argv[2])
