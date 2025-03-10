import re

# Input and output VCF file paths
input_vcf = "/home/chbope/extension/data/annotations/T1544_annotated_variants.vcf"
output_vcf = "/home/chbope/extension/data/annotations/T1544_annotated_variants_clean.vcf"


# Open the VCF files
with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
    for line in infile:
        if line.startswith("#"):  # Keep header lines unchanged
            outfile.write(line)
        else:
            fields = line.strip().split("\t")
            info_field = fields[7]  # INFO is the 8th field in VCF

            # Safely parse the INFO field
            info_items = info_field.split(";")
            info_dict = {}
            for item in info_items:
                if "=" in item:  # Only process key=value pairs
                    key, value = item.split("=", 1)  # Allow values to contain '='
                    info_dict[key] = value

            # Process the SVLEN field if it exists
            if "SVLEN" in info_dict:
                svlen_values = re.split(",|;", info_dict["SVLEN"])  # Split on commas or semicolons
                svlen_values = [int(v) for v in svlen_values if v.isdigit() or v.lstrip("-").isdigit()]  # Filter valid integers

                # Choose the cleaning strategy
                if svlen_values:
                    cleaned_svlen = svlen_values[0]  # Use the first value
                    # cleaned_svlen = max(svlen_values, key=abs)  # Use the absolute maximum
                    info_dict["SVLEN"] = str(cleaned_svlen)
                else:
                    del info_dict["SVLEN"]  # Remove invalid SVLEN field

            # Reconstruct the INFO field
            new_info_field = ";".join(f"{k}={v}" for k, v in info_dict.items())
            fields[7] = new_info_field
            outfile.write("\t".join(fields) + "\n")

