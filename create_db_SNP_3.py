import os
import subprocess
import urllib.request
import pandas as pd

def download_file(url, filename):
    print(f"Downloading {filename}...")
    urllib.request.urlretrieve(url, filename)
    print(f"Download of {filename} complete.")

# 1. Download archiew with VCF
vcf_url = "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
vcf_file = "GCF_000001405.40.gz"
if not os.path.exists(vcf_file):
    download_file(vcf_url, vcf_file)

# 2. Unzip VCF
subprocess.run(['gunzip', '-k', vcf_file], check=True)
unzipped_vcf_file = vcf_file.replace(".gz", "")
print(f"Unzipped: {unzipped_vcf_file}")


# 3. Compress file to BGZF
bgzipped_vcf_file = unzipped_vcf_file + ".gz"
subprocess.run(['bgzip', '-c', unzipped_vcf_file], stdout=open(bgzipped_vcf_file, 'wb'), check=True)
print(f"BGZF compressed: {bgzipped_vcf_file}")

# 4. Index VCF
subprocess.run(['bcftools', 'index', '--threads', '16', bgzipped_vcf_file], check=True)
print(f"Indexed: {bgzipped_vcf_file}")


# 5. Download Assembly Report
def download_file(url, filename):
    print(f"Downloading {filename}...")
    urllib.request.urlretrieve(url, filename)
    print(f"Download of {filename} complete.")


assembly_report_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
assembly_report_file = "GCF_000001405.40_GRCh38.p14_assembly_report.txt"
if not os.path.exists(assembly_report_file):
     download_file(assembly_report_url, assembly_report_file)

# 6. Creating a file with matching chromosome numbers

chrnames_file = "GCF_000001405.40_GRCh38.p14_assembly_report.chrnames"

final_chrnames_file = "GCF_000001405.40_GRCh38.p14_assembly_report.filtered.chrnames"

# 7. Extracting two columns from an Assembly Report
def extract_columns(report_file, output_file):
    """Extracting two columns from an Assembly Report."""
    with open(report_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
               continue
            columns = line.strip().split('\t')
            if len(columns) >= 7:
                outfile.write(f"{columns[6]} {columns[-1]}\n")
    print(f"Extracted columns to: {output_file}")

extracted_file = "GCF_000001405.40_GRCh38.p14_assembly_report_extracted.txt"
extract_columns(assembly_report_file, extracted_file)

# 8. Delete lines from NA and saving to a separate file
filtered_file_na = "GCF_000001405.40_GRCh38.p14_assembly_report_filtered_na.txt"

def filter_na(input_file, output_file):
    """Delete lines from NA and saving to a separate file."""
    df = pd.read_csv(input_file, sep=' ', header=None, na_values=['NA', 'na', 'Na', 'NA '])  
    df = df.dropna()  # Delete line with  NA
    df.to_csv(output_file, sep=' ', index=False, header=False)
    print(f"Filtered data (no NA) saved to: {output_file}")

filter_na(extracted_file, filtered_file_na)

# 9. Filtering empty columns
final_filtered_file = "GCF_000001405.40_GRCh38.p14_assembly_report_final_filtered.txt"

def filter_empty_columns(input_file, output_file):
    """Filtering empty columns and saving to a separate file."""
    df = pd.read_csv(input_file, sep=' ', header=None)
    df = df[df[1].notna()]  #Delete liness where the second column is empty
    df.to_csv(output_file, sep=' ', index=False, header=False)
    print(f"Filtered empty columns and saved to: {output_file}")

filter_empty_columns(filtered_file_na, final_filtered_file)

#Annotation vcf
annotation_file = "GCF_000001405.40_GRCh38.p14_assembly_report_final_filtered.txt"
annotated_vcf_file = "GRCh38.dbSNP.vcf.gz"
subprocess.run([
    'bcftools', 'annotate',
    '--rename-chrs', annotation_file,
    '--threads', '16', '-Oz',
    '-o', annotated_vcf_file,
    bgzipped_vcf_file
], check=True)
print(f"Annotated VCF saved as: {annotated_vcf_file}")

