
import os
import subprocess

from pathlib import Path
#Path to samples and reference.fasta

path_to_sample = '/data/students_projects/volkova/NIST-hg001-7001'
#path_to_sample = '/data/students_projects/volkova/35_L001/fastq_35_L001'
#path_to_sample = '/data/students_projects/volkova/35_L002/fastq_35_L002'
#path_to_sample = '/data/students_projects/volkova/86_L001/fastq_86_L001'
#path_to_sample = '/data/students_projects/volkova/86_L002/fastq_86_L002'

fq_1 = "/data/students_projects/volkova/NIST-hg001-7001/output_NIST-hg001-7001_R1_trimmed_2.fastq.gz"
fq_2 = "/data/students_projects/volkova/NIST-hg001-7001/output_NIST-hg001-7001_R2_trimmed_2.fastq.gz"


#fq_1 = "/data/students_projects/volkova/35_L001/fastq_35_L001/35_L001_R1.fastq.gz"
#fq_2 = "/data/students_projects/volkova/35_L001/fastq_35_L001/35_L001_R2.fastq.gz"

#fq_1 = "/data/students_projects/volkova/35_L002/fastq_35_L002/35_L002_R1.fastq.gz"
#fq_2 = "/data/students_projects/volkova/35_L002/fastq_35_L002/35_L002_R2.fastq.gz"

#fq_1 = "/data/students_projects/volkova/86_L001/fastq_86_L001/86_L001_R1.fastq.gz"
#fq_2 = "/data/students_projects/volkova/86_L001/fastq_86_L001/86_L001_R2.fastq.gz"

#fq_1 = "/data/students_projects/volkova/86_L002/fastq_86_L002/86_L002_R1.fastq.gz"
#fq_2 = "/data/students_projects/volkova/86_L002/fastq_86_L002/86_L002_R2.fastq.gz"

genome_path = '/data/students_projects/volkova/GRCh38_1/GRCh38_1.fasta'
threads = 16

# 3. Alignment fastq files

import shutil
 
bwa_path = "/data/programs/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
output_sam = "/data/students_projects/volkova/NIST-hg001-7001/NIST-hg001-7001_alignment_output.sam"
#output_sam = "/data/students_projects/volkova/35_L001/35_L001_alignment_output.sam"
#output_sam = "/data/students_projects/volkova/35_L002/35_L002_alignment_output.sam"
#output_sam = "/data/students_projects/volkova/86_L001/86_L001_alignment_output.sam"
#output_sam = "/data/students_projects/volkova/86_L002/86_L002_alignment_output.sam"

folder_name= "NIST-hg001-7001"
#folder_name= "35_L001"
#folder_name= "35_L002"
#folder_name= "86_L001"
#folder_name= "86_L002"


import os
import shutil

# Пути к файлам
WORKDIR = "interval_conversion"
#BED_FILE = "nexterarapidcapture_expandedexome_targetedregions.hg38.bed"
BED_FILE = "TruSeq_exome_targeted_regions.hg38.bed"
REFERENCE_FASTA =  genome_path
#OUTPUT_INTERVAL_LIST = os.path.join(WORKDIR, "nexterarapidcapture_expandedexome_targetedregions.hg38.interval_list")
OUTPUT_INTERVAL_LIST = os.path.join(WORKDIR, "TruSeq_exome_targeted_regions.hg38.interval_list")

# 1. Создаем рабочую папку
os.makedirs(WORKDIR, exist_ok=True)

# 2. Копируем файлы в рабочую папку
shutil.copy(BED_FILE, WORKDIR)
shutil.copy(REFERENCE_FASTA, WORKDIR)

# Обновляем пути к файлам
#BED_FILE = os.path.join(WORKDIR, "nexterarapidcapture_expandedexome_targetedregions.hg38.bed")
BED_FILE = os.path.join(WORKDIR, "TruSeq_exome_targeted_regions.hg38.bed")
REFERENCE_FASTA = os.path.join(WORKDIR, genome_path)

# Функция для получения длины хромосом из FASTA
def get_chromosome_lengths(fasta_file):
    chrom_lengths = {}
    with open(fasta_file, "r") as f:
        chrom_name = None
        length = 0
        for line in f:
            if line.startswith(">"):
                if chrom_name:
                    chrom_lengths[chrom_name] = length
                chrom_name = line.strip().split()[0][1:]  # Убираем '>'
                length = 0
            else:
                length += len(line.strip())
        if chrom_name:
            chrom_lengths[chrom_name] = length
    return chrom_lengths

# Получаем длины хромосом
chrom_sizes = get_chromosome_lengths(REFERENCE_FASTA)

# Читаем BED и записываем в Interval List
with open(BED_FILE, "r") as bed, open(OUTPUT_INTERVAL_LIST, "w") as interval_list:
    # Заголовки
    interval_list.write("@HD\tVN:1.6\tSO:coordinate\n")
    for chrom, length in chrom_sizes.items():
        interval_list.write(f"@SQ\tSN:{chrom}\tLN:{length}\n")

    # Запись интервалов из BED
    for line in bed:
        parts = line.strip().split()
        if len(parts) < 3:
            continue  # Пропускаем строки с ошибками
        chrom, start, end = parts[:3]
        interval_list.write(f"{chrom}\t{start}\t{end}\t+\tGRCh38\n")

print(f"Файл {OUTPUT_INTERVAL_LIST} успешно создан в папке {WORKDIR}!")
