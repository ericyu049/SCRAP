import sys
import getopt
import subprocess
import os
import glob
import concurrent.futures
import pandas as pd
import gzip
from datetime import datetime


# Parameters for filtering pre-miRNA BLAST results

minimum_pre_miRNA_evalue = .05
minimum_pre_miRNA_bit_score = 50

# Parameters for filtering tRNA BLAST results

minimum_tRNA_evalue = .05
minimum_tRNA_bit_score = 50

# Parameters for filtering BLAST results

minimum_evalue = .05
minimum_length = 14
maximum_mismatch = 2
maximum_gap = 1
maximum_gap_if_maximum_mismatch = 0

# Minimum length of sequence after the sncRNA within a read (after adapter and barcode removal)

minimum_length_after_sncRNA = 15

directory = ''
adapter_file = ''
paired_end = ''
pre_filtered = ''
reference_directory = ''
miRBase_species_abbreviation = ''
genome_species_abbreviation = ''

samples = []
samples_with_path = []
fastqc_report_folder = ''
multiqc_report_folder = ''

five_prime_adapter = ''
three_prime_adapter = ''
five_prime_barcode = ''
three_prime_barcode = ''


def fastqc():
    print('Start running FastQC -- ', datetime.now())
    try:
        os.mkdir(fastqc_report_folder)
    except:
        print('FastQC_Reports folder already exists.')

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(
            fastqc_worker, directory, sample) for sample in samples]
        concurrent.futures.wait(futures)

    print('Done running FastQC -- ', datetime.now())



def fastqc_worker(directory, sample):
    print('Running FastQC - Sample: ', sample, ' -- ', datetime.now())
    files = ' '.join(glob.glob(os.path.join(directory, sample) + '/*fastq*'))
    subprocess.run(['fastqc', '-o', fastqc_report_folder, files, '-t', '2'])
    print('Done running FastQC - Sample: ', sample, ' -- ', datetime.now())


def multiqc():
    print('Running MultiQC -- ', datetime.now())
    try:
        os.mkdir(multiqc_report_folder)
    except:
        print('MultiQC_Report folder already exists.')
    # subprocess.run(['multiqc', fastqc_report_folder,
    #                '-o', multiqc_report_folder])

    print('Done running MultiQC -- ', datetime.now())

def sampleWorker(df, directory, sample):
    countReads(directory, sample)
    # if (paired_end == 'yes'):
    #         flash()
    # if five_prime_adapter:
    #     print('five prime adapter')
    # if three_prime_adapter:
    #     print('three prime adapter')
    # if five_prime_barcode:
    #     print('five prime barcode')
    # if three_prime_barcode:
    #     print('three prime barcode')



def readAdapterFile(df, sample):
    print('Extract adapter and barcode sequences from adapter file')

    five_prime_adapter = df["5'Adapter"][sample]
    three_prime_adapter = df["3'Adapter"][sample]
    five_prime_barcode = df["5'Barcode"][sample]
    three_prime_barcode = df["3'Barcode"][sample]

    sample_summary_txt = os.path.join(directory, sample, sample + '.summary.txt')
    if os.path.exists(sample_summary_txt):
        os.remove(sample_summary_txt)

    with open(sample_summary_txt, 'a') as file:
        line1 = f"{sample}\n"
        line2 = f"5'-{five_prime_adapter}{five_prime_barcode}...sncRNA-targetRNA...{three_prime_adapter}{three_prime_barcode}\n"
        line3 = f"miRBase Species Abbreviation: {miRBase_species_abbreviation}\n"
        line4 = f"Genome Species Abbreviation: {genome_species_abbreviation}\n"
        line5 = f"Start: {datetime.now()}\n"

        #write to the file
        file.writelines([line1, line2, line3, line4, line5])


def countReads(directory, sample):
    files = glob.glob(os.path.join(directory, sample) + '/*fastq*')

    # Separate R1 and R2 into 2 subproccess workers for counting reads.

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(
            countReadsSubWorker, zipfile, directory, sample) for zipfile in files]
        concurrent.futures.wait(futures)
        

def countReadsSubWorker(zipfile, directory, sample):
    print('Counting reads for ', zipfile, ' -- ', datetime.now())
    with gzip.open(zipfile, 'rt') as f_in:
        num_reads = sum(1 for line in f_in) // 4

    sample_summary_txt = os.path.join(directory, sample, sample + '.summary.txt')

    with open(sample_summary_txt, 'a') as file:
        file.write(f"{zipfile} raw reads: {num_reads}\n")
    print('Done counting reads for ', zipfile, ' -- ', datetime.now())


def flash():
    print('Running FLASH -- ', datetime.now())


########## Pipeline begins here #############
if __name__ == '__main__':
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "hd:a:p:f:r:m:g:")
    for opt, arg in opts:
        match opt:
            case '-h':
                print('help')
                sys.exit()
            case '-d':
                directory = arg
            case '-a':
                adapter_file = arg
            case '-p':
                paired_end = arg
            case '-f':
                pre_filtered = arg
            case '-r':
                reference_directory = arg
            case '-m':
                miRBase_species_abbreviation = arg
            case '-g':
                genome_species_abbreviation = arg
            case default:
                print('invalid command')
                sys.exit()

    if directory == '':
        print("Error: Path to sample directories not provided [-d]")
        sys.exit(1)
    if adapter_file == '':
        print("Error: Path to adapter file not provided [-a]")
        sys.exit(1)
    if paired_end == '':
        print("Error: Are samples paired-end? [-p]")
        sys.exit(1)
    if pre_filtered == '':
        print("Error: Filter out pre-miRNAs and tRNAs? [-f]")
        sys.exit(1)
    if reference_directory == '':
        print("Error: Path to reference directory not provided [-r]")
        sys.exit(1)
    if miRBase_species_abbreviation == '':
        print("Error: miRBase species abbreviation not provided [-m]")
        sys.exit(1)
    if genome_species_abbreviation == '':
        print("Error: Species genome abbreviation not provided [-g]")
        sys.exit(1)

    df = pd.read_csv(adapter_file, sep='\t', skiprows=[0], index_col=0, names=[
        "5'Adapter", "3'Adapter", "5'Barcode", "3'Barcode"])
    samples = df.index.unique()

    fastqc_report_folder = os.path.join(directory + 'FastQC_Reports')
    multiqc_report_folder = os.path.join(directory + 'MultiQC_Report')

    # fastqc()
    # multiqc()

    for sample in samples:
        readAdapterFile(df, sample)



    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(
            sampleWorker, df, directory, sample) for sample in samples]
        concurrent.futures.wait(futures)
    

    
