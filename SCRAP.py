import sys
import getopt
import subprocess
import os
import concurrent.futures


def fastqc(report_folder, sample):
    print('Running FastQC - Sample: ', sample)
    files = " ".join(os.listdir(sample))
    subprocess.run(['fastqc', '-o', report_folder, files, '-t', '2'])


def multiqc():
    print('Running MultiQC')


def readAdapterFile(samples):
    print('Extract adapter and barcode sequences from adapter file')
    for sample in samples:
        print(sample)
        # Read lines in adapter file
        # write to summary.txt
        # List all of the files in the sample folder (1 file for single-end sequencing, 2 files for paired-end sequencing)
        # Print the number of reads in each of the files in the sample directory in the 'sample.summary.txt' file


def flash():
    print('Running FLASH')


def main(argv):
    directory = ''
    adapter_file = ''
    paired_end = ''
    pre_filtered = ''
    reference_directory = ''
    miRBase_species_abbreviation = ''
    genome_species_abbreviation = ''

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
                print(default)
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

    fastqc_report_folder = os.path.join(directory + 'FastQC_Reports')
    multiqc_report_folder = os.path.join(directory + 'MultiQC_Report')

    # fastqc

    try:
        os.mkdir(fastqc_report_folder)
    except:
        print('FastQC_Reports folder already exists.')

    samples = [f.path for f in os.scandir(directory) if f.is_dir()]

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(
            fastqc, fastqc_report_folder, sample) for sample in samples]
        concurrent.futures.wait(futures)

    # multiqc

    try:
        os.mkdir(multiqc_report_folder)
    except:
        print('MultiQC_Report folder already exists.')

    subprocess.run(['multiqc', fastqc_report_folder,
                   '-o', multiqc_report_folder])

    # readAdapterFile()

    # flash()


if __name__ == "__main__":
    main(sys.argv[1:])
