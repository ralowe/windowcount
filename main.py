"""
Get counts for windows across multiple BAM files
"""
import time
from sys import argv
from sys import stderr
import os
import pysam

class BamFile:
    """
    Class to process BAM file
    """
    def __init__(self, filename):
        self.itr = pysam.AlignmentFile(filename, "rb").fetch()
        self.eof = False
        self.read = None
        self.read_next()

    def current_read(self):
        """
        return the current read
        """
        return self.read

    def read_next(self):
        """
        Get the next read
        """
        while True and not self.eof:
            try:
                self.read = next(self.itr)
                if self.read.next_reference_start > self.read.reference_start:
                    self.frag = [self.read.reference_name, self.read.reference_start,
                                 self.read.reference_start + self.read.template_length,
                                 self.read.query_name]
                    break
            except StopIteration:
                self.eof = True

    def get_fragment(self):
        """
        Return the current fragment.
        Fragment is a list of:
        Reference name
        Reference start
        Fragment end
        Query Name
        """
        return self.frag


def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

def print_help():
    """
    Return command options
    """
    print("Usage: windowcount -g <chromosome sizes> -d <directory of BAMs>")
    print("Options:")
    print("\t-g\t\tSet the genome chromosome sizes")
    print("\t-d\t\tDirectory to search for BAM files")
    print("\t-w\t\tSize of window [Int: Default 100]")
    print("\t-s\t\tAmount of shift between each window [Int: Default 50]")
    exit()

def parse_chrom_sizes(filename):
    """
    Get chromosome sizes from fai file
    """
    chromosome_names = []
    chromsome_lengths = []
    for k in open(filename):
        chromosome_name = k.split("\t")[0]
        chromsome_length = int(k.split("\t")[1])
        chromosome_names.append(chromosome_name)
        chromsome_lengths.append(chromsome_length)
    return(chromosome_names, chromsome_lengths)

def extract_read_counts(chrom, spos, epos, sample, buffer):
    """
    For a window defined by chromosome, start position and
    end position obtain the counts from buffer.
    """
    #Check we havent gone to the next chromosome
    while (sample.frag[1] < epos and chrom == sample.frag[0]):
        buffer.append(sample.frag)
        sample.read_next()
        if sample.eof:
            break

    for k in range(len(buffer)-1, -1, -1):
        if buffer[k][2] < spos:
            del buffer[k]

    count = 0
    if buffer:
        count += len(buffer)
    return(buffer, count)

def window_count(chrom, chromsome_length, window_size, shift, files, outfile):
    """
    Method to produce counts in windows across a chromosome for multiple files.
    First set up the read buffers for each file.
    Break chromsome into multiple windows dependent on size and shift (sep).
    """
    position = 0
    end_position = position + window_size
    readbuffers = []
    for _ in range(len(files)):
        readbuffers.append([])
    #Only allow complete windows starting from 0
    while end_position < chromsome_length:
        counts = []
        for index, value in enumerate(files):
            readbuffers[index], file_count = extract_read_counts(chrom, position,
                                                                 end_position, value,
                                                                 readbuffers[index])
            counts.append(file_count)
        if sum(counts) > 0:
            string = f"{chrom}\t{position}\t{end_position}"
            for count in counts:
                string += f"\t{count}"
            outfile.write(f"{string}\n")
        position += shift
        end_position = position + window_size

def check_chrom(chromosome_names, current, files):
    """
    Skip to the right chromosome
    """
    for file in files:
        while chromosome_names.index(file.frag[0]) < current:
            file.read_next()
            if file.eof:
                break


def get_files(direc, out_file):
    """
    Get list of files from directory and write
    header of output file
    """
    files = []
    header = "chrom\tspos\tepos"
    for i in os.listdir(direc):
        if i.find(".bam") != -1 and i.find(".bai") == -1:
            files.append(BamFile(f"{direc}/{i}"))
            header += f"\t{i.split('.bam')[0]}"
    out_file.write(f"{header}\n")
    return files

def main():
    """
    Main function to perform the window count
    """
    myargs = getopts(argv)
    window_size = 100
    shift = 50
    output_filename = "counts.txt"
    if '-g' in myargs:
        ref = myargs["-g"]
    else:
        print_help()
    if '-d' in myargs:
        direc = myargs["-d"]
    else:
        print_help()
    if '-w' in myargs:
        window_size = int(myargs["-w"])
    if '-s' in myargs:
        shift = int(myargs["-s"])
    if '-o' in myargs:
        output_filename = myargs["-o"]
    chromosome_names, chromosome_lengths = parse_chrom_sizes(ref)
    out_file = open(output_filename, "w")
    files = get_files(direc, out_file)
    for index, chromosome_name in enumerate(chromosome_names):
        current_time = time.time()
        check_chrom(chromosome_names, index, files)
        window_count(chromosome_name, chromosome_lengths[index],
                     window_size, shift, files, out_file)
        stderr.write(f"Done {chromosome_name} in {time.time()-current_time}s\n")

if __name__ == "__main__":
    main()
