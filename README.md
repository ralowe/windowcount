# windowcount
Count reads from several BAMs in windows across the genome of interest.

## Install

Requires `pysam` (https://pysam.readthedocs.io/en/latest/)

This can be installed using pip

`pip install pysam`

Then clone repository

`git clone https://github.com/ralowe/windowcount.git`

## Running 

`windowcount` requires 2 inputs. A chromosome sizes file and a directory location of bam files.

To create the chromosome sizes file run samtools (http://www.htslib.org) on .fa file:

`samtools faidx ref.fasta`

Then run `windowcount`

`python3 main.py -g ref.fasta.fai -d directory/of/BAMS`

Full usage can be found by running

`python3 main.py`

```
Usage: windowcount -g <chromosome sizes> -d <directory of BAMs>
Options:
	-g		Set the genome chromosome sizes
	-d		Directory to search for BAM files
	-w		Size of window [Int: Default 100]
	-s		Amount of shift between each window [Int: Default 50]
```


