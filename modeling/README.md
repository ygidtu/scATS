# ATAMix

A sutie of scripts for infer alternative transcription start sites based on 5' single cell sequencing

## Installation

```bash
git clone 

cd

julia build.jl
```

> `samtools` maybe required, if there is not `.bai` index for input bam file

## Usage

### 1. prepare UTR region based on gtf file

```bash
> julia ./prepare.jl --help
usage: prepare.jl -g GTF -p PREFIX [-d DISTANCE] [-h]

optional arguments:
  -g, --gtf GTF         Path to gtf file
  -p, --prefix PREFIX   Prefix of output file
  -d, --distance DISTANCE
                        The distance between utr to merge. (type:
                        Int64, default: 0)
  -h, --help            show this help message and exit
```

### 2. infer ATS

```bash
> julia run.jl --help
usage: run.jl -i INPUT -b BAM -o OUTPUT [-d DISTANCE] [--using-R]
              [-p PROCESS] [--n-max-ats N-MAX-ATS]
              [--n-min-ats N-MIN-ATS] [--mu-f MU-F]
              [--sigma-f SIGMA-F] [--min-ws MIN-WS]
              [--max-beta MAX-BETA] [--min-reads MIN-READS]
              [--single-end] [--fixed-inference] [--verbose] [-h]

optional arguments:
  -i, --input INPUT     Path to utr bed
  -b, --bam BAM         Path to bam file
  -o, --output OUTPUT   Prefix of output file
  -d, --distance DISTANCE
                        The minimum distance of read in utr. (type:
                        Int64, default: 1500)
  --using-R             whether to use R version of model
  -p, --process PROCESS
                        How many processes to use (type: Int64,
                        default: 1)
  --n-max-ats N-MAX-ATS
                        the maximum of ats inside a utr (type: Int64,
                        default: 5)
  --n-min-ats N-MIN-ATS
                        the minimum of ats inside a utr (type: Int64,
                        default: 1)
  --mu-f MU-F           mu f (type: Int64, default: 300)
  --sigma-f SIGMA-F     sigma-f (type: Int64, default: 50)
  --min-ws MIN-WS       min ws (type: Float64, default: 0.01)
  --max-beta MAX-BETA   maximum beta (type: Float64, default: 50.0)
  --min-reads MIN-READS
                        minimum reads to construct ATS (type: Int64,
                        default: 0)
  --single-end          whether this is sinle-end sequencing
  --fixed-inference     inference with fixed parameters
  --verbose             whether to display additional message
  -h, --help            show this help message and exit
```
