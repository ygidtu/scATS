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

- default julia kernel support multi-threading, please set the number of threads to use by `julia -t n` or `export JULIA_NUM_THREADS=n`

- `--using-R` is conflist with Theads, therefore pleas using `--using-R n` to enable  multi processing to process data

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

### 3. merge output peaks

```bash
> julia merge.jl --help
usage: merge.jl -o OUTPUT [-e EXPAND] [-h] ats...

positional arguments:
  ats                  Path to atsmix output file

optional arguments:
  -o, --output OUTPUT  Path to output.
  -e, --expand EXPAND  How many bp to expand (type: Int64, default:
                       100)
  -h, --help           show this help message and exit
```

### 4. quantification

```bash
> julia quant.jl --help
usage: quant.jl -i INPUT -c CELLRANGER -o OUTPUT [-p PROCESS] [-h]

optional arguments:
  -i, --input INPUT     Path to merged peaks bed
  -c, --cellranger CELLRANGER
                        Path to cellranger outs directory
  -o, --output OUTPUT   Prefix of output file
  -p, --process PROCESS
                        How many processes to use (type: Int64,
                        default: 1)
  -h, --help            show this help message and exit
```
