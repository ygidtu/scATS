# scATS.jl

A sutie of scripts for infer alternative transcription start sites based on 5' single cell sequencing

This version should be faster than python version of scATS.

## Installation

```bash
git clone 

cd scATS
git checkout julia

julia -e "using Pkg; Pkg.resolve(); Pkg.instantiate()"
```

## Usage

```bash
❯ julia --project=@. ./main.jl
[2021-12-08 15:33:15 - info | root]: please set running subcommand: ats or count
usage: main.jl cmd
```

### 1. infer ATS

```bash
❯ julia --project=@. ./main.jl ats -h
usage: main.jl -b BAM -g GTF -o OUTPUT [-u UTR-LENGTH] [-t THREADS]
               [--n-max-ats N-MAX-ATS] [--n-min-ats N-MIN-ATS]
               [--min-ws MIN-WS] [--max-unif-ws MAX-UNIF-WS]
               [--max-beta MAX-BETA] [--step-size STEP-SIZE]
               [--nround NROUND] [--min-reads MIN-READS] [--seed SEED]
               [--fixed-inference] [-p PROCESS] [-h]

optional arguments:
  -b, --bam BAM         Path to bam list
  -g, --gtf GTF         Path to reference annotation file, GTF format
  -o, --output OUTPUT   The path to output file
  -u, --utr-length UTR-LENGTH
                        The length of UTR region (type: Int64,
                        default: 500)
  -t, --threads THREADS
                        How many threads to use (type: Int64, default:
                        1)
  --n-max-ats N-MAX-ATS
                        the maximum of ats inside a utr (type: Int64,
                        default: 5)
  --n-min-ats N-MIN-ATS
                        the minimum of ats inside a utr (type: Int64,
                        default: 1)
  --min-ws MIN-WS       min ws (type: Float64, default: 0.1)
  --max-unif-ws MAX-UNIF-WS
                        maximum uniform ws (type: Float64, default:
                        0.1)
  --max-beta MAX-BETA   maximum beta (type: Float64, default: 50.0)
  --step-size STEP-SIZE
                        step size (type: Int64, default: 5)
  --nround NROUND       number of round to test (type: Int64, default:
                        50)
  --min-reads MIN-READS
                        minimum reads to construct ATS (type: Int64,
                        default: 5)
  --seed SEED           seed for consistance results (type: Int64,
                        default: 42)
  --fixed-inference     inference with fixed parameters
  -p, --process PROCESS
                        Number of processes to used. (type: Int64,
                        default: 1)
  -h, --help            show this help message and exit
```

### 2. quantification

```bash
❯ julia --project=@. ./main.jl count -h
usage: main.jl -i INPUT -b BAM -o OUTPUT [-p PROCESS] [-h]

optional arguments:
  -i, --input INPUT     Path to merged peaks bed.
  -b, --bam BAM         Path to bam list.
  -o, --output OUTPUT   Prefix of output file.
  -p, --process PROCESS
                        Number of processes to used. (type: Int64,
                        default: 1)
  -h, --help            show this help message and exit
```

#### Run Example

```bash
cd simulation

python simulation.py /path/to/reference/gtf /path/to/reference/fasta --total 10000
```

Generate bam list file as follow

```bash
./simu.bam    simu
```

Run inference

```bash
julia --project=@.  ../main.jl ats -b ./bam.list  -g ./tss.gtf -o ./simu.txt -p 6
```
