# scATS

ATS idenfication and quantification tools

>This code still under development

## Installation

```bash
git clone git@github.com:ygidtu/scATS.git

cd scATS

# install scATS as command line tool
python3 setup.py install

# Or just run source code
python3 main.py
```

## Usage

```bash
Usage: main.py [OPTIONS] COMMAND [ARGS]...

  Welcome

  This function is used to test the function of sashimi plotting

  Created by ygidtu@gmail.com at 2018.12.19 :return:

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  ats      Inference :param debug: enable debug mode
  coexp    Co-expression Note: The input count file must...
  count    Postprocess: count ats and calculate psi
  isoform  Infer isoforms :param debug: enable debug...
```

### ats

Identify ATSs based on aligned BAM files.

```bash
➜  afe git:(master) ✗ python main.py ats --help
Usage: main.py ats [OPTIONS] BAMS...

  Inference

  :param debug: enable debug mode

Options:
  -g, --gtf PATH                 The path to gtf file.   [required]
  -o, --output PATH              The path to output file.   [required]
  -u, --utr-length INTEGER       The length of UTR.
  --n-max-ats INTEGER RANGE      The maximum number of ATSs in same UTR.
  --n-min-ats INTEGER RANGE      The minimum number of ATSs in same UTR.
  --min-ws FLOAT                 The minimum weight of ATSs.
  --min-reads INTEGER            The minimum number of reads in UTR.
  --max-unif-ws FLOAT            The maximum weight of uniform component.
  --max-beta INTEGER             The maximum std for ATSs.
  --fixed-inference              Inference with fixed parameters.
  -d, --debug                    Enable debug mode to get more debugging
                                 information, Never used this while
                                 running.
  -p, --processes INTEGER RANGE  How many cpu to use.
  --remove-duplicate-umi         Only kept reads with different UMIs for
                                 ATS inference.
  --strict                       Only kept reads with different UMIs for
                                 ATS inference.
  -h, --help                     Show this message and exit.
```

### isoform

Isoform assignment.

```bash
➜  afe git:(master) ✗ python main.py isoform --help
Usage: main.py isoform [OPTIONS] BAMS...

  Infer isoforms

Options:
  -i, --ats PATH                 The path to utr file, bed format.
                                 [required]
  -g, --gtf PATH                 The path to reference gtf file.   [required]
  -o, --output PATH              The path to output file.   [required]
  --mu-f INTEGER                 The mean of fragment length.
  --sigma-f INTEGER              The standard deviation of fragment length.
  --min-frag-length INTEGER      The minimum fragment length.
  -d, --debug                    Enable debug mode to get more debugging
                                 information.
  -p, --processes INTEGER RANGE  How many cpu to use.   [1<=x<=48]
  -h, --help                     Show this message and exit.
```

### count

Quantification of ATSs.

```bash
➜  afe git:(master) ✗ python main.py count --help
Usage: main.py count [OPTIONS]

  Postprocess: count ats and calculate psi

Options:
  -i, --ats PATH                 The path to inferred ats sites.
  -b, --bam PATH                 The file contains path to bams.
  --delimiter TEXT               The delimiter of input bam list
  -o, --output PATH              The path to output utr file, bed format.
  -p, --processes INTEGER RANGE  How many cpu to use.
  -c, --compress                 Wheter to save in gzip format.
  --bulk                         Wheter the input bam is Nanopore or
                                 PacBio.
  -h, --help                     Show this message and exit.
```

### Others

The go source code under others is used to format and filter the quantification and PSI matrix.
