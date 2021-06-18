# format

This is used to convert the count and psi file to csv file.

The row names of csv is ATS sites. The column names of csv is cell barcode.

```bash
go run main.go -h
Usage: main [global options] 

Global options:
        -i, --input  The path to count or psi mtx
        -o, --output The path to output file
        -t, --thread How many threads to use (default: 4)
        -h, --help   Show this help
```
