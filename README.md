# BKbiokit

`BKbiokit` is a collection of self-made function to process, analysis, and read biological data.

## Installation

```r
  # Install the the development version from GitHub:
  # install.packages("devtools") # if you do not install devtools please install it 
  devtools::install_github("airbj31/BKbiokit")
```
## Functions developed

* for usage information, plz see `??<command>`
* entiries without link are planning.

### FILE I/O

- [read_fastqc](./R/read_fastqc.R) is a function to load fastqc report file (.zip file) into R environment.
- [read_fastqcdir](./R/read_fastqcdir.R) is a function to read all fastqc report files (.zip files) and aggregate them into one object.
- [read_gff](./R/read_gff.R) is .parser of gff file into R list. attribute columns automatically divided into columns. so the make it easy to understand.
- [read_phylip_dist](./R/read_phylip_dist.R)` reads phylip produced distance matrix file into R environment.

### DATA manipulation.

- [BWTcode and vBWTcode](./R/BWT.R) are function inspired by Burrow-Wheeler Transform and one-hot encoding. The function is developed to analyze short tandem repeats with K-mer approach. The main purpose of the program is merging repeat elements into one K-mer. e.g. `TTAGG`, `GTTAG`, `GGTTA`, `AGGTT`, `TAGGT` are merged into `AGGTT`.**plz note that the function choose 1st entry from sorting out all the possible k-mers lexicographical order.** When Reverse Complement options is on, we also consider the reverse complement using Biostrings::reverseComplement function. `vBWTcode` is vectorized function to use function in tidyverse way.

```r
BWTcode("TTAGGG",rc=FALSE)
BWTcode("TTAGGG",rc=TRUE)

```


### Getting urls from accession numbers

- [bioLinkGen](./R/bioLinkGen.R)

  ```r 
  ## bioLinkGen Usage
  broweURL(bioLinkGen("SRR17041298"))     ## open SRR17041298 run
  bioLinkGen("SRR17041298",format="html") ## make link for SRR17041298
  ```

<!--


- read_plink_ped()      - read plink version 1 file.
- read_plink_genome()   - read `--genome` output of plink
- read_plink_miss()     - read `--miss` output of plink 
- read_plink_pca()      - 
- read_plink_mds()      -

## TODO
- [trimmomatic](./trim.R) call system command to do adapter clipping/low quality removal.
- [NCBI_search](./R/ncbi_search.R) search things in various 



## R pipelines

# NG

https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR11888826/SRR11888826.1

## References

[^1] [plink v1.9](https://www.cog-genomics.org/plink2/)

-->
