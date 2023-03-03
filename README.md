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
- [read_fastqcdir](./R/read_fastqcdir.R) is a funciton to read all fastqc report files (.zip files) and aggregate them into one object.
- [read_gff](./R/read_gff.R) is .parser of gff file into R list. attribute columns automatically divided into columns. so the make it easy to understand.
- [read_phylip_dist](./R/read_phylip_dist.R)` reads phylip produced distance matrix file into R environment.

### Getting urls from accession numbers

- [bioLinkGen](./R/bioLinkGen.R)

  ```r 
  ## bioLinkGen Usage
  broweURL(bioLinkGen("SRR17041298"))     ## open SRR17041298 run
  bioLinkGen("SRR17041298",format="html") ## make link for SRR17041298
  ```

<!--
### plink/gwas

plink output files were handled by followings[^1].

#### FILE I/O


- read_plink_ped()      - read plink version 1 file.
- read_plink_genome()   - read `--genome` output of plink
- read_plink_miss()     - read `--miss` output of plink 
- read_plink_pca()      - 
- read_plink_mds()      -
-->





### Visualization




## TODO
- [trimmomatic](./trim.R) call system command to do adapter clipping/low quality removal.
- [NCBI_search](./R/ncbi_search.R) search things in various 



## R pipelines

# NG

https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR11888826/SRR11888826.1

## References

[^1] [plink v1.9](https://www.cog-genomics.org/plink2/)

