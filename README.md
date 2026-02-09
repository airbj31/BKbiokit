# BKbiokit

`BKbiokit` is a collection of self-made function to process, analysis, and read biological data.

## Installation

```r
  # Install the the development version from GitHub:
  # install.packages("devtools") # if you do not install devtools please install it 
  devtools::install_github("airbj31/BKbiokit")
```

## Help

To learn how to use the package, please see the vignette.

```r
browseVignettes("BKbiokit")
```

## Functions developed

* for usage information, plz see `??<command>`
* please note that entiries without links are planning/currently under development.

### FILE I/O

for more example about the functions,please see the vignette.


- [read_fastqc](./R/read_fastqc.R) is a function to load fastqc report file (.zip file) into R environment.
- [read_fastqcdir](./R/read_fastqcdir.R) is a function to read all fastqc report files (.zip files) and aggregate them into one object.
- [read_gff](./R/read_gff.R) is a gff parser which convert gff file into R list. attribute columns automatically divided into columns. so the function make the file easy to understand. the lists were named as `Type` columns.
- [write_gbff](./R/write_gbff.R) write genbank file format file[^2].
- read_gbff read genbank flat format file[^3].
- [read_phylip_dist](./R/read_phylip_dist.R) read phylip produced distance matrix file into R environment. The loaded distance matrix was automatically converted into `dist` object

### DATA manipulation.

- [BWTcode and vBWTcode](./R/BWT.R) are function inspired by Burrow-Wheeler Transform and one-hot encoding. The function is developed to analyze short tandem repeats with K-mer approach. The main purpose of the program is merging repeat elements into one K-mer. e.g. `TTAGG`, `GTTAG`, `GGTTA`, `AGGTT`, `TAGGT` are merged into `AGGTT`.**plz note that the function choose 1st entry from sorting out all the possible k-mers lexicographical order.** When Reverse Complement options is on, we also consider the reverse complement using Biostrings::reverseComplement function. `vBWTcode` is vectorized function to use function in tidyverse way.  

  ```r
  BWTcode("TTAGGG",rc=FALSE)
  BWTcode("TTAGGG",rc=TRUE)

  ```
- [getKmer](./R/getKmer.R)  - R version of K-mer extraction


  
### NCBI E-utility and/or Related things

NCBI E-utility APIs are implementing into `BKbiokit` for easy NCBI search.

- [einfo()](./R/einfo.R).  
  ```r
  einfo()            ## get db information
  einfo(db="pubmed") ## get information about pubmed db and field/link information.
  einfo(db="nuccore) ## get information about nuccore (nucelotide) db and field/link information
  ```
- [esearch()](./R/einfo.R)
  ```r
  esearch("apis cerana")
  
  ```

- [download_SRA()](./R/download_SRA.R) this function downloads SRA file and decrypt it using fasterq-dump. sratookit and aws client are required to use this function.
  ```r
  download_SRA("SRR17041298")

  ```


### Getting urls from accession numbers

- [bioLinkGen](./R/bioLinkGen.R)

  ```{r cache=T} 
  ## bioLinkGen Usage
  browseURL(bioLinkGen("SRR17041298"))     ## open SRR17041298 run
  bioLinkGen("SRR17041298",format="html") ## make link for SRR17041298
  ```

### GEO

- [searchGEO](./R/GEO.R)
- [download_GEO](./R/GEO.R)
- [read_GSE](./R/GEO.R)



<!--


- read_plink_ped()      - read plink version 1 file[^1].
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

## References 

[^1]: [plink v1.9](https://www.cog-genomics.org/plink2/)

[^2]: [gb release file](https://ftp.ncbi.nih.gov/genbank/gbrel.txt)

[^3]: [gbff file format](https://www.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/file-formats/about-ncbi-gbff/)
