#' efetch
#' UNDER DEV
#' search NCBI db using esearch API. this function is under development
#' the function returns
#'      Returns formatted data records for a list of input UIDs
#'      Returns formatted data records for a set of UIDs stored on the Entrez History server
#'
#' @param db (required) Database name. NCBI DBs were described in details.
#' @param id (required) Unique identifier or accession number. NCBI uid. vector object.
#' @param input Read identifier(s) from file (one line).
#' @param api_key (optional) NCBI api_key. you can get your own api_key from NCBI setting page after log-in. To use the function easily, you need to declare "NCBI_API_KEY" in .Renviron file.
#' @param usehistory T/F. save result in your search history. you can also get WebEnv information. (default=T)
#' @param retstart (optional) number of skipping from the search  (Default=0)
#' @param retmax   (optional) total number of unique identifiers from the retrived set to be shown in the output (default=20)
#' @param rettype  (optional)  (optional) "uilist" or "count" (default="uilist"). if the value is count, the function only redirects the total number of query results.
#' @param retmode  (optional)  determine if the output is XML/plain text. the option will override user input if there is only one option is available.
#' @param datetype (optional) one from c("crdt","edat","pdat","mhda")
#' @return list or vector object.
#' @import tidyverse
#' @import xml2
#' @import stringi
#' @details
#' | Databases   | Record Type                     | retmode  | rettype  |
#' | ----------- | ------------------------------- | -------- | -------- |
#' | All	       | Document Summary	               | xml	    | docsum   |
#' | All	       | List of UIDs in XML             | xml	    | uilist   |
#' | All	       | List of UIDs in plain text      | text  	  | uilist   |
#' | MeSH	       | Full record in plain text	     | text	    | full     |
#' | NLM Catalog | Full record in plain text	     | text	    |          |
#' | NLM Catalog | XML	                           | xml	    |          |
#' | PubMed	     | ASN.1 (text)	                   | asn.1	  |          |
#' | PubMed	     | PubMed XML                      | xml	    |          |
#' | PubMed	     | MEDLINE (text)	                 | text	    | medline  |
#' | PubMed	     | PMID list (text)	               | text	    | uilist   |
#' | PubMed	     | Abstract (text)	               | text	    | abstract |
#' | PubMed	     |Summary (text)                   | text     | docsum   |
#'
#' | db         | format            | mode  |  Report Type                         |
#' | ---------- | ----------------- | ----- | ------------------------------------ |
#' |  (all)     |                   |       |                                      |
#' |            |     docsum        |       | DocumentSummarySet XML               |
#' |            |     docsum        | json  | DocumentSummarySet JSON.             |
#' |            |     full          |       | Same as native except for mesh       |
#' |            |     uid           |       | Unique Identifier List               |
#' |            |     url           |       | Entrez URL                           |
#' |            |     xml           |       | Same as -format full -mode xml       |
#' |            |                   |       |                                      |
#' | bioproject |                   |       |                                      |
#' |            |     native        |       | BioProject Report                    |
#' |            |     native        | xml   | RecordSet XML.                       |
#' |            |                   |       |                                      |
#' | biosample  |                   |       |                                      |
#' |            |     native        |       | BioSample Report                     |
#' |            |     native        | xml   | BioSampleSet XML                     |
#' | biosystems |                   |       |                                      |
#' |            |     native        | xml   | Sys-set XML                          |
#' |            |                   |       |                                      |
#' | clinvar    |                   |       |                                      |
#' |            |     variation     |       | Older Format                         |
#' |            |     variationid   |       | Transition Format                    |
#' |            |     vcv           |       | VCV Report                           |
#' |            |     clinvarset    |       | RCV Report                           |
#' | gds        |                   |       |                                      |
#' |            |     native        | xml   | RecordSet XML                        |
#' |            |     summary       |       | Summary                              |
#' | gene       |                   |       |                                      |
#' |            |    full_report    |       | Detailed Report                      |
#' |            |    gene_table     |       | Gene Table                           |
#' |            |    native         |       | Gene Report                          |
#' |            |    native         | asn.1 | Entrezgene ASN.1                     |
#' |            |    native         | xml   | Entrezgene-Set XML                   |
#' |            |    tabular        |       | Tabular Report                       |
#' | homologene |                   |       |                                      |
#' |            |   alignmentscores |       | Alignment Scores                     |
#' |            |   fasta           |       | FASTA                                |
#' |            |   homologene      |       | Homologene Report                    |
#' |            |   native          |       | Homologene List                      |
#' |            |   native          | asn.1 | HG-Entry ASN.1                       |
#' |            |   native          | xml   | Entrez-Homologene-Set XML            |
#' |            |                   |       |                                      |
#' | mesh       |                   |       |                                      |
#' |                   full                        Full Record
#' |                  native                      MeSH Report
#' |                 native             xml      RecordSet XML
#' |
#' | nlmcatalog
#' |                 native                      Full Record
#' |                  native             xml      NLMCatalogRecordSet XML
#' |
#' | pmc
#' |                bioc                        PubTator Central BioC XML
#' |                medline                     MEDLINE
#' |                native             xml      pmc-articleset XML
#' |
#' | pubmed
#' |                abstract                    Abstract
#' |                bioc                        PubTator Central BioC XML
#' |                medline                     MEDLINE
#' |                native             asn.1    Pubmed-entry ASN.1
#' |                native             xml      PubmedArticleSet XML
#' |
#' | (sequences)
#' |                acc                         Accession Number
#' |                est                         EST Report
#' |                fasta                       FASTA
#' |                fasta              xml      TinySeq XML
#' |                fasta_cds_aa                FASTA of CDS Products
#' |                fasta_cds_na                FASTA of Coding Regions
#' |                ft                          Feature Table
#' |                gb                          GenBank Flatfile
#' |                gb                 xml      GBSet XML
#' |                gbc                xml      INSDSet XML
#' |                gene_fasta                  FASTA of Gene
#' |                gp                          GenPept Flatfile
#' |                gp                 xml      GBSet XML
#' |                gpc                xml      INSDSet XML
#' |                gss                         GSS Report
#' |                ipg                         Identical Protein Report
#' |                ipg                xml      IPGReportSet XML
#' |                native             text     Seq-entry ASN.1
#' |                native             xml      Bioseq-set XML
#' |                seqid                       Seq-id ASN.1
#' | snp
#' |                  json                        Reference SNP Report
#' | sra
#' |                 native             xml      EXPERIMENT_PACKAGE_SET XML
#' |                  runinfo            xml      SraRunInfo XML
#' | structure
#' |                mmdb                        Ncbi-mime-asn1 strucseq ASN.1
#' |                  native                      MMDB Report
#' |                  native             xml      RecordSet XML
#' |
#' | taxonomy
#' |                  native                       Taxonomy List
#' |                  native             xml       TaxaSet XML

efetch <- function(id=NA,input=NA,
                   query_key,
                   WebEnv,
                   db="pubmed",
                   format=c("docsum","full","uilist","summary","gene_table",
                            "gb","gbc","ft","gbwithparts","fasta_cds_na",
                            "fasta_cds_aa","gp","gpc","ipg"),
                   retmode=c("xml","text","asn.1")) {
  message("temporary function")
  ## variables
  if(is.na(id) & is.na(input)) {message("you need id or input file");exit();}

  if(is.na(id)) {
  if(file.exists(input)) {
    if(is.na(id)) {id <-scan(input,what="character")} else {id<-c(id,scan(input,what="character"))}
  }}
  id<-paste0(id,collapse=",")
  baseurl <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch?"
  URL <- paste0(baseurl,"db=",db)
  URL <- paste0(URL,"&id=",id)
  URL <- paste0(URL,"&WebEnv=",id)

  if(db=="bioproject") {rettype <- c(rettype)}
  else if(db=="biosample") {psb_format <- "full";psb_retmode <-c("text","xml")}
  else if(db=="bioproject") {psb_format <- "full";psb_retmode <-c("xml")}
  else if(db=="gds") {psb_format <- "summary";psb_retmode <-"text"}
  else if(db=="gene") {
    psb_format<-c(NA,"gene_table")
    if(is.na(format)) {format <- psb_format[1] }
    if(is.na(format)==NA) {psb_retmode <- c("asn.1","xml")}
    else if(format=="gene_table") {psb_retmode <- c("text")}
  } else if(db=="homologene") {
    psb_format<-c(NA,"alignmentscore","fasta","homologene")
    if(is.na(format)) {format <- psb_format[1] }
    if(is.na(format)) {psb_retmode <- c("asn.1","xml")}
    else if(format=="alignmentscore") {psb_retmode <-"text"}
    else if(format=="fasta") {psb_retmode <-"text"}
    else if(format=="homologene") {psb_retmode <-"text"}
  } else if(db=="mesh") {psb_format <- "full";possibleretmode <-"text"}
  else if(db=="nlmcatalog") {psb_format="null";retmode="text"}
  else if(db=="nuccore") {
    psb_format <-c("null","native","acc","fasta","seqid", ## common in nuccore, protein or popset
                   "gb","gbc","ft","gbwithparts","fasta_cds_na","fasta_cds_aa")
    if(format=="null") {psb_retmode <-c("text","asn.1")}
    else if(format == 'native') {psb_retmode <- c("xml")}
    else if(format == 'acc') {psb_retmode <- c("text")}
    else if(format == 'fasta') {psb_retmode <- c("text","xml")}
    else if(format == 'seqid') {psb_retmode <- c("text")}
    else if(format == 'gb') {psb_retmod <- c("text","xml")}
    else if(format == 'gbc') {psb_retmod <- c("xml")}
  } else if(db=="protein"){
    psb_format <-c("null","native","acc","fasta","seqid","gb","gbc")
    if(format=="null") {psb_retmode <-c("text","asn.1")}
    else if(format == 'native') {psb_retmode <- c("xml")}
    else if(format == 'acc') {psb_retmode <- c("text")}
    else if(format == 'fasta') {psb_retmode <- c("text","xml")}
    else if(format == 'seqid') {psb_retmode <- c("text")}
    else if(format == 'gb') {psb_retmod <- c("text","xml")}

  }
  else if(db=="popset") {
    if(format=="null") {psb_retmode <-c("text","asn.1")}
    else if(format == 'native') {psb_retmode <- c("xml")}
    else if(format == 'acc') {psb_retmode <- c("text")}
    else if(format == 'fasta') {psb_retmode <- c("text","xml")}
    else if(format == 'seqid') {psb_retmode <- c("text")}
    else if(format == 'gb') {psb_retmod <- c("text","xml")}

  }

  else {if(is.na(psb_retmode)) { retmode <- psb_retmode[1]}}
  URL <- paste0(URL,"&format=",format)
  URL <- paste0(URL,"&retmod=",retmode)
  #if(retmode=="xml") {search_result<-xml2::read_xml(URL)}
  return(URL)
}


#' pubmed2tibble
#'
#' this function received pubmed ids to generated tibble df
