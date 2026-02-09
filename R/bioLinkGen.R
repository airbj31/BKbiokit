#' bioLinkGen
#'
#' this function generates links (bioproject, run, experiment, sample information for) for a given strings.
#'
#' @param x a vector file.
#' @param format text or html link tag (<a href="blah blah>x</a>).
#' @param target when format is "html", 4 target link action works. (_self, _blank, _parent, or _top). please see https://www.w3schools.com/tags/att_a_target.asp. if the format is text, the target option is not used.
#' @param desc TRUE or FALSE. describe the accession.
#' @return the function vector object which contains text or html link which redirect the accessions.
#' @import dplyr
#' @examples
#' get links
#' bioLinkGen("SRR17041298")
#' bioLinkGen("PRJNA216922")
#' browseURL(bioLinkGen("SRR17041298"))
#' @export

bioLinkGen <- function(x,format="text",target="_blank",class="auto",desc=FALSE) {
  ## References
  ## genbank: [ ] https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
  ## Refseq : [x] https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
  ## Uniprot: [x] https://www.uniprot.org/help/accession_numbers/
  ## geo:     [x] https://www.ncbi.nlm.nih.gov/geo/info/

#  ifelse(is.na(x),return(x),)

  if(class=="auto") {
    class <- case_when(
                    is.na(x)                                                 ~ "None",
                    grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$",x,perl=TRUE) ~ "uniprotkb",       ## uniprot
                    substr(x,1,3) %in% c("SRR","ERR")                        ~ "Run",             ## SRA Run
                    substr(x,1,3)=="SRX"                                     ~ "Experiment",      ## SRA Experiment
                    substr(x,1,3)=="PRJ"                                     ~ "BioProject",      ## SRA Bioproject
                    substr(x,1,4)=="SAMN"                                    ~ "biosample",       ## Biosample
                    substr(x,1,3) %in% c("GCA","GCF")                        ~ "assembly",        ## Assembly accession
                    substr(x,1,3) %in% c("GDS","GSM","GPL")                  ~ "GEO",             ## GEO accession
                    substr(x,1,3) %in% c("AC_","NC_","NG_")                  ~ "RefSeq",          ## Genomic molecule for alternative assembly, referencem and incomplete genomic region
                    substr(x,1,3) %in% c("NT_","NW_")                        ~ "RefSeq",          ## Contig or scaffold.
                    substr(x,1,3) %in% c("NZ_")                              ~ "RefSeq",          ## Complete genomes and unfinished WGS data
                    substr(x,1,3) %in% c("NM_","XM_")                        ~ "RefSeq",          ## Protein-coding transcripts for curation and predicted model respectly.
                    substr(x,1,3) %in% c("NR_","XR_")                        ~ "RefSeq",          ## Non-protein-coding transcripts/	Predicted model non-protein-coding transcript
                    substr(x,1,3) %in% c("AP_","NP_","YP_","XP_","WP_")      ~ "RefSeq",          ## Annotated on AC_ alternate assembly, Associated with an NM_ or NC_ accession, Annotated on genomic molecules without an instantiated, (Predicted model) associated with an XM_ accession, Non-redundant across multiple strains and species
                    grepl("^[a-zA-Z]{1}[0-9]{5}",x,perl=TRUE)                ~ "genbank",         ## genbank nucleotides
                    grepl("^[a-zA-Z]{2}[0-9]{6}",x,perl=TRUE)                ~ "genbank",         ## genbank nucleotides
                    grepl("^[a-zA-Z]{2}[0-9]{8}",x,perl=TRUE)                ~ "genbank",         ## genbank nucleotides
                    grepl("^[a-zA-Z]{3}[0-9]{5}",x,perl=TRUE)                ~ "genbank_protein", ## genbank protein
                    grepl("^[a-zA-Z]{3}[0-9]{7}",x,perl=TRUE)                ~ "genbank_protein", ## genbank protein
                    grepl("^[a-zA-Z]{4}[0-9]{8+}",x,perl=TRUE)               ~ "genbank_WGS",     ## genbank WGS
                    grepl("^[a-zA-Z]{6}[0-9]{9+}",x,perl=TRUE)               ~ "genbank_WGS",     ## genbank WGS
                    grepl("^[a-zA-Z]{5}[0-9]{7}",x,perl=TRUE)                ~ "genbank_MGA",     ## genbank MGA
                    grepl("^10\\.[0-9]{4,}/\\S+$",x,perl=TRUE)             ~ "doi"              ## doi
                  )
  }
  ## error handling

   stopifnot(class %in% c("Run","Experiment","BioProject","biosample","assembly",
                          "RefSeq","None","genbak","genbank_protein","genbank_WGS",
                          "genbank_MGA",
                          "uniprotkb","doi","GEO","pubmed"))
   html <- case_when(
     class == "None"                       ~ as.character(NA),
     class == "Run"                        ~ paste0('https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=',x),
     class == "Experiment"                 ~ paste0('https://www.ncbi.nlm.nih.gov/sra/',x),
     class == "BioProject"                 ~ paste0('https://www.ncbi.nlm.nih.gov/bioproject/',x),
     class == "biosample"                  ~ paste0('https://www.ncbi.nlm.nih.gov/biosample/',x),
     class == "assembly"                   ~ paste0('https://www.ncbi.nlm.nih.gov/assembly/',x),
     class == "RefSeq"                     ~ paste0('https://www.ncbi.nlm.nih.gov/refseq/',x),
     class == "uniprotkb"                  ~ paste0('https://www.uniprot.org/uniprotkb/',x),
     class %in% c("genbank","genbank_WGS") ~ paste0('https://www.ncbi.nlm.nih.gov/nuccore/',x),
     class == "genbank_protein"            ~ paste0('https://www.ncbi.nlm.nih.gov/nuccore/',x),
     class == "doi"                        ~ paste0('https://www.doi.org/',x),
     class == "pubmed"                     ~ paste0('https://pubmed.ncbi.nlm.nih.gov/',x),
     class %in% c("GSE","GSM","GPL","GEO")       ~ paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',x),
   )
  if(format=="text") {return(html)}

   x <- case_when(
       class == "None"       ~ as.character(NA),
#      class == "Run"        ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
#      class == "Experiment" ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
#      class == "BioProject" ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
#      class == "biosample"  ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
#      class == "assembly"   ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
#      class == "RefSeq"     ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
       TRUE        ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>')
    )

  if(format=="html") {return(x)}
}
