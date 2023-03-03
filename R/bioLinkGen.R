#' bioLinkGen
#'
#' this function generates links (bioproject, run, experiment, sample information for) for a given strings.
#'
#' @param x a vector file.
#' @param format text or html a tag.
#' @param when format is "html", 4 target link action works. (_self, _blank, _parent, or _top). please see https://www.w3schools.com/tags/att_a_target.asp
#' @param desc TRUE or FALSE. describe the accession.
#' @import tidyverse
#' @return list object which contains summary and reports for denoted qc reports.
#' @examples
#' bioLinkGen("SRR17041298")
#' browseURL(bioLinkGen("SRR17041298"))
#' @export

bioLinkGen <- function(x,format="text",target="_blank",class="auto",desc=FALSE) {
  ## References
  ## genbank: [ ] https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
  ## Refseq : [x] https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
  if(class=="auto") {
    class <- case_when(substr(x,1,3)=="SRR"                              ~ "Run",          ## SRA Run
                    substr(x,1,3)=="SRX"                                 ~ "Experiment",   ## SRA Experiment
                    substr(x,1,3)=="PRJ"                                 ~ "BioProject",   ## SRA Bioproject
                    substr(x,1,4)=="SAMN"                                ~ "biosample",    ## Biosample
                    substr(x,1,3) %in% c("GCA","GCF")                    ~ "assembly",     ## Assembly accession
                    substr(x,1,3) %in% c("AC_","NC_","NG_")              ~ "RefSeq",       ## Genomic molecule for alternative assembly, referencem and incomplete genomic region
                    substr(x,1,3) %in% c("NT_","NW_")                    ~ "RefSeq",       ## Contig or scaffold.
                    substr(x,1,3) %in% c("NZ_")                          ~ "RefSeq",       ## Complete genomes and unfinished WGS data
                    substr(x,1,3) %in% c("NM_","XM_")                    ~ "RefSeq",       ## Protein-coding transcripts for curation and predicted model respectly.
                    substr(x,1,3) %in% c("NR_","XR_")                    ~ "RefSeq",       ## Non-protein-coding transcripts/	Predicted model non-protein-coding transcript
                    substr(x,1,3) %in% c("AP_","NP_","YP_","XP_","WP_")  ~ "RefSeq",       ## Annotated on AC_ alternate assembly, Associated with an NM_ or NC_ accession, Annotated on genomic molecules without an instantiated, (Predicted model) associated with an XM_ accession, Non-redundant across multiple strains and species

                  )
  }
  ## error handling

   stopifnot(class %in% c("Run","Experiment","BioProject","biosample","assembly","RefSeq"))
   html <- case_when(
     class == "Run"        ~ paste0('https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=',x),
     class == "Experiment" ~ paste0('https://www.ncbi.nlm.nih.gov/sra/',x),
     class == "BioProject" ~ paste0('https://www.ncbi.nlm.nih.gov/bioproject/',x),
     class == "biosample"  ~ paste0('https://www.ncbi.nlm.nih.gov/biosample/',x),
     class == "assembly"   ~ paste0('https://www.ncbi.nlm.nih.gov/assembly/',x),
     class == "RefSeq"     ~ paste0('https://www.ncbi.nlm.nih.gov/refseq/',x)
   )
  if(format=="text") {return(html)}

   x <- case_when(
      class == "Run"        ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
      class == "Experiment" ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
      class == "BioProject" ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
      class == "biosample"  ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
      class == "assembly"   ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>'),
      class == "RefSeq"     ~ paste0('<a href="',html,'" target="',target,'">',x,'</a>')
    )

  if(format=="html") {return(x)}
}




#html_bioPrjLnkGen <- function(x) {
#  return(sprintf('<a href="https://www.ncbi.nlm.nih.gov/bioproject/%s" target="_blank">%s</a>',x,x))
#}

#html_bioSmpLnkGen <- function(x) {
#  x<-gsub("SAMN","",x)
#  return(sprintf('<a href="https://www.ncbi.nlm.nih.gov/biosample//%s" target="_blank">%s</a>',x,x))
#}

#html_RunLnkGen <- function(x) {
#  return(sprintf('<a href="https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=%s" target="_blank">%s</a>',x,x))
#}
