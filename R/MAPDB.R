#MAPDB.R
#' MAPDB: A package for mapping Pfam multiple sequence alignment (MSA in sth format) or Uniprot (MSA in fasta format)
#' to PDB residue coordinates.
#'
#' MAPDB provides an easy way to map all known PDB chains to an STH or fasta MSA file.
#'
#' using this package requires an interent connection. If more then several examples are required,
#' you should consider downloading first PDB and STH databases.
#'
#' This package can be used, for example, to calculate the accuracy of contact predictions based on MSA data.
#' If you use this package please cite:
#'
#' Etai Jacob, Ron Unger and Amnon Horovitz (2015). Codon-level information improves predictions of inter-residue contacts
#' in proteins by correlated mutation analysis.
#' eLife 2015;10.7554/eLife.08932 URL http://dx.doi.org/10.7554/eLife.08932.
#' It has two main functions:
#'
#' \itemize{
#' \item Map all known PDB residues' coordinates to a Pfam STH file
#' \item Map all known PDB residues' coordinates to a fasta file (partially supported at this version)
#' }
#'
#' To learn more about MAPDB, start with the vignettes:
#' \code{browseVignettes(package = "MAPDB")}
#'
#' @docType package
#' @name MAPDB
#' @import bio3d
#' @import utils
#' @import pdist
#' @importFrom seqinr read.alignment
#' @rawNamespace import(Biostrings, except = c("mask", "tail", "head"))

# @import import(Biostrings, except = c("mask", "tail", "head"))
NULL

#' Map all known PDB residues' coordinates to an Pfam STH file.
#'
#' @param pfamid A Pfam id string
#' @param pfamSthFile A pfam sth file for pfamid - if not given will be downloaded to sth.local.db.
#' @param max2take maximal number of pdb entries to associate with the sth sequences with the given resolution threshold.
#' @param resolution.th minimal resolution to include in analysis
#' @param sth.local.db A path to a local pfam database or temporary directory to save downloaded sth files.
#' @param pdb.local.db A path to a local pdb database or temporary directory to save downloaded pdb files.
#' @param pfamseqResource A sequence redundancy level as inidcated in Pfam - relevant when downloading the sth file.
#' @return A map of all known PDB residues to the msa file by PDB chain. TODO: add description
#' @examples
#' #In this example we download all relevant files for PF00075,
#' #including PDB files and sth alignment for the pfam domain.
#' #We limit the download to a total of 12 pdb instances with a resolution below 3.0 Ang.
#' #If any of the pdb files was already downloaded and located in the pdb.local.db,
#' #then it will be used to avoid
#' #unnecessary download.
#' myPfamMap <- PDBs2Pfam(pfamid = "PF00075", sth.local.db = "/tmp/",
#'                        max2take = 2, resolution.th = 3.0)
#'
#' #In this example we have already downloaded the sth file
#' #therefore including it in the call for the function in order
#' #to avoid unnecessary download
#'
#' myPfamMap <- PDBs2Pfam(pfamid = "PF00075", pfamSthFile = "/tmp/PF00075_rp75.sth",
#'                        sth.local.db = "/tmp/", max2take = 2, resolution.th = 3.0)
#' @export
PDBs2Pfam <- function(pfamid,
                      pfamSthFile = NA, max2take = 12, resolution.th = 3.0,
                      sth.local.db = "/tmp/", pdb.local.db = "/tmp/",
                      pfamseqResource = "rp75") {
  res <- map_all_PDB_structures_to_pfam_sth_file_using_anyAtom_Ca_Cb(pfamid,
                                                                     pfamSthFile = NA,
                                                                     sth.local.db = sth.local.db,
                                                                     pdb.local.db = pdb.local.db,
                                                                     resolution.th = resolution.th,
                                                                     pfamseqResource = pfamseqResource,
                                                                     max2take=max2take)
  return(res)
}


#' Mapping between Pfam domains and PDB chains.
#'
#' A dataset containing the mapping between Pfam domains (based on HMMER) and PDB for more than
#' 167177 PDB chains.
#'
#' @format A data frame with 130678 rows and 8 variables:
#' \describe{
#'   \item{PDB_ID}{4 letters PDB code}
#'   \item{CHAIN_ID}{PDB chain id}
#'   \item{PdbResNumStart}{PDB residue number where the pfam domain starts}
#'   \item{PdbResNumEnd}{PDB residue number where the pfam domain ends}
#'   \item{PFAM_ACC}{Pfam id}
#' }
#' @source \url{http://www.rcsb.org/pdb/rest/hmmer?file=hmmer_pdb_all.txt}
"hmmer_pdb_all"


#' Mapping between PDB chains and swissprot accs and ids.
#'
#' A processed dataset of almost 230,000 rows containing the mapping between PDB chains and swissprot entries.
#' The list is based on merging pdb_chain_uniprot.lst file from EBI and Uniprot. For more inforamtion see pdbUniprot.R in this
#' Package.
#'
#' @format A data frame with 229185 rows and 15 variables:
#' \describe{
#'   \item{SP_PRIMARY}{id}
#'   \item{PDB}{4 letters PDB code}
#'   \item{CHAIN}{PDB chain id}
#'   \item{RES_BEG}{Uniprot position where PDB chain starts}
#'   \item{RES_END}{Uniprot position where PDB chain ends}
#'   \item{PDB_BEG}{PDB residue number where the uniprot sequence starts}
#'   \item{PDB_END}{PDB residue number where the uniprot sequence ends}
#'   \item{SP_BEG}{start position}
#'   \item{SP_END}{end position}
#'   \item{Entry.name}{Entry name}
#' }
#' @source \url{http://www.uniprot.org/uniprot}
"pdbUniprot"

if(getRversion() >= "2.15.1")  utils::globalVariables(c("pdbUniprot", "hmmer_pdb_all", "pfmmpdb", "aaa2a", "a2aaa"))
