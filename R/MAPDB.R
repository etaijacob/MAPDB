#MAPDB.R
# MAPDB - An R package
# Copyright (C) 2015  Etai Jacob
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#TODO: edit all descriptions.

#' MAPDB: A package for mapping Pfam multiple sequence alignment (MSA in sth format) or Uniprot (MSA in fasta format)
#' to PDB residue coordinates.
#'
#' MAPDB provides an easy way to map all known PDB chains to an STH or fasta MSA file.
#'
#' This package is used to calculate the accuracy of contact predictions based on MSA data.
#' If you use this package please cite:
#'
#' Etai Jacob, Ron Unger and Amnon Horovitz (2015). Codon-level information improves predictions of inter-residue contacts
#' in proteins by correlated mutation analysis.
#' eLife 2015;10.7554/eLife.08932 URL http://dx.doi.org/10.7554/eLife.08932.
#' It has two main functions:
#'
#' \itemize{
#' \item Map all known PDB residues' coordinates to an Pfam STH file
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
#' @import Biostrings
NULL

#' Map all known PDB residues' coordinates to an Pfam STH file.
#'
#' @param pfamid A Pfam id string
#' @param pfamSthFile A pfam sth file for pfamid.
#' @param sth.local.db A path to a local pfam database or temporary directory to save downloaded sth files.
#' @param pfamseqResource A sequence redundancy level as inidcated in Pfam - relevant when downloading the sth file.
#' @return A map of all known PDB residues to the msa file by PDB chain.
#' @examples
#' myPfamMap <- PDBs2Pfam("PF00075")
PDBs2Pfam <- function(pfamid,
                      pfamSthFile = NA,
                      sth.local.db = "data-raw/STH/", pfamseqResource = "rp75", ...) {
  res <- map_all_PDB_structures_to_pfam_sth_file_using_anyAtom_Ca_Cb(pfamid,
                                                                     pfamSthFile = NA,
                                                                     sth.local.db = sth.local.db,
                                                                     pfamseqResource = pfamseqResource,
                                                                     max2take=1000)
  return(res)
}


#' Mapping between Pfam domains and PDB chains.
#'
#' A dataset containing the mapping between Pfam domains (based on HMMER) and PDB for more than
#' 130,000 PDB chains.
#'
#' @format A data frame with 130678 rows and 8 variables:
#' \describe{
#'   \item{PDB_ID}{4 letters PDB code}
#'   \item{CHAIN_ID}{PDB chain id}
#'   \item{PdbResNumStart}{PDB residue number where the pfam domain starts}
#'   \item{PdbResNumEnd}{PDB residue number where the pfam domain ends}
#'   \item{PFAM_ACC}{Pfam id}
#'   ...
#' }
#' @source \url{http://www.rcsb.org/pdb/rest/hmmer?file=hmmer_pdb_all.txt}
"hmmer_pdb_all"

#' Mapping between PDB chains and uniprot entries.
#'
#' A dataset containing the mapping between PDB chains and Uniprot entries.
#'
#' @format A data frame with 294960 rows and 9 variables:
#' \describe{
#'   \item{PDB}{4 letters PDB code}
#'   \item{CHAIN}{PDB chain id}
#'   \item{RES_BEG}{Uniprot position where PDB chain starts}
#'   \item{RES_END}{Uniprot position where PDB chain ends}
#'   \item{PDB_BEG}{PDB residue number where the uniprot sequence starts}
#'   \item{PDB_END}{PDB residue number where the uniprot sequence ends}
#'   ...
#' }
#' @source \url{ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst}
"pdb_chain_uniprot"

#' Mapping between PDB chains and swissprot accs and ids.
#'
#' A processed dataset of almost 230,000 rows containing the mapping between PDB chains and swissprot entries.
#' The list is based on merging pdb_chain_uniprot.lst file from EBI and Uniprot. For more inforamtion see pdbUniprot.R in this
#' Package.
#'
#' @format A data frame with 229185 rows and 15 variables:
#' \describe{
#'   \item{SP_PRIMARY}{price, in US dollars}
#'   \item{PDB}{4 letters PDB code}
#'   \item{CHAIN}{PDB chain id}
#'   \item{RES_BEG}{Uniprot position where PDB chain starts}
#'   \item{RES_END}{Uniprot position where PDB chain ends}
#'   \item{PDB_BEG}{PDB residue number where the uniprot sequence starts}
#'   \item{PDB_END}{PDB residue number where the uniprot sequence ends}
#'   \item{Entry.name}{Entry name}
#'   \item{Gene.names}{List of genes of the entry (space separated)}
#'   ...
#' }
#' @source \url{http://www.uniprot.org/uniprot}
"pdbUniprot"


