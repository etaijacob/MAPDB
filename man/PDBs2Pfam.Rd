% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAPDB.R
\name{PDBs2Pfam}
\alias{PDBs2Pfam}
\title{Map all known PDB residues' coordinates to an Pfam STH file.}
\usage{
PDBs2Pfam(
  pfamid,
  pfamSthFile = NA,
  max2take = 12,
  resolution.th = 3,
  sth.local.db = "/tmp/",
  pdb.local.db = "/tmp/",
  pfamseqResource = "rp75"
)
}
\arguments{
\item{pfamid}{A Pfam id string}

\item{pfamSthFile}{A pfam sth file for pfamid - if not given will be downloaded to sth.local.db.}

\item{max2take}{maximal number of pdb entries to associate with the sth sequences with the given resolution threshold.}

\item{resolution.th}{minimal resolution to include in analysis}

\item{sth.local.db}{A path to a local pfam database or temporary directory to save downloaded sth files.}

\item{pdb.local.db}{A path to a local pdb database or temporary directory to save downloaded pdb files.}

\item{pfamseqResource}{A sequence redundancy level as inidcated in Pfam - relevant when downloading the sth file.}
}
\value{
A map of all known PDB residues to the msa file by PDB chain. TODO: add description
}
\description{
Map all known PDB residues' coordinates to an Pfam STH file.
}
\examples{
#In this example we download all relevant files for PF00075,
#including PDB files and sth alignment for the pfam domain.
#We limit the download to a total of 12 pdb instances with a resolution below 3.0 Ang.
#If any of the pdb files was already downloaded and located in the pdb.local.db,
#then it will be used to avoid
#unnecessary download.
myPfamMap <- PDBs2Pfam(pfamid = "PF00075", sth.local.db = "/tmp/",
                       max2take = 2, resolution.th = 3.0)

#In this example we have already downloaded the sth file
#therefore including it in the call for the function in order
#to avoid unnecessary download

myPfamMap <- PDBs2Pfam(pfamid = "PF00075", pfamSthFile = "/tmp/PF00075_rp75.sth",
                       sth.local.db = "/tmp/", max2take = 2, resolution.th = 3.0)
}
