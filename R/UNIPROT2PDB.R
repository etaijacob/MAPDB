#UNIPROT2PDB.R

updateSwissprotPDBmaps <- function(datadir = "data/") {
  destfile1 <- sprintf("%s/pdb_chain_uniprot.lst", datadir)
  download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst",
                destfile = destfile1)
  pdb_chain_uniprot <- read.table(file = destfile, sep = "\t", comment.char = '#', header = T, stringsAsFactors = F)
  save(pdb_chain_uniprot, file = sprintf("%s/pdb_chain_uniprot.RData", datadir))

  IDACCmappingsURL <- "http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length&format=tab"
  destfile2 <- sprintf("%s/swissprot_idaccs.td", datadir)
  download.file(url = IDACCmappingsURL, destfile = destfile2)
  swissprot_idaccs <- read.table(file = destfile2, sep = "\t", comment.char = '#', header = T, stringsAsFactors = F, check.names = T)
  save(swissprot_idaccs, file = sprintf("%s/swissprot_idaccs.RData", datadir))

}
