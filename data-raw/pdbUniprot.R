#pdbUniprot.R

updateSwissprotPDBmaps()

updateSwissprotPDBmaps <- function(rawdatadir = "data-raw/", datadir = "data/") {
  destfile1 <- sprintf("%s/pdb_chain_uniprot.lst", rawdatadir)
  download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst",
                destfile = destfile1)
  pdb_chain_uniprot <- read.table(file = destfile, sep = "\t", comment.char = '#', header = T, stringsAsFactors = F)
  save(pdb_chain_uniprot, file = sprintf("%s/pdb_chain_uniprot.RData", rawdatadir))

  IDACCmappingsURL <- "http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length&format=tab"
  destfile2 <- sprintf("%s/swissprot_idaccs.td", rawdatadir)
  download.file(url = IDACCmappingsURL, destfile = destfile2)
  swissprot_idaccs <- read.csv2(file = destfile2, sep = "\t", comment.char = '#', header = T,
                                stringsAsFactors = F, check.names = T, quote = "")
  save(swissprot_idaccs, file = sprintf("%s/swissprot_idaccs.RData", rawdatadir))
  pdbUniprot <- merge(pdb_chain_uniprot, swissprot_idaccs, by.x="SP_PRIMARY", by.y="Entry")
  save(pdbUniprot, file = sprintf("%s/pdbUniprot.RData", datadir))

}
