#UNIPROT2PDB.R

updateSwissprotPDBmaps <- function(datadir = "data/") {
  destfile <- sprintf("%s/pdb_chain_uniprot.lst", datadir)
  download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst",
                destfile = destfile)
  pdb_chain_uniprot <- read.table(file = destfile, sep = "\t", comment.char = '#', header = T, stringsAsFactors = F)
  save(pdb_chain_uniprot, file = sprintf("%s/pdb_chain_uniprot.RData", datadir))
}
