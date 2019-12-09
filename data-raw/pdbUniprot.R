#pdbUniprot.R

#This source code should be executed only when an updated of the package is required.
#Currently for swissprot
updateUniprotPDBmaps <- function(rawdatadir = "data-raw/", datadir = "data/") {
  require(data.table)
  destfile1 <- sprintf("%s/pdb_chain_uniprot.lst", rawdatadir)
  download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst",
                destfile = destfile1)
  pdb_chain_uniprot <- fread(file = destfile1, sep = "\t", header = T, stringsAsFactors = F)

  IDACCmappingsURL <- "http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length&format=tab"
  destfile2 <- sprintf("%s/swissprot_idaccs.td", rawdatadir)
  download.file(url = IDACCmappingsURL, destfile = destfile2)
  swissprot_idaccs <- fread(file = destfile2, sep = "\t", header = T,
                            stringsAsFactors = F, check.names = T, quote = "")

  save(swissprot_idaccs, file = sprintf("%s/swissprot_idaccs.RData", rawdatadir))
  pdbUniprot <- merge(pdb_chain_uniprot, swissprot_idaccs, by.x="SP_PRIMARY", by.y="Entry")

  pdbUniprot <- pdbUniprot[, c("SP_PRIMARY", "PDB", "CHAIN", "RES_BEG", "RES_END", "PDB_BEG", "PDB_END", "SP_BEG", "SP_END", "Entry.name")]
  save(pdbUniprot, file = sprintf("%s/pdbUniprot.RData", datadir), compress = "xz")

}

updateUniprotPDBmaps()
