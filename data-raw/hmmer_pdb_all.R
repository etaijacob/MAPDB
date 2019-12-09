

#This source file should be executed only when update to package is required.
updatePFAMPDBmaps <- function(rawdatadir = "data-raw/", datadir = "data/") {
  require(data.table)
  destfile1 <- sprintf("%s/hmmer_pdb_all.txt", rawdatadir)
  download.file(url = "http://www.rcsb.org/pdb/rest/hmmer?file=hmmer_pdb_all.txt",
                destfile = destfile1)
  hmmer_pdb_all <- data.table::fread(file = destfile1, sep = "\t", header = T, stringsAsFactors = F)
  hmmer_pdb_all <- hmmer_pdb_all[, c("PDB_ID", "CHAIN_ID", "PdbResNumStart", "PdbResNumEnd", "PFAM_ACC")]
  save(hmmer_pdb_all, file = sprintf("%s/hmmer_pdb_all.RData", datadir), compress = "xz")

}

updatePFAMPDBmaps()
