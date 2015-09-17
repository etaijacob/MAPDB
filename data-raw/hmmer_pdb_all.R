


updatePFAMPDBmaps <- function(rawdatadir = "data-raw/", datadir = "data/") {
  destfile1 <- sprintf("%s/hmmer_pdb_all.txt", rawdatadir)
  download.file(url = "http://www.rcsb.org/pdb/rest/hmmer?file=hmmer_pdb_all.txt",
                destfile = destfile1)
  hmmer_pdb_all <- read.table(file = destfile1, sep = "\t", header = T, stringsAsFactors = F)
  save(hmmer_pdb_all, file = sprintf("%s/hmmer_pdb_all.RData", datadir))

}

updatePFAMPDBmaps()
