#PFAM2PDB.R
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


map_all_PDB_structures_to_pfam_sth_file_using_anyAtom_Ca_Cb <- function(pfamid,
                                                                        pfamSthFile = NA,
                                                                        sth.local.db = "data-raw/STH/", pfamseqResource = "rp75",
                                                                        max2take=1000) {
  if(is.na(pfamSthFile)) {
    pfamSthFile = sprintf("%s/%s_%s.sth", sth.local.db, pfamid, pfamseqResource)
    cat(sprintf("STH file was not given. Downloading to %s.\n", pfamSthFile))
    download.file(url = sprintf("http://pfam.xfam.org/family/%s/alignment/%s", pfamid, pfamseqResource),
                  destfile = pfamSthFile)
  }

  pfmmpdb <- choose_PDBs_for_a_pfamid(pfamid = pfamid, pfamSthFile = pfamSthFile)
  if(class(pfmmpdb) == "data.frame") {
    if(dim(pfmmpdb)[1] == 0) {
      return(NA)
    }
  }

  if(length(pfmmpdb) == 0) {
    return(NA)
  } else if(is.na(pfmmpdb)) {
    return(NA)
  }
  sthpdbmaps <- lapply(1:min(dim(pfmmpdb)[1], max2take),
                       function(x) {
                         cat(sprintf("Mapping pdb number: %d...\n", x))
                         get_PFAM2PDB_pairwise_distances_between_anyAtom_Ca_Cb(pfamid = pfamid, pfamSthFile = pfamSthFile,
                                                                               pfmmpdb = pfmmpdb[x,])
                       }
  )
  return(list(sthpdbmaps = sthpdbmaps, pfmmpdb = pfmmpdb))

}

get_PFAM2PDB_pairwise_distances_between_anyAtom_Ca_Cb <- function(pfamid,
                                                                  pfamSthFile = NA,
                                                                  pfmmpdb,
                                                                  sth.local.db = "data-raw/STH/", pfamseqResource = "rp75",
                                                                  minimalOverlap=0.80) {

  if(is.na(pfamSthFile)) {
    pfamSthFile = sprintf("%s/%s_%s.sth", sth.local.db, pfamid, pfamseqResource)
    cat(sprintf("STH file was not given. Downloading to %s.\n", pfamSthFile))
    download.file(url = sprintf("http://pfam.xfam.org/family/%s/alignment/%s", pfamid, pfamseqResource),
                  destfile = pfamSthFile)
  }


  sth <- read.stockholm.alignment(pfamSthFile)
  pdbEntry <- get_pdb_atom_coordinates(pfmmpdb$PDB_ID[1], pfmmpdb$CHAIN_ID[1])
  if(length(pdbEntry) == 0) {
    cat("ERROR - pdb entry is NA.\n")
    return(NA)
  } else if(is.na(pdbEntry)) {
    cat("ERROR - pdb entry is NA.\n")
    return(NA)
  }

  atoms <- pdbEntry$atoms
  #return(atoms)
  #INDEXING:
  atomsidx <- which(as.numeric(names(atoms)) %in% pfmmpdb$PdbResNumStart[1]:pfmmpdb$PdbResNumEnd[1])

  if(length(atomsidx) == 0) {
    cat("ERROR - pdb entry overlap is zero.\n")
    return(NA)
  }
  ols <- length(atomsidx)/length(pfmmpdb$PdbResNumStart[1]:pfmmpdb$PdbResNumEnd[1])
  cat(sprintf("Overlap fraction is: %f.\n", ols))
  if(ols < minimalOverlap) {
    cat("ERROR - pdb entry overlap is less than minimal.\n")
    return(NA)
  }
  resnoidx <- as.numeric(names(atoms))[atomsidx] #atoms$resno[atomsidx]
  s.pdb <- paste(aaa2a[sapply(atomsidx, function(x) atoms[[x]]$resid[1]),], collapse = "")
  #pdb inx: 34 56789 10 <--- pdbidx <- data.frame(alnidx = which(pdbres != '-'), pdbidx = resnoidx, pdbres = pdbres[pdbres != '-'])
  #pdb:     YL-EEAGAY
  #aln idx: 123456789
  #sth:     YLPE--GAY <--- sthidx <- data.frame(alnidx = which(sthres != '-'), sthidx = 1:nchar(s.sth), sthres = sthres[sthres != '-'])
  #dca idx: 1234  789

  # Preparation Pfam MSA domain sequence - s.sth:
  # This as in DCA definition excludes only '.' and lower case letters
  ####################################################################
  s.pfam <- sth[pfmmpdb$id, "seq"]
  s.pfam <- gsub("-", "*", s.pfam) #differeing Hmmer gap from pairwise alignment gap
  chars <- strsplit(s.pfam, "")[[1]]
  #excluded.idxs <- which(!(grepl("[[:upper:]]+", chars) & chars != '.')) # tmp!!!! this one is false
  excluded.idxs <- which(grepl("[[:lower:]]+", chars) | chars == '.') # | chars == '-')
  #idxs=which(!(grepl("[[:upper:]]+", chars) & chars != '.'))
  #return(list(idxs.org=idxs, idxs.new=excluded.idxs, chars=chars))
  s.sth <- paste(chars[-excluded.idxs], collapse = "")

  # Preparation of the assigned domain indices seq to the actual pdb - sthdomainidx:
  # This include all sequence letters (even lower case)
  ##############################################################################
  #               12345  6789
  # STH SEQ:      AB**A  RRRKKE*EE  S
  # DOMAIN SEQ:   AB  AwqRRRKKE EEppS
  #               12  3456789

  sthdomainidx <- rep(NA, length(chars[-excluded.idxs])) #length(chars[-excluded.idxs]) equals nchar(s.sth)
  gapsidxs <- which(chars[-excluded.idxs] == '*')
  domainseq <- chars[-which(chars == '.' | chars == '*')] #includes lower case letters but no gaps or dots
  domainidx <- (1:length(domainseq))[which(!grepl("[[:lower:]]+", domainseq))]
  if(length(gapsidxs) > 0) {
    sthdomainidx[-gapsidxs] <- domainidx
  } else {
    sthdomainidx <- domainidx
  }

  #return(list(sthdomainidx=sthdomainidx, domainseq=domainseq, sthseq = chars[-excluded.idxs],
  #            gapsidxs=gapsidxs, domainidx=domainidx, sthidx=(1:nchar(s.sth))))


  #do pairwise alignment:
  #######################

  #X Padding a prior to alignment in order to avoid tails clipping in global MSA:
  aln <- Biostrings::pairwiseAlignment(paste("XXXXXX", s.pdb, "XXXXXX", sep = "", collapse = ""),
                           paste("XXXXXX", s.sth, "XXXXXX", sep = "", collapse = ""),
                           type="global")

  pattern <- strsplit(as.character(aln@pattern), "")[[1]]
  pattern <- pattern[-c(1:6, length(pattern):(length(pattern)-5))]

  subject <- strsplit(as.character(aln@subject), "")[[1]]
  subject <- subject[-c(1:6, length(subject):(length(subject)-5))]

  pdbres <- pattern
  pdbidx <- data.frame(alnidx = which(pdbres != '-'), pdbidx = resnoidx, pdbres = pdbres[pdbres != '-'])
  sthres <- subject
  sthidx <- data.frame(alnidx = which(sthres != '-'), sthidx = 1:length(which(sthres != '-')),
                       domainidx = sthdomainidx, sthres = sthres[sthres != '-'])

  sthpdbidx <- merge(sthidx, pdbidx)
  #return(list(gapsidxs=gapsidxs, sthdomainidx=sthdomainidx, sthidx=sthidx, pdbidx=pdbidx, sthpdbidx=sthpdbidx, domainidx=domainidx, resnoidx=resnoidx,
  #            aln=aln, pfmmpdb=pfmmpdb, s.pfam=s.pfam, s.sth=s.sth, s.pdb = s.pdb, atoms=atoms, pattern = pattern, subject = subject))
  cat("Doing pairwise minimal atom distance calculations.\n")
  pdm <- get_minimal_pairwise_dist_between_atoms_Ca_Cb(atoms[atomsidx])
  rownames(sthpdbidx) <- sthpdbidx$pdbidx
  pdm$domainseq.idx.i <- sthpdbidx[ as.character(pdm$pdb.resno.i), "domainidx"]
  pdm$domainseq.idx.j <- sthpdbidx[ as.character(pdm$pdb.resno.j), "domainidx"]
  pdm$sth.idx.i <- sthpdbidx[ as.character(pdm$pdb.resno.i), "sthidx"]
  pdm$sth.idx.j <- sthpdbidx[ as.character(pdm$pdb.resno.j), "sthidx"]

  return(list(pdm = pdm, sthpdbidx = sthpdbidx, aln = aln,
              annot = list(pdb = pdbEntry$annot, uniprot = pfmmpdb)))

  #   map <- merge(sthpdbidx, atoms, by.x="pdbidx", by.y="resno")
  #   return(list(map = map, annot = list(pdb = pdbEntry$annot, uniprot = pfmmpdb)))
}


choose_PDBs_for_a_pfamid <- function(pfamid = "PF00075", pfamSthFile = NA,
                                     resolution.th = 3.0, experimentalTech = "X-RAY DIFFRACTION",
                                     pfamseqResource = "rp75",
                                     pdb.local.db = "data-raw/PDB/", sth.local.db = "data-raw/STH/") {

  #Get PDB assigned to Pfam:
  myPFmap <- hmmer_pdb_all[grep(pfamid, hmmer_pdb_all$PFAM_ACC), ]
  if(dim(myPFmap)[1] == 0) {
    cat(sprintf("%s has 0 PDBs. Quiting.\n", pfamid))
    return(NA)
  }
  myPFmap$PDB_ID <- tolower(myPFmap$PDB_ID)
  mm <- merge(myPFmap, pdb_chain_uniprot, by.x=c("PDB_ID", "CHAIN_ID"), by.y = c("PDB", "CHAIN"))

  if(is.na(pfamSthFile)) {
    pfamSthFile = sprintf("%s/%s_%s.sth", sth.local.db, pfamid, pfamseqResource)
    cat(sprintf("STH file was not given. Downloading to %s.\n", pfamSthFile))

    download.file(url = sprintf("http://pfam.xfam.org/family/%s/alignment/%s", pfamid, pfamseqResource),
                  destfile = pfamSthFile)

  }
  pf <- get_uniprot_accs_and_ids_from_pfam_sth_file(pfamSthFile = pfamSthFile)
  pfmm <- merge(pf, mm, by.x="acc", by.y="SP_PRIMARY")
  if(dim(pfmm)[1] == 0) {
    cat(sprintf("%s has 0 PDBs with uniprot mapping. Quiting.\n", pfamid))
    return(NA)
  }

  pdbAnnots <- NULL
  cntr <- 0
  #This is done in a for loop since the main delay is the download process
  for(pdbcode in unique(pfmm$PDB_ID)) {
    cntr <- cntr + 1
    fname <- bio3d::get.pdb(as.character(sprintf("%s", pdbcode)),
                     path=pdb.local.db, overwrite=F, gzip=T)
    if(fname == "1") {
      cat(sprintf("Could not download file for pdbcode: %s.\n", pdbcode))
      next
    }
    if(!file.exists(fname)) {
      cat(sprintf("Cannot find file: %s.\n", fname))
      next
    }
    cat(sprintf("Processing PDB file: %s (%d/%d).\n", fname, cntr, length(unique(pfmm$PDB_ID))))
    pdb <- bio3d::read.pdb(fname)
    pdb.annot <- bio3d::pdb.annotate(pdbcode)
    pdbAnnots <- rbind(pdbAnnots, pdb.annot)
  }
  if(length(pdbAnnots) == 0 | is.null(pdbAnnots)) {
    cat(sprintf("%s has 0 PDB annotations. Quiting.\n", pfamid))
    return(NA)
  }
  pdbAnnots$structureId <- tolower(pdbAnnots$structureId)
  pfmmpdb <- merge(pfmm, pdbAnnots, by.x = c("PDB_ID", "CHAIN_ID"), by.y = c("structureId", "chainId"))
  pfmmpdb <- pfmmpdb[which(as.numeric(pfmmpdb$resolution) <= resolution.th & pfmmpdb$experimentalTechnique == experimentalTech), ]

  return(pfmmpdb)
}

get_uniprot_accs_and_ids_from_pfam_sth_file <- function(pfamSthFile) {
  fname <- pfamSthFile
  text <- readLines(fname)
  accs <- unlist(lapply(text[grep("#=GS.+AC", x=text, perl=T)],
                        function(x) gsub("\\.\\d+", "", tail(strsplit(x, " ")[[1]], n=1))))
  ids <- unlist(lapply(text[grep("#=GS.+AC", x=text, perl=T)],
                       function(x) strsplit(x, " ")[[1]][2]))
  return(data.frame(acc=accs, id=ids, stringsAsFactors = F))
}


