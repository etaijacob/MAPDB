#PFAM2PDB.R

# tmp2 <- map_all_PDB_structures_to_pfam_sth_file_using_anyAtom_Ca_Cb(pfamid = "PF00075", max2take = 3,
#                                                                     pfamSthFile = "/tmp/PF00075_rp75.sth",
#                                                                     sth.local.db = "/tmp/", pdb.local.db = "/tmp/")
map_all_PDB_structures_to_pfam_sth_file_using_anyAtom_Ca_Cb <- function(pfamid,
                                                                        pfamSthFile = NA,
                                                                        sth.local.db, pdb.local.db,
                                                                        pfamseqResource = "rp75",
                                                                        resolution.th = 3.0,
                                                                        max2take=12) {
  pfmmpdb <- NULL
  if(is.na(pfamSthFile)) {
    pfamSthFile = sprintf("%s/%s_%s.sth", sth.local.db, pfamid, pfamseqResource)
    cat(sprintf("STH file was not given. Downloading to %s.\n", pfamSthFile))
    download.file(url = sprintf("http://pfam.xfam.org/family/%s/alignment/%s", pfamid, pfamseqResource),
                  destfile = pfamSthFile)
  }


  try(
    pfmmpdb <- choose_PDBs_for_a_pfamid(pfamid = pfamid, pfamSthFile = pfamSthFile, max2take = max2take,
                                        resolution.th = resolution.th,
                                        sth.local.db = sth.local.db, pdb.local.db = pdb.local.db))
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
  sthpdbmaps <- lapply(1:min(nrow(pfmmpdb), max2take),
                       function(x) {
                         cat(sprintf("Mapping pdb number: %d...\n", x))
                         get_PFAM2PDB_pairwise_distances_between_anyAtom_Ca_Cb(pfamid = pfamid, pfamSthFile = pfamSthFile,
                                                                               pdbLocalPath = pdb.local.db,
                                                                               pfmmpdb = pfmmpdb[x,])
                       }
  )
  return(list(sthpdbmaps = sthpdbmaps, pfmmpdb = pfmmpdb))

}

get_PFAM2PDB_pairwise_distances_between_anyAtom_Ca_Cb <- function(pfamid,
                                                                  pfamSthFile = NA,
                                                                  pfmmpdb, pdbLocalPath,
                                                                  sth.local.db = "data-raw/STH/", pfamseqResource = "rp75",
                                                                  minimalOverlap=0.80) {

  if(is.na(pfamSthFile)) {
    pfamSthFile = sprintf("%s/%s_%s.sth", sth.local.db, pfamid, pfamseqResource)
    cat(sprintf("STH file was not given. Downloading to %s.\n", pfamSthFile))
    download.file(url = sprintf("http://pfam.xfam.org/family/%s/alignment/%s", pfamid, pfamseqResource),
                  destfile = pfamSthFile)
  }


  sth <- read.stockholm.alignment(pfamSthFile)
  pdbEntry <- get_pdb_atom_coordinates(pfmmpdb$PDB_ID[1], pfmmpdb$CHAIN_ID[1], pdbLocalPath =  pdbLocalPath)
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

#' Download and annotate the PDB files associated with the PFAM domain sequences in the sth file
#'
#' @param pfamid A Pfam id string.
#' @param pfamSthFile A pfam sth file for pfamid.
#' @param max2take The maximal number of PDB files to download with the given resolution (default = 3.0).
#' @param resolution.th Minimal resoultion required for a solved structure to be considered.
#' @param experimentalTech Experimental method required for the structre (default is X-RAY).
#' @param pfamseqResource A sequence redundancy level as inidcated in Pfam - relevant when downloading the sth file.
#' @param pdb.local.db A directory path to a local pdb database or temporary directory to save downloaded sth files.
#' Can be used to accomulate files and avoid repeated downloads.
#' @param sth.local.db A directory path to a local pfam database or temporary directory to save downloaded sth files.
#' Can be used to accomulate files and avoid repeated downloads.
#' @return A map of all known PDB residues to the msa file by PDB chain. TODO: add description
#' @examples
#' #In this example we download all relevant relevant PDB files pfam domain PF00075.
#' #We limit the download to a total of 12 pdb instances with a resolution below 3.0 Ang.
#' #If any of the pdb files was already downloaded and located in the pdb.local.db,
#' #then it will be used to avoid
#' #unnecessary download.
#' tmp <- choose_PDBs_for_a_pfamid(pfamid = "PF00075", pfamSthFile = "/tmp/PF00075_rp75.sth",
#'                                 pdb.local.db = "/tmp/", sth.local.db = "/tmp/", max2take = 3)
#' @export
choose_PDBs_for_a_pfamid <- function(pfamid = "PF00075", pfamSthFile = NA, max2take = 12,
                                     resolution.th = 3.0, experimentalTech = "X-RAY DIFFRACTION",
                                     pfamseqResource = "rp75",
                                     pdb.local.db, sth.local.db) {

  #Get PDB assigned to Pfam:
  myPFmap <- hmmer_pdb_all[grep(pfamid, hmmer_pdb_all$PFAM_ACC), ]
  if(dim(myPFmap)[1] == 0) {
    cat(sprintf("%s has 0 PDBs. Quiting.\n", pfamid))
    return(NA)
  }
  myPFmap$PDB_ID <- tolower(myPFmap$PDB_ID)
  mm <- merge(myPFmap, pdbUniprot, by.x=c("PDB_ID", "CHAIN_ID"), by.y = c("PDB", "CHAIN"))

  if(is.na(pfamSthFile)) {
    pfamSthFile = sprintf("%s/%s_%s.sth", sth.local.db, pfamid, pfamseqResource)
    cat(sprintf("STH file was not given. Downloading to %s.\n", pfamSthFile))

    download.file(url = sprintf("http://pfam.xfam.org/family/%s/alignment/%s", pfamid, pfamseqResource),
                  destfile = pfamSthFile)

  }
  pf <- get_uniprot_entry_names_and_ids_from_pfam_sth_file(pfamSthFile = pfamSthFile)
  #return(list(pf = pf, mm = mm))
  #pfmm <- merge(pf, mm, by.x="acc", by.y="SP_PRIMARY")
  pfmm <- merge(pf, mm, by.x="prot_name", by.y="SP_PRIMARY")
  if(dim(pfmm)[1] == 0) {
    cat(sprintf("%s has 0 PDBs with uniprot mapping. Quiting.\n", pfamid))
    return(NA)
  }

  pdbAnnots <- NULL
  cntr <- 0
  pdbCounter <- 0
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
    try_ret1 <- try(pdb <- bio3d::read.pdb(fname))
    try_ret2 <- try(pdb.annot <- bio3d::pdb.annotate(pdbcode))
    if(try_ret1[1] != "try-error" & try_ret2[1] != "try-error") {
      pdbAnnots <- rbind(pdbAnnots, pdb.annot)
      message("Resolution: ", min(pdb.annot$resolution))
      if(sum(pdb.annot$resolution < resolution.th) > 0)
        pdbCounter <- pdbCounter + 1

    } else {
      message("Error in retrieving ", pdbcode, ".")
    }
    if(pdbCounter >= max2take) {
      message("Reached the max number of pdb entries to download (", max2take, ").")
      break
    }

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

get_uniprot_entry_names_and_ids_from_pfam_sth_file <- function(pfamSthFile) {
  fname <- pfamSthFile
  text <- readLines(fname)


  ids <- unname(sapply(text[-grep("^#=|^#", x = text)], function(x) strsplit(x = x, " +")[[1]][1]))
  prot_names <- unname(sapply(ids, function(x) strsplit(x, "/")[[1]][1]))
  prot_names <- unname(sapply(prot_names, function(x) strsplit(x, split = "\\.")[[1]][1]))
  # accs <- unlist(lapply(text[grep("#=GS.+AC", x=text, perl=T)],
  #                       function(x) gsub("\\.\\d+", "", tail(strsplit(x, " ")[[1]], n=1))))
  # ids <- unlist(lapply(text[grep("#=GS.+AC", x=text, perl=T)],
  #                      function(x) strsplit(x, " ")[[1]][2]))
  return(data.frame(prot_name = prot_names, id = ids, stringsAsFactors = F))
}


