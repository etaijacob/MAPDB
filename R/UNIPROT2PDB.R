#UNIPROT2PDB.R

get_UNIPROT2PDB_pairwise_distances_between_anyAtom_Ca_Cb <- function(uniprotseq,
                                                                     mypdb,
                                                                     minimalOverlap=0.80) {

  uniprotseq <- gsub(pattern = "\n", "", uniprotseq)
  #Get PDB sequences:
  pdbEntry <- get_pdb_atom_coordinates(mypdb$PDB[1], mypdb$CHAIN[1])
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
  atomsidx <- which(as.numeric(names(atoms)) %in% mypdb$PDB_BEG[1]:mypdb$PDB_END[1])

  if(length(atomsidx) == 0) {
    cat("ERROR - pdb entry overlap is zero.\n")
    return(NA)
  }
  ols <- length(atomsidx)/length(mypdb$PDB_BEG[1]:mypdb$PDB_END[1])
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

  # Preparation of sequence:
  # This as in DCA definition excludes only '.' and lower case letters
  ####################################################################

  uniprotseq <- gsub("-", "*", uniprotseq) #differeing input gaps from pairwise alignment gaps
  uniprotseqnogaps <- gsub("\\*", "", uniprotseq)
  chars <- strsplit(uniprotseq, "")[[1]]
  gapsidxs <- which(chars == '*')

  alignmentidx <- 1:length(chars)
  uniprotseqidx <- rep(NA, length(alignmentidx))
  if(length(gapsidxs) > 0) {
    uniprotseqidx[ -gapsidxs] <- 1:length(chars[-gapsidxs])

  } else {
    uniprotseqidx <- alignmentidx
  }


  #do pairwise alignment:
  #######################

  #X Padding a prior to alignment in order to avoid tails clipping in global MSA:
  aln <- Biostrings::pairwiseAlignment(paste("XXXXXX", s.pdb, "XXXXXX", sep = "", collapse = ""),
                           paste("XXXXXX", uniprotseqnogaps, "XXXXXX", sep = "", collapse = ""),
                           type="global")

  pattern <- strsplit(as.character(aln@pattern), "")[[1]]
  pattern <- pattern[-c(1:6, length(pattern):(length(pattern)-5))]

  subject <- strsplit(as.character(aln@subject), "")[[1]]
  subject <- subject[-c(1:6, length(subject):(length(subject)-5))]

  pdbres <- pattern
  pdbidx <- data.frame(alnidx = which(pdbres != '-'), pdbidx = resnoidx, pdbres = pdbres[pdbres != '-'])
  seqres <- subject
  seqidx <- data.frame(alnidx = which(seqres != '-'), uniprotidx = 1:length(which(seqres != '-')),
                       alignmentidx = alignmentidx[-gapsidxs], uniprotres = seqres[seqres != '-'])

  uniprotpdbidx <- merge(seqidx, pdbidx)

  #return(list(aln=aln, uniprotpdbidx = uniprotpdbidx))
  #return(list(gapsidxs=gapsidxs, sthdomainidx=sthdomainidx, sthidx=sthidx, pdbidx=pdbidx, sthpdbidx=sthpdbidx, domainidx=domainidx, resnoidx=resnoidx,
  #            aln=aln, pfmmpdb=pfmmpdb, s.pfam=s.pfam, s.sth=s.sth, s.pdb = s.pdb, atoms=atoms, pattern = pattern, subject = subject))
  cat("Doing pairwise minimal atom distance calculations.\n")
  pdm <- get_minimal_pairwise_dist_between_atoms_Ca_Cb(atoms[atomsidx])
  rownames(uniprotpdbidx) <- uniprotpdbidx$pdbidx
  pdm$alignmentidx.idx.i <- uniprotpdbidx[ as.character(pdm$pdb.resno.i), "alignmentidx"]
  pdm$alignmentidx.idx.j <- uniprotpdbidx[ as.character(pdm$pdb.resno.j), "alignmentidx"]
  pdm$uniprot.idx.i <- uniprotpdbidx[ as.character(pdm$pdb.resno.i), "seqidx"]
  pdm$uniprot.idx.j <- uniprotpdbidx[ as.character(pdm$pdb.resno.j), "seqidx"]

  return(list(pdm = pdm, uniprotpdbidx = uniprotpdbidx, aln = aln,
              annot = list(pdb = pdbEntry$annot, uniprot = uniprotseq)))

}


get_pdbs_for_uniprotID <- function(uniprotID = "FIBA_CHICK",
                                   pdbLocalPath = "data-raw/PDB/",
                                   resolution.th = 3.0, experimentalTech = "X-RAY DIFFRACTION") {
  mypdbs <- pdbUniprot[(pdbUniprot$Entry.name == uniprotID),]
  if(dim(mypdbs)[1] == 0) {
    stop(sprintf("No entries found for %s", uniprotID))
  }

  pdbAnnots <- NULL
  cntr <- 0
  for(pdbcode in unique(mypdbs$PDB)) {
    cntr <- cntr + 1
    fname <- get.pdb(as.character(sprintf("%s", pdbcode)),
                     path=pdbLocalPath, overwrite=F, gzip=T)
    if(fname == "1") {
      cat(sprintf("Could not download file for pdbcode: %s.\n", pdbcode))
      next
    }
    if(!file.exists(fname)) {
      cat(sprintf("Cannot find file: %s.\n", fname))
      next
    }
    cat(sprintf("Processing PDB file: %s (%d/%d).\n", fname, cntr, length(unique(mypdbs$PDB))))
    pdb <- read.pdb(fname)
    pdb.annot <- pdb.annotate(pdbcode)
    pdbAnnots <- rbind(pdbAnnots, pdb.annot)
  }
  if(length(pdbAnnots) == 0 | is.null(pdbAnnots)) {
    stop(sprintf("%s has 0 PDB annotations. Quiting.\n", uniprotID))
    return(NA)
  }
  pdbAnnots$structureId <- tolower(pdbAnnots$structureId)
  pdbs <- merge(mypdbs, pdbAnnots, by.x = c("PDB", "CHAIN"), by.y = c("structureId", "chainId"))
  if(as.numeric(pdbs$resolution) > resolution.th | pdbs$experimentalTechnique != experimentalTech) {
    cat("PDBs that will be excluded are:\n")
    print(pdbs[as.numeric(pdbs$resolution) > resolution.th | pdbs$experimentalTechnique != experimentalTech,
               c("PDB", "CHAIN", "resolution", "experimentalTechnique")])

  }
  pdbs <- pdbs[which(as.numeric(pdbs$resolution) <= resolution.th & pdbs$experimentalTechnique == experimentalTech), ]

  return(pdbs)
}


