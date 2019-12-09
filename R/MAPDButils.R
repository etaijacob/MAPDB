#MAPDButils.R

# Calculates lower triangle of i,j matrix:
get_minimal_pairwise_dist_between_atoms_Ca_Cb <- function(atoms) {

  pdm <- NULL
  nns <- combn(x = length(atoms), 2)
  cat(sprintf("Calculating %d distances.\n", length(nns)))
  get_my_dists <- function(mi) {
    #print(mi)
    i <- names(atoms)[mi[1]]
    j <- names(atoms)[mi[2]]

    anyAtom = min(as.matrix(pdist::pdist(atoms[[i]][, c("x", "y", "z")],
                                  atoms[[j]][, c("x", "y", "z")])))
    if(atoms[[i]]$resid[1] == "GLY") {
      cbi <- atoms[[i]][ atoms[[i]]$elety == "CA", c("x", "y", "z")]
    } else {
      cbi <- atoms[[i]][ atoms[[i]]$elety == "CB", c("x", "y", "z")]
    }
    if(atoms[[j]]$resid[1] == "GLY") {
      cbj <- atoms[[j]][ atoms[[j]]$elety == "CA", c("x", "y", "z")]
    } else {
      cbj <- atoms[[j]][ atoms[[j]]$elety == "CB", c("x", "y", "z")]
    }

    cai <- atoms[[i]][ atoms[[i]]$elety == "CA", c("x", "y", "z")]
    caj <- atoms[[j]][ atoms[[j]]$elety == "CA", c("x", "y", "z")]

    row <- data.frame(pdb.resno.i = atoms[[i]]$resno[1], pdb.resno.j = atoms[[j]]$resno[1],
                      aai = aaa2a[atoms[[i]]$resid[1], ], aaj = aaa2a[atoms[[j]]$resid[1], ],
                      dist.anyAtom = anyAtom,
                      dist.ca = min(as.matrix(pdist::pdist(cai, caj))),
                      dist.cb = min(as.matrix(pdist::pdist(cbi, cbj))),
                      stringsAsFactors=F)
    return(row)
  }
  pdm <- do.call(rbind, apply(nns, MARGIN = 2, get_my_dists))
  return(pdm)
}

read.fasta.alignment <- function(inputMSAfile, inputMSAfileCommentSep = "|") {

  msa <- seqinr::read.alignment(inputMSAfile, format = "fasta", forceToLower = F)
  if(!is.na(inputMSAfileCommentSep)) {
    mynames <- as.vector(sapply(msa$nam,
                                function(s) gsub(sprintf("^\\s+%s\\s+$", inputMSAfileCommentSep), "",
                                                 strsplit(x = s, sprintf("\\%s", inputMSAfileCommentSep))[[1]][1])))
    msa$nam <- mynames
  }
  return(msa)
}
read.stockholm.alignment <- function(file) {
  mm <- readLines(file)
  mm <- do.call(rbind, strsplit(mm[-grep("^#", mm)], "\\s+", perl = T))
  mm <- data.frame(name=mm[,1], seq=mm[,2], row.names = mm[,1], stringsAsFactors = F)
  return(mm)
}

get_pdb_atom_coordinates <- function(pdbcode, chain, pdbLocalPath = "data-raw/PDB/") {
  fname <- bio3d::get.pdb(as.character(sprintf("%s", pdbcode)),
                   path=pdbLocalPath, overwrite=F, gzip=T)
  if(!file.exists(fname)) {
    stop(sprintf("Cannot find file: %s.\n", fname))
  }

  print(fname)
  pdb <- bio3d::read.pdb(fname)
  pdb.annot <- bio3d::pdb.annotate(pdbcode)

  cbetas <- data.frame(pdb$atom, stringsAsFactors = F)
  cbetas$resno <- as.numeric(as.character(cbetas$resno))

  caidx <- bio3d::atom.select(pdb, string = "calpha", chain = chain)

  ca <- pdb$atom[ caidx$atom, ] #atoms[as.character(atoms$elety) == "CA", ]

  names(ca)[match(c("x", "y", "z"), names(ca))] <- c("PDB.ca_pos_x", "PDB.ca_pos_y", "PDB.ca_pos_z")
  cbidx <- bio3d::atom.select(pdb, elety = "CB", chain = chain)
  cb <- pdb$atom[ cbidx$atom, ] #cb <- atoms[as.character(atoms$elety) == "CB", ]
  cb <- cb[, c("resno", "x", "y", "z")]

  names(cb) <- c("resno", "PDB.cb_pos_x", "PDB.cb_pos_y", "PDB.cb_pos_z")
  cbetas <- merge(ca, cb, by = "resno", all.x  = T)
  cbetas$aatype <- aaa2a[ cbetas$resid,]
  atoms <- getAtomListByResno(pdb, chain)
  return(list(cbetas = cbetas, atoms = atoms, annot = pdb.annot))
}

getAtomListByResno <- function(pdb, chain) {
  myAtoms <- bio3d::atom.select(pdb, chain = chain)
  atomdf <- pdb$atom[myAtoms$atom, ]
  atomdf <- atomdf[ atomdf$type == "ATOM", ]
  myRess <- split(atomdf, atomdf$resno)
  return(myRess)
}

my3lettersAA <- c("STP", "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
                  "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR")

my1letterAA <- c("*", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")


aaa2a <- data.frame(my1letterAA, row.names = my3lettersAA, stringsAsFactors = F)
a2aaa <- data.frame(my3lettersAA, row.names = my1letterAA, stringsAsFactors = F)

