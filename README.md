# MAPDB

An R package for mapping Pfam multiple sequence alignment (MSA in sth format) 
or Uniprot (MSA in fasta format) to PDB residue coordinates.


Using this package requires an interent connection. If more then several examples are required,
you should consider downloading first PDB and STH databases.

This package can be used, for example, to calculate the accuracy of contact predictions based on MSA data (co-evolution).

Please see an example of an application @ Codon-level information improves predictions of inter-residue contacts in proteins by correlated mutation analysis - See more at: http://elifesciences.org/content/early/2015/09/14/eLife.08932#.dpuf

If you use this package please cite:

Etai Jacob, Ron Unger and Amnon Horovitz (2015). Codon-level information improves predictions of inter-residue contacts in proteins by correlated mutation analysis. eLife 2015;10.7554/eLife.08932 http://dx.doi.org/10.7554/eLife.08932.

# To install:

devtools::install_github("etaijacob/MAPDB")

# Example:

In this example we download all relevant files for PF00075,
including PDB files and sth alignment for the pfam domain.
We limit the download to a total of 12 pdb instances with a resolution below 3.0 Ang.
If any of the pdb files was already downloaded and located in the pdb.local.db,
then it will be used to avoid unnecessary download.
```
myPfamMap <- PDBs2Pfam(pfamid = "PF00075", sth.local.db = "/tmp/",
                      max2take = 2, resolution.th = 3.0)
``` 
In this example we have already downloaded the sth file
therefore including it in the call for the function in order
to avoid unnecessary download

```
myPfamMap <- PDBs2Pfam(pfamid = "PF00075", pfamSthFile = "/tmp/PF00075_rp75.sth",
                       sth.local.db = "/tmp/", max2take = 2, resolution.th = 3.0)
```


