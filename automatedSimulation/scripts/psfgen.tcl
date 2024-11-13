package require psfgen

topology top_all27_prot_lipid.inp

pdbalias residue HIS HSE ;# we define some aliases histidine is called HIS in p>
pdbalias atom ILE CD1 CD

segment U {pdb protein.pdb} ;# We generate a segment called U using atoms form >
coordpdb protein.pdb U

guesscoord ;# we guess the coordines of missing atoms (mainly H, but also some >

writepdb protein/system.pdb
writepsf protein/system.psf
