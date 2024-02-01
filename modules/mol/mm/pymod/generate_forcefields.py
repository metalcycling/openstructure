import os
from ost.mol import mm

path_to_gromacs = "xyz"
charmm_ff = "charmm27.ff"
amber_ff = "amber03.ff"

charmm_path = os.path.join(path_to_gromacs, "share", "top", charmm_ff)
amber_path = os.path.join(path_to_gromacs, "share", "top", amber_ff)

# Let's first do CHARMM
reader = mm.FFReader(charmm_path)
reader.ReadGromacsForcefield()
reader.ReadResidueDatabase("aminoacids")
reader.ReadResidueDatabase("rna")
reader.ReadResidueDatabase("dna")
reader.ReadResidueDatabase("lipids")
reader.ReadITP("tip3p")
reader.ReadITP("ions")
ff = reader.GetForcefield()
ff.Save("CHARMM27.dat")

# Let's do amber
reader = mm.FFReader(amber_path)
reader.ReadGromacsForcefield()
reader.ReadResidueDatabase("aminoacids")
reader.ReadResidueDatabase("rna")
reader.ReadResidueDatabase("dna")
reader.ReadITP("tip3p")
reader.ReadITP("ions")
ff = reader.GetForcefield()
ff.Save("AMBER03.dat")

