# remove all objects from scene, just in case
scene.RemoveAll()
# Load chain A of SDH
ent = io.LoadPDB('data/sdh.pdb', restrict_chains='A')

# create graphics object
go = gfx.Entity('SDH', ent)
# set rendermode of full graphical object to interpolated spline
go.SetRenderMode(gfx.SLINE)
# select sidechains (including bond to C-alpha and N in case of proline)
sidechains = ent.Select('aname!=CA,C,N,O', mol.EXCLUSIVE_BONDS)
# set rendering of sidechains to simple
go.SetRenderMode(gfx.SIMPLE, sidechains)

# add object to scene
scene.Add(go)

# center the scene on the geometric center of the SDH
scene.CenterOn(go)

print('Demo 1: loading and displaying a pdb file')
