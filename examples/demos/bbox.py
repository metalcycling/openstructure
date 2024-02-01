# This code example shows how to calculate and display an oriented bounding box 
# that contains all the atoms of an entity.

# remove all objects from scene, just in case
scene.RemoveAll()

# Load sdh molecule
sdh=io.LoadPDB('data/sdh.pdb')
# Select backbone atoms of residues 99-128, which form a helix.
helix=sdh.Select('rnum=99:128 and cname=A and aname=CA,C,N,O')
# create a graphical representation of the helix and render it with simple 
# lines.
go=gfx.Entity('helix', gfx.SIMPLE, helix)
# add it to the scene
scene.Add(go)

# calculate the bounding box for the helix
bbox=mol.BoundingBoxFromEntity(helix)

# create a graphical representation of the bounding box and add it to the scene
bb=gfx.Cuboid('xxx', bbox)
bb.SetFillColor(gfx.Color(0.5, 0.8, 0.5, 0.2))
scene.Add(bb)

# center the camera on the graphical object
scene.center=go.center
