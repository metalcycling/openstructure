from ost.mol import alg
# tabula rasa first
scene.RemoveAll()

# load the fragment
frag=io.LoadPDB('data/1crn.pdb')

# now let's calculate the size of the map by taking the bounding box of
# the fragment and adding 5 Angstrom on each side.
frag_bounds=frag.bounds
real_size=frag_bounds.max-frag_bounds.min+geom.Vec3(10,10,10)
# We use 2 pixels per Angstroem, so we have to multiply the size by 2.
map_size=img.Size(int(real_size.x)*2, int(real_size.y)*2, int(real_size.z)*2)
den_map=img.CreateImage(map_size)
# we use 2 pixels per Angstrom
den_map.SetSpatialSampling(2)
den_map.SetAbsoluteOrigin(frag_bounds.min-geom.Vec3(5, 5, 5))


# Calculate a density map from the fragment using a simple real-space method:
#
# This function is based on Dimaio et al.,J. Mol. Biol(2009)
# In essence every atom gives rise to a Gaussian centered at the atom's
# position. Sigma and Magnitude of the Gaussian are calculated from the
# atomic mass.
alg.EntityToDensityRosetta(frag.CreateFullView(), den_map,
                           alg.HIGH_RESOLUTION, 5.0)

# Create Graphical Representations and add to the Scene
gfx_frag=gfx.Entity('Fragment', gfx.HSC, frag)
gfx_frag.SetColor(gfx.YELLOW)
gfx_map=gfx.MapIso('Density', den_map, 250.0)
scene.Add(gfx_frag)
scene.Add(gfx_map)
scene.CenterOn(gfx_frag)
