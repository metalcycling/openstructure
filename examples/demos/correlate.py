# this script shows how to score a set of candidate fragments according to their 
# correlation with a density map. fragments with a higher correlation are ranked 
# better than fragments with a low correlation.
#

# tabula rasa first
scene.RemoveAll()

import functools
import ost.img.alg as img_alg
from ost.mol.alg import EntityToDensityRosetta, HIGH_RESOLUTION
sdh=io.LoadPDB('data/sdh-complete.pdb')
sdh_within=sdh.Select('rnum=94:101 and aname=CA,C,N,O')

padding=4
# we use a isotropic sampling of 1Anstroem
bounds=sdh_within.bounds
actual_size=bounds.max-bounds.min
size=img.Size(int(actual_size.x+2*padding), 
              int(actual_size.y+2*padding), 
              int(actual_size.x+2*padding))

def CreateDensityMap(size, origin):
  # create 3D image to hold the density information. Pixel sampling is 1A
  # by default. The image is empty
  dmap=img.CreateImage(size)
  dmap.SetAbsoluteOrigin(origin)  
  return dmap
origin=bounds.min-geom.Vec3(padding, padding, padding)
dmap=CreateDensityMap(size, origin)
cmap=CreateDensityMap(size, origin)

# create density for residue 94 to residue 101
EntityToDensityRosetta(sdh_within, dmap, HIGH_RESOLUTION, 5.0)
map_go=gfx.MapIso('Density', dmap, 0.0)
map_go.SetLevel(map_go.GetMean()*2.0)
scene.Add(map_go)
scene.CenterOn(map_go)

# load all loop candidates stored in the PDB file. note the use of the 
# load_multi parameter
candidates=io.LoadPDB('data/loop-candidates.pdb', 
                      load_multi=True)

print(len(candidates))
#-------------------------------------------------------------------------------
def Correlate(candidates):
  # this function converts the candidate loop into a density and calculates the 
  # real-spatial cross correlation with the actual density. The correlation 
  # coefficient is stored as the "correl" property
  
  for index, candidate in enumerate(candidates):
    EntityToDensityRosetta(candidate.CreateFullView(), cmap, 
                           HIGH_RESOLUTION, 5.0, True)
  
    correl=img_alg.RealSpatialCrossCorrelation(dmap, cmap, dmap.GetExtent())
    candidate.SetFloatProp('correl', correl)
#-------------------------------------------------------------------------------  
  
  
def Visualize(candidates):  
  # finally, sort the candidates by correlation and display them. 
  min_correl=min([c.GetFloatProp('correl') for c in candidates])
  max_correl=max([c.GetFloatProp('correl') for c in candidates])

  gradient=gfx.Gradient()
  gradient.SetColorAt(0.0, gfx.RED)
  gradient.SetColorAt(1.0, gfx.GREEN)

  # sort candidates by correlation
  def cmp_correl(lhs, rhs):
    lhs_correl=lhs.GetFloatProp('correl')
    rhs_correl=rhs.GetFloatProp('correl')
    if lhs_correl>rhs_correl:
      return -1
    elif lhs_correl==rhs_correl:
      return 0
    else:
      return 1
  candidates=sorted(candidates, key=functools.cmp_to_key(cmp_correl))

  # create render objects and add to scene
  for index, candidate in enumerate(candidates):
    candidate_go=gfx.Entity('candidate %03d' % index, candidate)
    norm=candidate.GetFloatProp('correl')/(max_correl-min_correl)
    candidate_go.SetColor(gradient.GetColorAt(norm))
    scene.Add(candidate_go)


Correlate(candidates)
Visualize(candidates)
