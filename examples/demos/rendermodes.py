# This script shows how to switch between different render modes 
# programmatically and change render options.


# remove all objects from scene, just in case
scene.RemoveAll()

# load pdb structure
eh=io.LoadPDB("data/sdh.pdb")


sdh_go=gfx.Entity("SDH2", eh.Select("cname=A"))
# create graphical object and add it to the scene
scene.Add(sdh_go)
scene.CenterOn(sdh_go)

def Tube():
  sdh_go.SetRenderMode(gfx.TUBE)
  sdh_go.tube_options.SetTubeRadius(0.8)
  # apply color gradient for the atomic bfactors
  sdh_go.ColorBy("abfac",gfx.BLUE,gfx.RED)
  # and apply radius gradient as well 
  sdh_go.RadiusBy("abfac",0.8,2.2)
  
def Trace():
  sdh_go.SetRenderMode(gfx.LINE_TRACE)
  sdh_go.GetOptions(gfx.LINE_TRACE).SetLineWidth(2.5)

def Spheres():
  sdh_go.SetRenderMode(gfx.CPK)
  sdh_go.GetOptions(gfx.CPK).SetSphereDetail(4)

def Cartoon():  
  sdh_go.SetRenderMode(gfx.RenderMode.HSC)
  sdh_go.SetColor(gfx.Color(0.5, 0.5, 0.5), '')
  sdh_go.SetColor(gfx.Color(0,1,0),"rtype=H")
  sdh_go.SetDetailColor(gfx.Color(0.7,1.0,0.8),"rtype=H")
  # strands
  sdh_go.SetColor(gfx.Color(1,0,0),"rtype=E")
  sdh_go.SetDetailColor(gfx.Color(1.0,0.8,0.2),"rtype=E")
  # these are the default params for the rendering
  #sdh_go.cartoon_options.SetArcDetail(4)    # circular profile detail
  #sdh_go.cartoon_options.SetSphereDetail(4) # sphere detail
  #sdh_go.cartoon_options.SetSplineDetail(8) # segments between backbone atoms
  #sdh_go.cartoon_options.SetTubeRadius(0.4) # coil radius 
  #sdh_go.cartoon_options.SetTubeRatio(1.0)  # coil axial ratio
  #sdh_go.cartoon_options.SetHelixWidth(1.1)   # helical width
  #sdh_go.cartoon_options.SetHelixThickness(0.2)   # helical thickness
  #sdh_go.cartoon_options.SetStrandWidth(1.2)  # strand width
  #sdh_go.cartoon_options.SetStrandThickness(0.2)  # strand thickness

Cartoon()

