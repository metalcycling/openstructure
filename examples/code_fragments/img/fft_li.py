import sys
import ost.img.alg
if len(sys.argv)<2:
  images=io.LoadImageList(['square.png','circle.png'])  
else:
  images=io.LoadImageList(sys.argv[1:])
viewers=[]
i=1
main_area=gui.GostyApp.Instance().perspective.main_area
for im in images:
   im.CenterSpatialOrigin()
   im.ApplyIP(img.alg.DFT()) 
   v=gui.CreateDataViewer(im,"Image")
   viewers.append(v)
   i+=1