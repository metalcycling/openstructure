import sys
import math
import ost.img.alg

if len(sys.argv)==2:
  image1=io.LoadImage(sys.argv[1])
  image2=io.LoadImage(sys.argv[2])
else:
  image1=io.LoadImage('square.png')
  image2=io.LoadImage('circle.png')
if image1.GetExtent() != image2.GetExtent():
  raise RuntimeError('The input images should have the same size.')
image1.CenterSpatialOrigin()
image2.CenterSpatialOrigin()
image1.ApplyIP(ost.img.alg.DFT())
image2.ApplyIP(ost.img.alg.DFT())
ex_it=img.ExtentIterator(image1.GetExtent())
diff_image=img.CreateImage(image1.GetExtent())
for pixel in ex_it:
  phase1=img.Phase(image1.GetComplex(pixel))
  phase2=img.Phase(image2.GetComplex(pixel))
  phase_diff=phase1-phase2
  diff_image.SetReal(pixel,180.0*float(phase_diff)/math.pi)
v=gui.CreateDataViewer(diff_image,"Phase difference (in degrees)")
