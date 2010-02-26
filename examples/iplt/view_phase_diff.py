import sys
import math
import ost.iplt.alg

image1=io.LoadImage(sys.argv[1])
image2=io.LoadImage(sys.argv[2])
if image1.GetExtent() != image2.GetExtent():
  raise RuntimeError('The input images should have the same size.')
image1.CenterSpatialOrigin()
image2.CenterSpatialOrigin()
image1.ApplyIP(ost.iplt.alg.DFT())
image2.ApplyIP(ost.iplt.alg.DFT())
ex_it=iplt.ExtentIterator(image1.GetExtent())
diff_image=iplt.CreateImage(image1.GetExtent())
for pixel in ex_it:
  phase1=iplt.Phase(image1.GetComplex(pixel))
  phase2=iplt.Phase(image2.GetComplex(pixel))
  phase_diff=phase1-phase2
  diff_image.SetReal(pixel,180.0*float(phase_diff)/math.pi)
v=gui.CreateDataViewer(diff_image)
v.SetName("Phase difference (in degrees)")
