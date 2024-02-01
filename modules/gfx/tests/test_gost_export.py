import desost_util as util

e=io.LoadPDB("../../../examples/demos/data/sdh.pdb")

s=util.CreateSurf(e.Select("chain=A"),density=4)

scene.Add(gfx.Entity("trace",gfx.TUBE,e.Select("chain=A")))
scene["trace"].tube_options.arc_detail=8
scene["trace"].tube_options.spline_detail=8
scene["trace"].tube_options.tube_radius=1.0
grad=gfx.Gradient()
grad.SetColorAt(0.0,gfx.RED)
grad.SetColorAt(0.5,gfx.GREEN)
grad.SetColorAt(1.0,gfx.BLUE)
scene["trace"].ColorBy("rnum",grad)

scene.Add(gfx.Surface("surf",s))
scene["surf"].ColorBy("rnum",grad)

if not gui_mode:
  ge = gfx.GostExporter("test.gost")
  scene.Export(ge)
