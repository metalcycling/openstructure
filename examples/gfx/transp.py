e=io.LoadPDB("../demos/data/sh2.pdb")
s=io.LoadSurface("../demos/data/sh2","msms")

scene.Add(gfx.Entity("x",gfx.HSC,e.Select("chain=A")))
scene["x"].SetColor(gfx.Color(0.9,0.95,1.0))
scene["x"].SetColor(gfx.Color(1.0,0.9,0.7),"rtype=helix")
scene["x"].SetColor(gfx.Color(0.5,1.0,0.7),"rtype=strand")

scene.Add(gfx.Surface("s",s))
scene["s"].SetOpacity(0.4)
scene["s"].SetMat(0,1,0,0)

scene.bg=gfx.WHITE
