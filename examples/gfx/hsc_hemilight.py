e=io.LoadPDB("../demos/data/sdh.pdb")
scene.Add(gfx.Entity("x",gfx.HSC,e.Select("chain=A")))
scene["x"].GetOptions(gfx.HSC).SetTubeRadius(0.6)
# spiral helices (default)
scene["x"].GetOptions(gfx.HSC).SetHelixMode(0)
# straight strands
scene["x"].GetOptions(gfx.HSC).SetStrandMode(1)
# colors
scene["x"].SetColor(gfx.Color(0.9,0.95,1.0))
scene["x"].SetColor(gfx.Color(1.0,0.9,0.7),"rtype=helix")
scene["x"].SetDetailColor(gfx.Color(1.0,0.95,0.78),"rtype=helix")
scene["x"].SetColor(gfx.Color(0.5,1.0,0.7),"rtype=strand")
scene["x"].SetDetailColor(gfx.Color(0.7,1.0,0.79),"rtype=strand")

# switch to hemilighting mode
scene.SetShadingMode("hf")
# add outline to improved visibility
scene["x"].Outline(True)
scene["x"].SetOutlineMode(1)
scene["x"].SetOutlineWidth(4.5)

# looks best on white bg
scene.bg=gfx.WHITE

scene.SetRTC(geom.Mat4(-0.589889,-0.594876,0.546034,28.2245,
                        0.806165,-0.472502,0.356147,20.1685,
                        0.0461388,0.65028,0.758291,73.5665,
                        -3.52709,1.72253,-120,1))


