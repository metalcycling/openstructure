s=io.LoadSurface("../demos/data/sh2","msms")
scene.Add(gfx.Surface("s",s))

# alternative shading mode
scene.SetShadingMode("hf")

# add outlines to surface
scene["s"].SetOutline(True)
scene["s"].SetOutlineMode(3)

# turn on solid rendering
scene["s"].solid=True
scene["s"].solid_color=gfx.RGB(0.8,0.8,0.7)

scene.AutoAutoslab(False)
scene.far=170
scene.near=100
scene.fog=False
