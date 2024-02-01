s=io.LoadSurface("../demos/data/sh2","msms")
scene.Add(gfx.Surface("s",s))

# alternative shading mode
scene.SetShadingMode("hf")

# add outlines to surface
scene["s"].Outline(True)
scene["s"].SetOutlineMode(3)

# turn on realtime shadow mapping
scene.SetShadow(True)

scene.SetShadowQuality(3) # default

