s=io.LoadSurface("../demos/data/sh2","msms")
scene.Add(gfx.Surface("s",s))

# alternative shading mode
scene.SetShadingMode("hf")

# add outlines to surface
scene["s"].SetOutline(True)
scene["s"].SetOutlineMode(3)

# turn on realtime ambient occlusion
scene.SetAmbientOcclusion(True)

# strength of the AO effect, default is 1.0
scene.SetAmbientOcclusionWeight(1.0)

# AO algorithms 0-3, default is 1
scene.SetAmbientOcclusionMode(1)

