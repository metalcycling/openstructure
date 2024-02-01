
e=io.LoadPDB("../demos/data/sdh.pdb")
scene.Add(gfx.Entity("x",gfx.CPK,e.Select("chain=A")))

# fast spheres with 3D sprites, default
scene["x"].GetOptions(gfx.CPK).SetSphereMode(1)

# turn on realtime ambient occlusion
scene.SetAmbientOcclusion(True)

# strength of the AO effect, default is 1.0
scene.SetAmbientOcclusionWeight(1.0)

# AO algorithms 0-3, default is 1
scene.SetAmbientOcclusionMode(1)

# orientation
scene.SetRTC(geom.Mat4(0.827929,-0.460893,-0.319547,29.688,
                       0.174973,0.753608,-0.633609,19.9405,
                       0.532839,0.468671,0.704578,73.0075,
                       -5.4704,1.54823,-116,1))
