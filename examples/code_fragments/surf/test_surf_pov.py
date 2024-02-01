# run with ost_cl
struct=io.LoadEntity("sdh.pdb")
surf=io.LoadSurface("sdh","msms")
surf.Attach(struct,5.0)
go=gfx.Surface("surf",surf)
scene.Add(go)
scene.SetCenter(go.GetCenter())
scene.Autoslab()
scene.SetFog(True)
scene.Autoslab(True)
scene.fog_near_offset=10
scene.SetBackground(gfx.WHITE)

# fast version, via attached atoms
go.ColorBy("abfac",gfx.Color(0,0,1),gfx.Color(1,0,0),30,80)

scene.SetRTC(geom.Mat4(0.656168,-0.753826,-0.0341944,29.688,0.471087,0.444621,-0.761822,19.9405,0.589488,0.483781,0.646863,73.0075,-1.51607,3.36172,-100,1))

scene.ExportPov("test")

