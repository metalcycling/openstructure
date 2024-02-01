struct=io.LoadEntity("sdh.pdb")
surf=io.LoadSurface("sdh","msms")
surf.Attach(struct,5.0)
go=gfx.Surface("surf",surf)
# this is broken
#go.clip=True
#go.clip_plane=geom.Vec4(0,0,-1,-130)
# this works, z offset to near clipping plane, must be >0
go.clip_offset=1
go.solid=True
go.solid_color=gfx.RGB(1,1,0)

scene.rtc=geom.Mat4(0.914075, -0.355397, 0.195293, 29.688, 0.405471, 0.7934, -0.45397, 19.9405, 0.00639453, 0.494153, 0.869337, 73.0075, 0, 0, -143.077, 1)

scene.Add(go)
scene.SetCenter(go.GetCenter())
scene.Autoslab()
scene.auto_autoslab=False
