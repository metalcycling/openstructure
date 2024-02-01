struct=io.LoadEntity("sdh.pdb")
surf=io.LoadSurface("sdh","msms")
surf.Attach(struct,5.0)
go=gfx.Surface("surf",surf)
scene.Add(go)
scene.SetCenter(go.GetCenter())
scene.Autoslab()

# color surface by atomic properties

# slow version, via positional search
#gradient=gfx.Gradient()
#gradient.SetColorAt(0.0, gfx.Color(0.0, 0.0, 1.0, 1.0))  # blue
#gradient.SetColorAt(0.5, gfx.Color(1.0, 1.0, 1.0, 1.0))  # white
#gradient.SetColorAt(1.0, gfx.Color(1.0, 0.0, 0.0, 1.0))  # red
#ev=struct.Select("")
#go.ColorBy(ev,"abfac",gradient,30,80)

# fast version, via attached atoms
go.ColorBy("abfac",gfx.Color(0,0,1),gfx.Color(1,0,0),30,80)

# clear
go.SetColor(gfx.Color(1,1,1))

# color by selection, again based on the attached structure
go.SetColor(gfx.Color(1,1,0),"rname=PHE,TYR,TRP")
