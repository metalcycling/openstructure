pdb_name='sdh.pdb'
if len(sys.argv)>1:
  pdb_name=sys.argv[1]
ent=io.LoadEntity(pdb_name)

vv=ent.Select('peptide=true and aname=CA,N,C,O')
gradient=gfx.Gradient()
gradient.SetColorAt(1.0, gfx.Color(1.0, 0.0, 0.0, 1.0))
gradient.SetColorAt(0.5, gfx.Color(0.0, 1.0, 0.0, 1.0))
gradient.SetColorAt(0.0, gfx.Color(0.0, 0.0, 1.0, 1.0))

go=gfx.Entity(os.path.splitext(os.path.basename(pdb_name))[0], vv)
go.ColorBy('abfac', gradient)
scene.Add(go)
scene.center=go.center