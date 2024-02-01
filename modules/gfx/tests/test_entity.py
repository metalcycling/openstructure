eh=io.LoadPDB("1crn",remote=True)

scene.Add(gfx.Entity("simple_default",gfx.SIMPLE,eh))
scene.Add(gfx.Entity("simple_mod",gfx.SIMPLE,eh))
scene["simple_mod"].ColorBy("abfac",gfx.Gradient([gfx.BLUE,gfx.RED]))
scene["simple_mod"].simple_options.line_width=3.5
scene["simple_mod"].simple_options.aa_lines=True

scene.Add(gfx.Entity("custom_default",gfx.CUSTOM,eh))
scene.Add(gfx.Entity("custom_mod",gfx.CUSTOM,eh))
scene["custom_mod"].custom_options.bond_rad=0.45
scene["custom_mod"].custom_options.sphere_rad=0.55
scene["custom_mod"].SetColor(gfx.ORANGE,"ele=C")

scene.Add(gfx.Entity("cpk_default", gfx.CPK,eh))
scene["cpk_default"].cpk_options.sphere_detail=2 # should have no influence
scene.Add(gfx.Entity("cpk_mod", gfx.CPK,eh))
scene["cpk_mod"].cpk_options.sphere_mode=0
scene["cpk_mod"].cpk_options.sphere_detail=2

scene.Add(gfx.Entity("ltrace_default", gfx.LINE_TRACE,eh))
scene.Add(gfx.Entity("ltrace_mod", gfx.LINE_TRACE,eh))
scene["ltrace_mod"].ColorBy("rnum",gfx.Gradient([gfx.PURPLE,gfx.YELLOW,gfx.CYAN]))
scene["ltrace_mod"].line_trace_options.aa_lines=True
scene["ltrace_mod"].line_trace_options.line_width=3.5

scene.Add(gfx.Entity("trace_default", gfx.TRACE,eh))
scene.Add(gfx.Entity("trace_mod", gfx.TRACE,eh))
scene["trace_mod"].ColorBy("rnum",gfx.Gradient([gfx.RED,gfx.GREEN,gfx.BLUE]))
scene["trace_mod"].trace_options.arc_detail=6
scene["trace_mod"].trace_options.tube_radius=1.3

scene.Add(gfx.Entity("tube_default", gfx.TUBE,eh))
eh2=eh.Copy()
scene.Add(gfx.Entity("tube_mod", gfx.TUBE,eh2))
bfmin=1e-6
bfmax=-1e-6
for a in eh2.Select("peptide=1 and aname=CA").atoms:
  bfmin=min(bfmin,a.b_factor)
  bfmax=max(bfmax,a.b_factor)
for a in eh2.Select("peptide=1 and aname=CA").atoms:
  a.handle.SetFloatProp("trace_rad",5.0*min(1.0,max(0.0,(a.b_factor-bfmin)/(bfmax-bfmin)))+0.3)
scene["tube_mod"].ColorBy("rbfac",gfx.Gradient([gfx.BLUE,gfx.RED]))
scene["tube_mod"].tube_options.spline_detail=8
scene["tube_mod"].tube_options.arc_detail=8
scene["trace_mod"].tube_options.tube_radius=1.0

scene.Add(gfx.Entity("hsc_default", gfx.HSC,eh))
scene.Add(gfx.Entity("hsc_mod", gfx.HSC,eh))
scene["hsc_mod"].SetColor(gfx.PURPLE,"rtype=H")
scene["hsc_mod"].SetColor(gfx.ORANGE,"rtype=E")
scene["hsc_mod"].cartoon_options.helix_width=1.3
scene["hsc_mod"].cartoon_options.helix_thickness=0.6
scene["hsc_mod"].cartoon_options.strand_width=1.6
scene["hsc_mod"].cartoon_options.strand_thickness=0.4
scene["hsc_mod"].cartoon_options.tube_radius=0.5



scene.Match("*").Hide()

