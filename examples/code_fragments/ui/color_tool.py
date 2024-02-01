"""
This script illustrates how to implement a custom tool for interacting with the
3D scene.
"""
from ost import gui


class ColorTool(gui.Tool):
  def __init__(self):
    gui.Tool.__init__(self, "Color Tool")   
    self.tool_options.AddOption(gui.ToolOptionFloat("r", "Red", 0.5, 0.0, 1.0))
    self.tool_options.AddOption(gui.ToolOptionFloat("g", "Green", 0.5, 0.0, 1.0))
    self.tool_options.AddOption(gui.ToolOptionFloat("b", "Blue", 0.5, 0.0, 1.0))
    apply_to=gui.ToolOptionEnum("apply_to", "Apply To")
    apply_to.Add("Atom", 0)
    apply_to.Add("Residue", 1)
    apply_to.Add("Chain", 2)
    apply_to.SetIndex(0)
    self.tool_options.AddOption(apply_to)

  def GetColor(self):
    return gfx.Color(self.tool_options.GetOption('r').value,
                     self.tool_options.GetOption('g').value,
                     self.tool_options.GetOption('b').value)
                     
  def CanOperateOn(self, obj):
    return True
  def GetIconPath(self):
    return ''
  def GetAtomsToColor(self, atom):
    apply_to=self.tool_options.GetOption('apply_to').value
    if apply_to==0:
      return [atom]
    if apply_to==1:
      return atom.residue.atoms
    if apply_to==2:
      l=[]
      for r in atom.residue.chain.residues:
        for a in r.atoms:
          l.append(a)
      return l

  def Click(self, event):
    obj, atom=gfx.PickAtom(scene, event.pos.x(), event.pos.y())
    if atom.IsValid():
      atoms=self.GetAtomsToColor(atom)
      color=self.GetColor()
      for atom in atoms:
        obj.SetColorForAtom(color, atom.handle)

color_tool=ColorTool()

gui.ToolManager.Instance().AddTool(color_tool)

e=io.LoadEntity(len(sys.argv)>1 and sys.argv[1] or '../entity/fragment.pdb')
g=gfx.Entity('e', e)
g.SetRenderMode(gfx.RenderMode.CUSTOM)
scene.Add(g)
scene.center=g.center
