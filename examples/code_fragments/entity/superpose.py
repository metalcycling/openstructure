import math
import ost.mol.alg as a

#loading the two entities to be superposed
one=io.LoadEntity('pdb1ake.ent')
two=io.LoadEntity('pdb4ake.ent')

#generate a set of exactly corresponding entity views
ev_a=one.Select('cname=A and (rnum=1:100) and aname=C,CA,N,O')
ev_b=two.Select('cname=A and (rnum=1:100) and aname=C,CA,N,O')

#superposing the first entity view onto the second
#the transformation is directly applied to ev_a
res=a.SuperposeSVD(ev_a, ev_b)
#the superposition result res contains the rmsd between the two entity views
# and additionally more information in case IterativeSuperposition was called

#add the two entities to the scene, color them individually and center
g=gfx.Entity("second", ev_b)
scene.Add(g)
h=gfx.Entity("first",ev_a)
scene.Add(h)
h.SetColor(gfx.RED)
g.SetColor(gfx.BLUE)
scene.CenterOn(g)

