"""
Functions to calculate helix kinks: bend, face shift and wobble angles

Author: Niklaus Johner
"""

import os
import ost

def __FindProline(sele,proline):
  if not sele.IsValid():
    print('selection is not valid')
    raise RuntimeError
  if proline==False:
    proline=sele.Select('rname=PRO')
  if not proline.GetResidueCount()==1:
    print(proline.GetResidueCount(),'prolines in the selection. One proline is needed')
    raise RuntimeError
  proline=proline.residues[0]
  proline_ca=proline.FindAtom('CA')
  if not proline_ca.IsValid():
    print('proline has no CA atom')
    raise RuntimeError
  return (proline,proline_ca)
  proline.GetNumber().num
  
def __SelectPreAndPostProline(sele,proline_num):
  pre_proline=sele.Select('rnum<'+str(proline_num))
  post_proline=sele.Select('rnum>'+str(proline_num))
  print('pre-proline residues')
  for res in pre_proline.residues:
    print(res)
  print('post-proline residues')
  for res in post_proline.residues:
    print(res)
  if pre_proline.GetResidueCount()<4 or post_proline.GetResidueCount()<4:
    print('pre and post proline helices should be at least 4 residues long, 7 for better stability')
    raise RuntimeError
  return (pre_proline,post_proline)

def __FindCa3AndCa4(sele,proline_ca,proline_num):
  ca_3=sele.FindAtom(proline_ca.GetHandle().GetChain().GetName(),proline_num-3,'CA')
  ca_4=sele.FindAtom(proline_ca.GetHandle().GetChain().GetName(),proline_num-4,'CA')
  if not (ca_3.IsValid() and ca_4.IsValid()):
    print('CA not found in (i-4) or (i-3) residue')
    raise RuntimeError
  return (ca_3,ca_4)


def __CalculateBendAngle(pre_proline_axis,post_proline_axis):
  return ost.geom.Angle(pre_proline_axis,post_proline_axis)
  
def __CalculateWobbleAngle(pre_proline_axis,post_proline_axis,post_proline_centers,proline_pos):
  p1=proline_pos-post_proline_centers
  n1=p1-ost.geom.Dot(p1,post_proline_axis)*post_proline_axis
  p2=-pre_proline_axis
  n2=p2-ost.geom.Dot(p2,post_proline_axis)*post_proline_axis
  sign=ost.geom.Dot(ost.geom.Cross(pre_proline_axis,n2),n2-n1)
  sign=sign/abs(sign)
  return sign*ost.geom.Angle(n1,n2)

def __CalculateFaceShift(pre_proline_axis,post_proline_axis,pre_proline_centers,post_proline_centers,proline_pos,ca3_pos,ca4_pos,bend_angle):
  p1=proline_pos-post_proline_centers
  n1=p1-ost.geom.Dot(p1,post_proline_axis)*post_proline_axis
  p2=(ca3_pos+ca4_pos)/2.0-pre_proline_centers
  n2=p2-ost.geom.Dot(p2,pre_proline_axis)*pre_proline_axis
  #Here we want to apply a rotation to n1, to align the post-proline axis on the pre-proline axis
  R=ost.geom.AxisRotation(ost.geom.Cross(post_proline_axis,pre_proline_axis),bend_angle)
  n1=R*n1
  #We also need to determine the sign of the angle
  sign=ost.geom.Dot(ost.geom.Cross(pre_proline_axis,n2),n2-n1)
  sign=sign/abs(sign)
  return sign*ost.geom.Angle(n1,n2)


def AnalyzeHelixKink(t, sele, proline=False):
  """
  This function calculates the bend, wobble and face-shift angles in an alpha-
  helix over a trajectory. The determination is more stable if there are at
  least 4 residues on each side (8 is even better) of the proline around which
  the helix is kinked. The selection should contain all residues in the correct
  order and with no gaps and no missing C-alphas.

  :param t: The trajectory to be analyzed
  :type t: :class:`~ost.mol.CoordGroup`
  :param sele: A selection containing the alpha helix to be analyzed
  :type sele: :class:`~ost.mol.EntityView`
  :param proline: A selection containing only the proline (or another residue)
                  around which the helix is kinked. If False, the proline will
                  be searched for automatically
  :type proline: :class:`ost.mol.EntityView`

  :return: A tuple (bend_angle, face_shift, wobble_angle).
  :rtype: (FloatList, FLoatList, FloatList)
  """
  n_frames=t.GetFrameCount()
  (proline,proline_ca)=__FindProline(sele,proline)
  proline_num=proline.GetNumber().num
  (pre_proline,post_proline)=__SelectPreAndPostProline(sele,proline_num)
  (ca_3,ca_4)=__FindCa3AndCa4(sele,proline_ca,proline_num)
  #Here we extract the necessary information from the trajectory
  pre_proline_axis=ost.geom.Vec3List()
  post_proline_axis=ost.geom.Vec3List()
  pre_proline_centers=ost.geom.Vec3List()
  post_proline_centers=ost.geom.Vec3List()
  ost.mol.alg.AnalyzeAlphaHelixAxis(t,pre_proline,pre_proline_axis,pre_proline_centers)
  ost.mol.alg.AnalyzeAlphaHelixAxis(t,post_proline,post_proline_axis,post_proline_centers)
  proline_pos=ost.mol.alg.AnalyzeAtomPos(t,proline_ca.GetHandle())
  ca3_pos=ost.mol.alg.AnalyzeAtomPos(t,ca_3.GetHandle())
  ca4_pos=ost.mol.alg.AnalyzeAtomPos(t,ca_4.GetHandle())
  #Now we calculate the bend angle
  bend_angle=ost.FloatList()
  face_shift=ost.FloatList()
  wobble_angle=ost.FloatList()
  for i in range(n_frames):
    bend_angle.append(__CalculateBendAngle(pre_proline_axis[i],post_proline_axis[i]))
    face_shift.append(__CalculateFaceShift(pre_proline_axis[i],post_proline_axis[i],pre_proline_centers[i],post_proline_centers[i],proline_pos[i],ca3_pos[i],ca4_pos[i],bend_angle[i]))
    wobble_angle.append(__CalculateWobbleAngle(pre_proline_axis[i],post_proline_axis[i],post_proline_centers[i],proline_pos[i]))
  return (bend_angle,face_shift,wobble_angle)


def CalculateHelixKink(sele, proline=False):
  """
  This function calculates the bend, wobble and face-shift angles in an alpha-
  helix of an EntityView. The determination is more stable if there are at least
  4 residues on each side (8 is even better) of the proline around which the
  helix is kinked. The selection should contain all residues in the correct
  order and with no gaps and no missing C-alphas.
  
  :param sele: A selection containing the alpha helix to be analyzed
  :type sele:  :class:`~ost.mol.EntityView`
  :param proline: A selection containing only the proline (or another residue)
                  around which the helix is kinked. If False, the proline will
                  be searched for automatically
  :type proline:  :class:`ost.mol.EntityView`

  :return: A tuple (bend_angle, face_shift, wobble_angle).
  :rtype: (float, float, float)
  """
  (proline,proline_ca)=__FindProline(sele,proline)
  proline_num=proline.GetNumber().num
  (pre_proline,post_proline)=__SelectPreAndPostProline(sele,proline_num)
  (ca_3,ca_4)=__FindCa3AndCa4(sele,proline_ca,proline_num)  
  #Here we extract the necessary information from the structure
  pre_proline_axis=ost.mol.alg.structure_analysis.CalculateHelixAxis(pre_proline)
  post_proline_axis=ost.mol.alg.structure_analysis.CalculateHelixAxis(post_proline)
  prepa=pre_proline_axis.GetDirection()
  prepc=pre_proline_axis.GetOrigin()
  postpa=post_proline_axis.GetDirection()
  postpc=post_proline_axis.GetOrigin()
  #calculate the angles
  bend=__CalculateBendAngle(prepa,postpa)
  wobble=__CalculateWobbleAngle(prepa,postpa,postpc,proline_ca.pos)
  shift=__CalculateFaceShift(prepa,postpa,prepc,postpc,proline_ca.pos,ca_3.pos,ca_4.pos,bend)
  return (bend,shift,wobble)
  
  
  
  
