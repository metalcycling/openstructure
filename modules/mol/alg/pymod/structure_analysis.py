"""
Some functions for analyzing structures

Author: Niklaus Johner (Niklaus.Johner@unibas.ch)
"""
import os
import ost

def GetFrameFromEntity(eh):
  """
  This function returns a CoordFrame from an EntityHandle
  
  :param eh: 
  :type eh: :class:`~ost.mol.EntityHandle`

  :return: :class:`ost.mol.CoordFrame`
  """
  return ost.mol.CreateCoordFrame(eh.GetAtomPosList(ordered_by_index=True))
  
def GetDistanceBetwCenterOfMass(sele1,sele2):
  """
  This function calculates the distance between the centers of mass
  of **sele1** and **sele2**, two selections from the same Entity.

  :param sele1:
  :param sele2:
  :type sele1: :class:`~ost.mol.EntityView`
  :type sele2: :class:`~ost.mol.EntityView`

  :return: :class:`float`
  """
  if not sele1.IsValid() and sele2.IsValid():
    print('invalid view')
    return
  eh=sele1.GetHandle()
  if not eh==sele2.GetHandle():
    print('The two views must be from the same entity')
    return
  f=GetFrameFromEntity(eh)
  return f.GetDistanceBetwCenterOfMass(sele1,sele2)

def GetMinDistanceBetweenViews(sele1,sele2):
  """
  This function calculates the minimal distance between
  **sele1** and **sele2**, two selections from the same Entity.
  
  :param sele1:
  :param sele2:
  :type sele1: :class:`~ost.mol.EntityView`
  :type sele2: :class:`~ost.mol.EntityView`

  :return: :class:`float`
  """
  if not sele1.IsValid() and sele2.IsValid():
    print('invalid view')
    return
  eh=sele1.GetHandle()
  if not eh==sele2.GetHandle():
    print('The two views must be from the same entity')
    return
  f=GetFrameFromEntity(eh)
  return f.GetMinDistance(sele1,sele2)

def GetMinDistBetwCenterOfMassAndView(sele1,sele2):
  """
  This function calculates the minimal distance between **sele2** and
  the center of mass of **sele1**, two selections from the same Entity.

  :param sele1: The selection from which the center of mass is taken
  :param sele2:
  :type sele1: :class:`~ost.mol.EntityView`
  :type sele2: :class:`~ost.mol.EntityView`

  :return: distance (\\ :class:`float`\\ )
  """
  if not sele1.IsValid() and sele2.IsValid():
    print('invalid view')
    return
  eh=sele1.GetHandle()
  if not eh==sele2.GetHandle():
    print('The two views must be from the same entity')
    return
  f=GetFrameFromEntity(eh)
  return f.GetMinDistBetwCenterOfMassAndView(sele1,sele2)
  

def GetAlphaHelixContent(sele1):
  """
  This function calculates the content of alpha helix in a view.
  All residues in the view have to ordered and adjacent (no gaps allowed)
  
  :param sele1:
  :type sele1: :class:`~ost.mol.EntityView`

  :return: :class:`float`
  """
  if not sele1.IsValid():
    print('invalid view')
    return
  eh=sele1.GetHandle()
  f=GetFrameFromEntity(eh)
  return f.GetAlphaHelixContent(sele1)


def CalculateBestFitLine(sele1):
  """
  This function calculates the best fit line to the atoms in **sele1**.
  
  :param sele1:
  :type sele1: :class:`~ost.mol.EntityView`

  :return: :class:`~ost.geom.Line3`
  """
  if not sele1.IsValid():
    print('invalid view')
    return
  eh=sele1.GetHandle()
  f=GetFrameFromEntity(eh)
  return f.GetODRLine(sele1)

def CalculateBestFitPlane(sele1):
  """
  This function calculates the best fit plane to the atoms in **sele1**.
  
  :param sele1:
  :type sele1: :class:`~ost.mol.EntityView`

  :return: :class:`~ost.geom.Plane`
  """
  if not sele1.IsValid():
    print('invalid view')
    return
  eh=sele1.GetHandle()
  f=GetFrameFromEntity(eh)
  return f.GetODRPlane(sele1)

def CalculateHelixAxis(sele1):
  """
  This function calculates the best fit cylinder to the CA atoms in **sele1**,
  and returns its axis.  Residues should be ordered correctly
  in **sele1**.
  
  :param sele1:
  :type sele1: :class:`~ost.mol.EntityView`

  :return: :class:`~ost.geom.Line3`
  """
  if not sele1.IsValid():
    print('invalid view')
    return
  eh=sele1.GetHandle()
  f=GetFrameFromEntity(eh)
  return f.FitCylinder(sele1)[0]


def CalculateDistanceDifferenceMatrix(sele1,sele2):
  """
  This function calculates the pairwise distance differences between two selections (\\ :class:`~ost.mol.EntityView`\\ ).
  The two selections should have the same number of atoms
  It returns an NxN DistanceDifferenceMatrix M (where N is the number of atoms in sele1)
  where M[i,j]=||(sele2.atoms[i].pos-sele2.atoms[j].pos)||-||(sele1.atoms[i].pos-sele1.atoms[j].pos)||

  :param sele1:
  :param sele2:
  :type sele1: :class:`~ost.mol.EntityView`
  :type sele2: :class:`~ost.mol.EntityView`

  :return: NxN numpy matrix
  """
  try:import numpy as npy
  except ImportError:
    LogError("Function needs numpy, but I could not import it.")
    raise
  if not sele1.IsValid() and sele2.IsValid():
    print('invalid view')
    return
  if not sele1.GetAtomCount()==sele2.GetAtomCount():
    print('The two views must have the same number of atoms')
    return
  n_atoms=sele1.GetAtomCount()
  M=npy.zeros([n_atoms,n_atoms])
  for i,a1 in enumerate(sele1.atoms):
    for j,a2 in enumerate(sele1.atoms):
      if i>=j:continue
      d1=ost.geom.Distance(a1.pos,a2.pos)
      d2=ost.geom.Distance(sele2.atoms[i].pos,sele2.atoms[j].pos)
      M[i,j]=d2-d1
      M[j,i]=d2-d1
  return M


