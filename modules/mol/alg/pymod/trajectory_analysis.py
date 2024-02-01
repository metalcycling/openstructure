"""
**This Module requires numpy**

This module contains functions to analyze trajectories, mainly 
similiraty measures baed on RMSDS and pairwise distances.

Author: Niklaus Johner (niklaus.johner@unibas.ch)
"""

import ost.mol.alg
import ost.geom
from ost import LogError
import os

def smooth(vec,n):
#Function to smooth a vector or a list of floats
#for each element it takes the average over itself and the
#n elements on each side, so over (2n+1) elements
  try:
    vec2=vec.copy()
  except:
    vec2=vec[:]
  for i in range(n):
    v=0.0
    count=1.0
    v+=vec[i]
    for j in range(n):
      count+=1
      v+=vec[i+j+1]
    for j in range(i):
      count+=1
      v+=vec[i-(j+1)]
    vec2[i]=v/float(count)
  for i in range(1,n+1):
    v=0.0
    count=1.0
    v+=vec[-i]
    for j in range(n):
      count+=1
      v+=vec[-(i+j+1)]
    for j in range(i-1):
      count+=1
      v+=vec[-i+j+1]
    vec2[-i]=v/float(count)
  for i in range(n,len(vec2)-n):
    v=vec[i]
    for j in range(n):
      v+=vec[i+j+1]
      v+=vec[i-j-1]
    vec2[i]=v/float(2.*n+1.)
  return vec2


"""
From here on the module needs numpy
"""

def RMSD_Matrix_From_Traj(t,sele,first=0,last=-1,align=True,align_sele=None):
  """
  This function calculates a matrix M such that M[i,j] is the
  RMSD (calculated on **sele**) between frames i and j of the trajectory **t**
  aligned on sele.

  :param t: the trajectory
  :param sele: the selection used for alignment and RMSD calculation
  :param first: the first frame of t to be used
  :param last: the last frame of t to be used
  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type sele: :class:`~ost.mol.EntityView`
  :type first: :class:`int`
  :type last: :class:`int`

  :return: Returns a numpy N\\ :subscript:`frames`\\ xN\\ :subscript:`frames` matrix,
   where N\\ :subscript:`frames` is the number of frames.
  """
  if not align_sele:align_sele=sele
  try:
    import numpy as npy
    if last==-1:last=t.GetFrameCount()
    n_frames=last-first
    rmsd_matrix=npy.identity(n_frames)
    for i in range(n_frames):
      if align:
        t=ost.mol.alg.SuperposeFrames(t,align_sele,begin=first,end=last,ref=i)
        eh=t.GetEntity()
      t.CopyFrame(i)
      rmsd_matrix[i,:]=ost.mol.alg.AnalyzeRMSD(t,sele,sele)
      if i==0:
        last=last-first
        first=0
    return rmsd_matrix
  except ImportError:
    LogError("Function needs numpy, but I could not import it.")
    raise


def PairwiseDistancesFromTraj(t,sele,first=0,last=-1,seq_sep=1):
  """
  This function calculates the distances between any pair of atoms in **sele**  
  with sequence separation larger than **seq_sep** from a trajectory **t**.
  It return a matrix containing one line for each atom pair and N\\ :subscript:`frames` columns, where
  N\\ :subscript:`frames` is the number of frames in the trajectory.
  
  :param t: the trajectory
  :param sele: the selection used to determine the atom pairs
  :param first: the first frame of t to be used
  :param last: the last frame of t to be used
  :param seq_sep: The minimal sequence separation between atom pairs
  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type sele: :class:`~ost.mol.EntityView`
  :type first: :class:`int`
  :type last: :class:`int`
  :type seq_sep: :class:`int`

  :return: a numpy N\\ :subscript:`pairs`\\ xN\\ :subscript:`frames` matrix.
  """
  try:
    import numpy as npy
    if last==-1:last=t.GetFrameCount()
    n_frames=last-first
    n_var=0
    for i,a1 in enumerate(sele.atoms):
      for j,a2 in enumerate(sele.atoms):
        if not j-i<seq_sep:n_var+=1
    #n_var=sele.GetAtomCount()
    #n_var=(n_var-1)*(n_var)/2.
    dist_matrix=npy.zeros(n_frames*n_var)
    dist_matrix=dist_matrix.reshape(n_var,n_frames)
    k=0
    for i,a1 in enumerate(sele.atoms):
      for j,a2 in enumerate(sele.atoms):
        if j-i<seq_sep:continue
        dist_matrix[k]=ost.mol.alg.AnalyzeDistanceBetwAtoms(t,a1.GetHandle(),a2.GetHandle())[first:last]
        k+=1
    return dist_matrix
  except ImportError:
    LogError("Function needs numpy, but I could not import it.")
    raise
    
def DistanceMatrixFromPairwiseDistances(distances,p=2):
  """
  This function calculates an distance matrix M(N\\ :subscript:`frames`\\ xN\\ :subscript:`frames`\\ ) from
  the pairwise distances matrix D(N\\ :subscript:`pairs`\\ xN\\ :subscript:`frames`\\ ), where
  N\\ :subscript:`frames` is the number of frames in the trajectory
  and N\\ :subscript:`pairs` the number of atom pairs.
  M[i,j] is the distance between frame i and frame j
  calculated as a p-norm of the differences in distances
  from the two frames (distance-RMSD for p=2).

  :param distances: a pairwise distance matrix as obtained from 
   :py:func:`~mol.alg.trajectory_analysis.PairwiseDistancesFromTraj`
  :param p: exponent used for the p-norm.

  :return: a numpy N\\ :subscript:`frames`\\ xN\\ :subscript:`frames` matrix, where N\\ :subscript:`frames`
   is the number of frames.
  """
  try:
    import numpy as npy
    n1=distances.shape[0]
    n2=distances.shape[1]
    dist_mat=npy.identity(n2)
    for i in range(n2):
      for j in range(n2):
        if j<=i:continue
        d=(((abs(distances[:,i]-distances[:,j])**p).sum())/float(n1))**(1./p)
        dist_mat[i,j]=d
        dist_mat[j,i]=d
    return dist_mat
  except ImportError:
    LogError("Function needs numpy, but I could not import it.")
    raise

def DistRMSDFromTraj(t,sele,ref_sele,radius=7.0,average=False,seq_sep=4,first=0,last=-1):
  """
  This function calculates the distance RMSD from a trajectory.
  The distances selected for the calculation are all the distances
  between pair of atoms from residues that are at least **seq_sep** apart
  in the sequence and that are smaller than **radius** in **ref_sel**.
  The number and order of atoms in **ref_sele** and **sele** should be the same.

  :param t: the trajectory
  :param sele: the selection used to calculate the distance RMSD
  :param ref_sele: the reference selection used to determine the atom pairs and reference distances
  :param radius: the upper limit of distances in ref_sele considered for the calculation
  :param seq_sep: the minimal sequence separation between atom pairs considered for the calculation 
  :param average: use the average distance in the trajectory as reference instead of the distance obtained from ref_sele
  :param first: the first frame of t to be used
  :param last: the last frame of t to be used
  
  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type sele: :class:`~ost.mol.EntityView`
  :type ref_sele: :class:`~ost.mol.EntityView`
  :type radius: :class:`float`
  :type average: :class:`bool`
  :type first: :class:`int`
  :type last: :class:`int`
  :type seq_sep: :class:`int`

  :return: a numpy vecor dist_rmsd(N\\ :subscript:`frames`).
  """
  if not sele.GetAtomCount()==ref_sele.GetAtomCount():
    print('Not same number of atoms in the two views')
    return
  try:
    import numpy as npy
    if last==-1:last=t.GetFrameCount()
    n_frames=last-first
    dist_rmsd=npy.zeros(n_frames)
    pair_count=0.0
    for i,a1 in enumerate(ref_sele.atoms):
      for j,a2 in enumerate(ref_sele.atoms):
        if j<=i:continue
        r1=a1.GetResidue()
        c1=r1.GetChain()
        r2=a2.GetResidue()
        c2=r2.GetChain()      
        if c1==c2 and abs(r2.GetNumber().num-r1.GetNumber().num)<seq_sep:continue
        d=ost.geom.Distance(a1.pos,a2.pos)
        if d<radius:
          a3=sele.atoms[i]
          a4=sele.atoms[j]
          d_traj=ost.mol.alg.AnalyzeDistanceBetwAtoms(t,a3.GetHandle(),a4.GetHandle())[first:last]
          if average:d=npy.mean(d_traj)
          for k,el in enumerate(d_traj):
            dist_rmsd[k]+=(el-d)**2.0
          pair_count+=1.0
    return (dist_rmsd/float(pair_count))**0.5
  except ImportError:
    LogError("Function needs numpy, but I could not import it.")
    raise
    
def AverageDistanceMatrixFromTraj(t,sele,first=0,last=-1):
  """
  This function calcultes the distance between each pair of atoms
  in **sele**, averaged over the trajectory **t**.

  :param t: the trajectory
  :param sele: the selection used to determine the atom pairs
  :param first: the first frame of t to be used
  :param last: the last frame of t to be used
  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type sele: :class:`~ost.mol.EntityView`
  :type first: :class:`int`
  :type last: :class:`int`

  :return: a numpy N\\ :subscript:`pairs`\\ xN\\ :subscript:`pairs` matrix, where N\\ :subscript:`pairs`
   is the number of atom pairs in **sele**.
  """
  try:
    import numpy as npy
  except ImportError:
    LogError("Function needs numpy, but I could not import it.")
    raise
  n_atoms=sele.GetAtomCount()
  M=npy.zeros([n_atoms,n_atoms])
  for i,a1 in enumerate(sele.atoms):
    for j,a2 in enumerate(sele.atoms):
      if j>i:continue
      d=ost.mol.alg.AnalyzeDistanceBetwAtoms(t,a1.GetHandle(),a2.GetHandle())[first:last]
      M[i,j]=npy.mean(d)
      M[j,i]=npy.mean(d)
  return M

def AnalyzeDistanceFluctuationMatrix(t,sele,first=0,last=-1):
  try:
    import numpy as npy
  except ImportError:
    LogError("Function needs numpy, but I could not import it.")
    raise
  n_atoms=sele.GetAtomCount()
  M=npy.zeros([n_atoms,n_atoms])
  for i,a1 in enumerate(sele.atoms):
    for j,a2 in enumerate(sele.atoms):
      if i>j:continue
      d=ost.mol.alg.AnalyzeDistanceBetwAtoms(t,a1.GetHandle(),a2.GetHandle())[first:last]
      M[j,i]=npy.std(d)
      M[i,j]=npy.std(d)
  return M
  
def IterativeSuperposition(t,sele,threshold=1.0,initial_sele=None,iterations=5,ref_frame=0):
  if initial_sele:current_sele=initial_sele
  else: current_sele=sele
  for i in range(iterations):
    t=ost.mol.alg.SuperposeFrames(t,current_sele,ref=ref_frame)
    al=[a for a in sele.atoms if ost.mol.alg.AnalyzeRMSF(t,ost.mol.CreateViewFromAtoms([a]))<threshold]
    if len(al)==0:return
    current_sele=ost.mol.CreateViewFromAtoms(al)
  return current_sele
    