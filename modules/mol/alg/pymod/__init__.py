import os.path
from ._ost_mol_alg import *
from ost.mol.alg.superpose import *
import ost.mol.alg.trajectory_analysis
import ost.mol.alg.structure_analysis
import ost.mol.alg.helix_kinks
import ost.mol.alg.hbond

# Fills a list of reference clashing distances from a file (requires a path to the file)
def FillClashingDistancesFromFile(filename):
  fh=open(filename,'r')
  lines=fh.readlines()
  fh.close()
  return FillClashingDistances(lines)

# Fills a list of bond stereo-chemical statistics from a file (requires a path to the file)
def FillBondStereoChemicalParamsFromFile(filename):
  fh=open(filename,'r')
  lines=fh.readlines()
  fh.close()
  return FillBondStereoChemicalParams("Bond",lines)

# Fills a list of angle stereo-chemical statistics from a file (requires a path to the file)
def FillAngleStereoChemicalParamsFromFile(filename):
  fh=open(filename,'r')
  lines=fh.readlines()
  fh.close()
  return FillAngleStereoChemicalParams("Angle",lines)

# Returns the default list of reference clashing distances (from the default OpenStructure parameter file)
def DefaultClashingDistances():
  shared_path=ost.GetSharedDataPath()	
  filename=os.path.join(shared_path,'stereo_chemical_props.txt')
  fh=open(filename,'r')
  lines=fh.readlines()
  fh.close()
  return FillClashingDistances(lines)

# Returns the default list of bond stereo-chemical statistics (from the default OpenStructure parameter file)
def DefaultBondStereoChemicalParams():
  shared_path=ost.GetSharedDataPath()	
  filename=os.path.join(shared_path,'stereo_chemical_props.txt')
  fh=open(filename,'r')
  lines=fh.readlines()
  fh.close()
  return FillStereoChemicalParams("Bond",lines)

# Returns the default list of angle stereo-chemical statistics (from the default OpenStructure parameter file)
def DefaultAngleStereoChemicalParams():
  shared_path=ost.GetSharedDataPath()	
  filename=os.path.join(shared_path,'stereo_chemical_props.txt')
  fh=open(filename,'r')
  lines=fh.readlines()
  fh.close()
  return FillStereoChemicalParams("Angle",lines)

