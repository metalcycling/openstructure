from ost.seq.alg import SubstWeightMatrix

def _InitMatrix(preset_enum):
  mat = SubstWeightMatrix()
  mat.AssignPreset(preset_enum)
  return mat

BLOSUM45 = _InitMatrix(SubstWeightMatrix.Preset.BLOSUM45)
BLOSUM62 = _InitMatrix(SubstWeightMatrix.Preset.BLOSUM62)
BLOSUM80 = _InitMatrix(SubstWeightMatrix.Preset.BLOSUM80)
BLOSUM100 = _InitMatrix(SubstWeightMatrix.Preset.BLOSUM100)
IDENTITY = _InitMatrix(SubstWeightMatrix.Preset.IDENTITY)
MATCH = _InitMatrix(SubstWeightMatrix.Preset.MATCH)
NUC44 = _InitMatrix(SubstWeightMatrix.Preset.NUC44)

__all__=['BLOSUM45','BLOSUM62','BLOSUM80','BLOSUM100', 'IDENTITY', 'MATCH',
         'NUC44']
