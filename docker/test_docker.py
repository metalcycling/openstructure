import ost
from ost.mol.alg import scoring, qsscoring
from ost import conop

lib = conop.GetDefaultLib()
if lib is not None:
    print("You have a valid compound library, last updated on " +
            lib.GetCreationDate())
else:
    ost.LogError("No compound library set as default!")
    print("The compound library is not working properly!")

# load two biounits to compare
ent_full = ost.io.LoadPDB('3ia3', remote=True)
ent_1 = ent_full.Select('cname=A,D')
ent_2 = ent_full.Select('cname=B,C')

# get scores
ost.PushVerbosityLevel(3)    

try:
    scorer = scoring.Scorer(ent_1, ent_2)
    ost.LogScript('lDDT:', str(scorer.lddt))
    ost.LogScript('QSscore:', str(scorer.qs_global))
    ost.LogScript('Chain mapping used:', str(scorer.mapping.GetFlatMapping()))
except qsscoring.QSscoreError as ex:
    # default handling: report failure and set score to 0
    ost.LogError('QSscore failed:', str(ex))
    qs_score = 0
    print("OST is not working properly!")
else:
    print("OST is working!")

