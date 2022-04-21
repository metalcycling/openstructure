..  Note on large code blocks: keep max. width to 100 or it will look bad
                               on webpage!
..  TODO: look at argparse directive to autogenerate --help output!

.. ost-actions:

OST Actions
================================================================================

A pure command line interface of OST is provided by actions.
You can execute ``ost -h`` for a list of possible actions and for every action,
you can type ``ost <ACTION> -h`` to get a description on its usage.

Here we list the most prominent actions with simple examples.

.. _ost compare structures:

Comparing two structures
--------------------------------------------------------------------------------

You can compare two structures in terms of quaternary structure score and
lDDT scores between two complexes from the command line with the
``ost compare-structures`` action.

In summary it performs the following steps:

- Read structures (PDB or mmCIF format, can be gzipped) and split into
  biological assemblies (all possible pairs are scored).
- Mandatory cleanup of hydrogen atoms, ligand and water chains, small
  (< 20 residues) peptides or chains with no amino acids as described in
  :attr:`QSscoreEntity.ent <ost.mol.alg.qsscoring.QSscoreEntity.ent>` and
  :attr:`QSscoreEntity.removed_chains <ost.mol.alg.qsscoring.QSscoreEntity.removed_chains>`.
- Optional cleanup of structures with :func:`~ost.mol.alg.Molck`.
- Optional structural checks with :func:`~ost.mol.alg.CheckStructure`.
- Unless user-provided, find chain mapping between complexes (see
  :attr:`here <ost.mol.alg.qsscoring.QSscorer.chain_mapping>` for details)
- Perform sequence alignment of chain pairs (unless user asks for alignment
  based on residue numbers). Alignment can optionally checked for consistency
  if 100% sequence identity is expected.
- Compute scores requested by user (CA-RMSD of oligomer,
  :mod:`QS scores of oligomer <ost.mol.alg.qsscoring>`,
  :class:`single chain lDDT scores <ost.mol.alg.qsscoring.OligoLDDTScorer.sc_lddt>`,
  :attr:`lDDT score of oligomer <ost.mol.alg.qsscoring.OligoLDDTScorer.oligo_lddt>`).
  Note that while the QS score is symmetric (same result when swapping reference
  and model), the lDDT scores are not. Extra atoms in the model for mapped
  chains have no effect on the score, while extra atoms in the reference reduce
  the score. For oligo_lddt, we do
  :attr:`penalize for extra chains <ost.mol.alg.qsscoring.OligoLDDTScorer.penalize_extra_chains>`
  in both reference and model.

.. note ::

  The action relies on OST's :mod:`~ost.mol.alg.qsscoring` module and has the
  same requirements on your OST installation (needs compound library, ClustalW,
  numpy and scipy).

Details on the usage (output of ``ost compare-structures --help``):

.. code-block:: console

  usage: ost compare-structures [-h] -m MODEL -r REFERENCE [-v VERBOSITY]
                                [-o OUTPUT] [-d] [-ds DUMP_SUFFIX]
                                [-rs REFERENCE_SELECTION] [-ms MODEL_SELECTION]
                                [-ca] [-ft] [-cl COMPOUND_LIBRARY] [-ml]
                                [-rm REMOVE [REMOVE ...]] [-ce] [-mn] [-sc]
                                [-p PARAMETER_FILE] [-bt BOND_TOLERANCE]
                                [-at ANGLE_TOLERANCE]
                                [-c CHAIN_MAPPING [CHAIN_MAPPING ...]]
                                [--qs-max-mappings-extensive QS_MAX_MAPPINGS_EXTENSIVE]
                                [-cc] [-rna] [-qs] [--qs-rmsd] [-l]
                                [-ir INCLUSION_RADIUS] [-ss SEQUENCE_SEPARATION]
                                [-spr]

  Evaluate model structure against reference.

  eg.

    ost compare-structures \
        --model <MODEL> \
        --reference <REF> \
        --output output.json \
        --lddt \
        --structural-checks \
        --consistency-checks \
        --molck \
        --remove oxt hyd \
        --map-nonstandard-residues

  Here we describe how the parameters can be set to mimick a CAMEO evaluation
  (as of August 2018).

  CAMEO calls the lddt binary as follows:

    lddt \
        -p <PARAMETER FILE> \
        -f \
        -a 15 \
        -b 15 \
        -r 15 \
        <MODEL> \
        <REF>

  Only model structures are "Molck-ed" in CAMEO. The call to molck is as follows:

    molck \
        --complib=<COMPOUND LIB> \
        --rm=hyd,oxt,unk,nonstd \
        --fix-ele \
        --map-nonstd \
        --out=<OUTPUT> \
        <FILEPATH>

  To be as much compatible with with CAMEO as possible one should call
  compare-structures as follows:

    ost compare-structures \
        --model <MODEL> \
        --reference <REF> \
        --output output.json \
        --molck \
        --remove oxt hyd unk nonstd \
        --clean-element-column \
        --map-nonstandard-residues \
        --structural-checks \
        --bond-tolerance 15.0 \
        --angle-tolerance 15.0 \
        --residue-number-alignment \
        --consistency-checks \
        --qs-score \
        --lddt \
        --inclusion-radius 15.0

  optional arguments:
    -h, --help            show this help message and exit

  required arguments:
    -m MODEL, --model MODEL
                          Path to the model file.
    -r REFERENCE, --reference REFERENCE
                          Path to the reference file.

  general arguments:
    -v VERBOSITY, --verbosity VERBOSITY
                          Set verbosity level. Defaults to 3.
    -o OUTPUT, --output OUTPUT
                          Output file name. The output will be saved as a JSON file.
    -d, --dump-structures
                          Dump cleaned structures used to calculate all the scores as
                          PDB files using specified suffix. Files will be dumped to the
                          same location as original files.
    -ds DUMP_SUFFIX, --dump-suffix DUMP_SUFFIX
                          Use this suffix to dump structures.
                          Defaults to .compare.structures.pdb.
    -rs REFERENCE_SELECTION, --reference-selection REFERENCE_SELECTION
                          Selection performed on reference structures.
    -ms MODEL_SELECTION, --model-selection MODEL_SELECTION
                          Selection performed on model structures.
    -ca, --c-alpha-only   Use C-alpha atoms only. Equivalent of calling the action with
                          '--model-selection="aname=CA" --reference-selection="aname=CA"'
                          options.
    -ft, --fault-tolerant
                          Fault tolerant parsing.
    -cl COMPOUND_LIBRARY, --compound-library COMPOUND_LIBRARY
                          Location of the compound library file (compounds.chemlib).
                          If not provided, the following locations are searched in this
                          order: 1. Working directory, 2. OpenStructure standard library
                          location.

  molecular check arguments:
    -ml, --molck          Run molecular checker to clean up input.
    -rm REMOVE [REMOVE ...], --remove REMOVE [REMOVE ...]
                          Remove atoms and residues matching some criteria:
                           * zeroocc - Remove atoms with zero occupancy
                           * hyd - remove hydrogen atoms
                           * oxt - remove terminal oxygens
                           * nonstd - remove all residues not one of the 20
                                      standard amino acids
                           * unk - Remove unknown and atoms not following the
                                   nomenclature
                          Defaults to hyd.
    -ce, --clean-element-column
                          Clean up element column
    -mn, --map-nonstandard-residues
                          Map modified residues back to the parent amino acid, for
                          example MSE -> MET, SEP -> SER.

  structural check arguments:
    -sc, --structural-checks
                          Perform structural checks and filter input data.
    -p PARAMETER_FILE, --parameter-file PARAMETER_FILE
                          Location of the stereochemical parameter file
                          (stereo_chemical_props.txt).
                          If not provided, the following locations are searched in this
                          order: 1. Working directory, 2. OpenStructure standard library
                          location.
    -bt BOND_TOLERANCE, --bond-tolerance BOND_TOLERANCE
                          Tolerance in STD for bonds. Defaults to 12.
    -at ANGLE_TOLERANCE, --angle-tolerance ANGLE_TOLERANCE
                          Tolerance in STD for angles. Defaults to 12.

  chain mapping arguments:
    -c CHAIN_MAPPING [CHAIN_MAPPING ...], --chain-mapping CHAIN_MAPPING [CHAIN_MAPPING ...]
                          Mapping of chains between the reference and the model.
                          Each separate mapping consist of key:value pairs where key
                          is the chain name in reference and value is the chain name in
                          model.
    --qs-max-mappings-extensive QS_MAX_MAPPINGS_EXTENSIVE
                          Maximal number of chain mappings to test for 'extensive'
                          chain mapping scheme which is used as a last resort if
                          other schemes failed. The extensive chain mapping search
                          must in the worst case check O(N!) possible mappings for
                          complexes with N chains. Two octamers without symmetry
                          would require 322560 mappings to be checked. To limit
                          computations, no scores are computed if we try more than
                          the maximal number of chain mappings. Defaults to 1000000.

  sequence alignment arguments:
    -cc, --consistency-checks
                          Take consistency checks into account. By default residue name
                          consistency between a model-reference pair would be checked
                          but only a warning message will be displayed and the script
                          will continue to calculate scores. If this flag is ON, checks
                          will not be ignored and if the pair does not pass the test
                          all the scores for that pair will be marked as a FAILURE.
    -rna, --residue-number-alignment
                          Make alignment based on residue number instead of using
                          a global BLOSUM62-based alignment.

  QS score arguments:
    -qs, --qs-score       Calculate QS-score.
    --qs-rmsd             Calculate CA RMSD between shared CA atoms of mapped chains.
                          This uses a superposition using all mapped chains which
                          minimizes the CA RMSD.

  lDDT score arguments:
    -l, --lddt            Calculate lDDT.
    -ir INCLUSION_RADIUS, --inclusion-radius INCLUSION_RADIUS
                          Distance inclusion radius for lDDT. Defaults to 15 A.
    -ss SEQUENCE_SEPARATION, --sequence-separation SEQUENCE_SEPARATION
                          Sequence separation. Only distances between residues whose
                          separation is higher than the provided parameter are
                          considered when computing the score. Defaults to 0.
    -spr, --save-per-residue-scores


By default the verbosity is set to 3 which will result in the informations
being shown in the console. The result can be (optionally) saved as JSON file
which is the preferred way of parsing it as the log output might change in the
future. Optionally, the local scores for lDDT can also be dumped to the output
file. Additionally, cleaned up structures can be saved to the disk.
The output file has following format:

.. code-block:: none

  {
    "options": { ... },  # Options used to run the script
    "result": {
      "<MODEL NAME>": { # Model name extracted from the file name
        "<REFERENCE NAME>": { # Reference name extracted from the file name
          "info": {
            "mapping": {
              "alignments": <list of chain-chain alignments in FASTA format>,
              "chain_mapping": <Mapping of chains eg. {"A": "B", "B": "A"}>,
              "chain_mapping_scheme": <Scheme used to get mapping, check mapping manually
                                       if "permissive" or "extensive">
            },
            "residue_names_consistent": <Are the residue numbers consistent? true or false>
          },
          "lddt": {
            # calculated when --lddt (-l) option is selected
            "oligo_lddt": {
              "error": <ERROR message if any>,
              "global_score": <calculated oligomeric lDDT score>,
              "status": <SUCCESS or FAILURE>
            },
            "single_chain_lddt": [
              # a list of chain-chain lDDTs
              {
                "conserved_contacts": <number of conserved contacts between model & reference>,
                "error": <ERROR message if any>,
                "global_score": <calculated single-chain lDDT score>,
                "model_chain": <name of the chain in model>,
                "reference_chain": <name of the chain in reference>,
                "status": <SUCCESS or FAILURE>,
                "total_contacts": <total number of contacts in reference>,
                "per_residue_scores": [
                  # per-residue lDDT scores
                  # only calculated when --save-per-residue-scores (-spr) option is selected
                  {
                    "residue_name": <three letter code of the residue in reference chain>,
                    "residue_number": <residue number in reference chain>,
                    "lddt": <residue lDDT score>,
                    "conserved_contacts": <conserved_contacts for given residue>,
                    "total_contacts": <total_contacts for given residue>
                  },
                  .
                  .
                  .
                ]
              }
            ],
            "weighted_lddt": {
              "error": <ERROR message if any>,
              "global_score": <calculated weighted lDDT score>,
              "status": <SUCCESS or FAILURE>
            }
          },
          "qs_score": {
            # calculated when --qs-score (-q) option is selected
            "best_score": <Best QS-score>,
            "error": <ERROR message if any>,
            "global_score": <Global QS-score>,
            "status": <SUCCESS or FAILURE>
          }
        }
      }
    }
  }

The "result" filed is a dictionary mapping from model to reference as eg. in
mmCIF file there can be many entities and the script will compare all
combinations.

Example usage:

.. code-block:: console

  $ CAMEO_TARGET_URL=https://www.cameo3d.org/static/data/modeling/2019.07.13/6PO4_F
  $ curl $CAMEO_TARGET_URL/bu_target_01.pdb > reference.pdb
  $ curl $CAMEO_TARGET_URL/servers/server20/oligomodel-1/oligomodel-1.pdb > model.pdb
  $ $OST_ROOT/bin/ost compare-structures \
        --model model.pdb --reference reference.pdb --output output.json \
        --qs-score --residue-number-alignment --lddt --structural-checks \
        --consistency-checks --inclusion-radius 15.0 --bond-tolerance 15.0 \
        --angle-tolerance 15.0 --molck --remove oxt hyd unk nonstd \
        --clean-element-column --map-nonstandard-residues

  ################################################################################
  Reading input files (fault_tolerant=False)
   --> reading model from model.pdb
  imported 2 chains, 462 residues, 3400 atoms; with 0 helices and 0 strands
   --> reading reference from reference.pdb
  imported 3 chains, 471 residues, 3465 atoms; with 0 helices and 0 strands
  ################################################################################
  Cleaning up input with Molck
  removing hydrogen atoms
   --> removed 0 hydrogen atoms
  removing OXT atoms
   --> removed 3 OXT atoms
  _.HCS1 is not a standard amino acid --> removed 
  _.ADE2 is not a standard amino acid --> removed 
  _.BO33 is not a standard amino acid --> removed 
  _.ADE4 is not a standard amino acid --> removed 
  _.HCS5 is not a standard amino acid --> removed 
  _.BO36 is not a standard amino acid --> removed 
  removing hydrogen atoms
   --> removed 0 hydrogen atoms
  removing OXT atoms
   --> removed 0 OXT atoms
  ################################################################################
  Performing structural checks
   --> for reference(s)
  Checking reference.pdb
  Checking stereo-chemistry
  Average Z-Score for bond lengths: 0.33163
  Bonds outside of tolerance range: 0 out of 2993
  Bond  Avg Length  Avg zscore  Num Bonds
  C-C 1.51236     0.03971     1682
  C-N 1.46198     0.96819     603
  C-O 1.25794     0.49967     674
  C-S 1.80242     0.15292     34
  Average Z-Score angle widths: -0.12077
  Angles outside of tolerance range: 0 out of 3260
  Filtering non-bonded clashes
  0 non-bonded short-range distances shorter than tolerance distance
  Distances shorter than tolerance are on average shorter by: 0.00000
   --> for model(s)
  Checking model.pdb
  Checking stereo-chemistry
  Average Z-Score for bond lengths: 0.23693
  Bonds outside of tolerance range: 0 out of 2976
  Bond  Avg Length  Avg zscore  Num Bonds
  C-C 1.52020     0.40359     1674
  C-N 1.43936     -0.19949    598
  C-O 1.25221     0.20230     670
  C-S 1.81182     0.38936     34
  Average Z-Score angle widths: 0.04946
  Angles outside of tolerance range: 0 out of 3241
  Filtering non-bonded clashes
  0 non-bonded short-range distances shorter than tolerance distance
  Distances shorter than tolerance are on average shorter by: 0.00000
  ################################################################################
  Comparing model.pdb to reference.pdb
  Chains in reference.pdb: AB
  Chains in model.pdb: AB
  Chemically equivalent chain-groups in reference.pdb: [['A', 'B']]
  Chemically equivalent chain-groups in model.pdb: [['A', 'B']]
  Chemical chain-groups mapping: {('A', 'B'): ('A', 'B')}
  Identifying Symmetry Groups...
  Symmetry threshold 0.1 used for angles of reference.pdb
  Symmetry threshold 0.1 used for axis of reference.pdb
  Symmetry threshold 0.1 used for angles of model.pdb
  Symmetry threshold 0.1 used for axis of model.pdb
  Selecting Symmetry Groups...
  Symmetry-groups used in reference.pdb: [('A',), ('B',)]
  Symmetry-groups used in model.pdb: [('A',), ('B',)]
  Closed Symmetry with strict parameters
  Mapping found: {'A': 'A', 'B': 'B'}
  --------------------------------------------------------------------------------
  Checking consistency between model.pdb and reference.pdb
  Consistency check: OK
  --------------------------------------------------------------------------------
  Computing QS-score
  QSscore reference.pdb, model.pdb: best: 0.96, global: 0.96
  --------------------------------------------------------------------------------
  Computing lDDT scores
  lDDT settings: 
  Inclusion Radius: 15
  Sequence separation: 0
  Cutoffs: 0.5, 1, 2, 4
  Residue properties label: lddt
  ===
   --> Computing lDDT between model chain A and reference chain A
  Coverage: 0.991416 (231 out of 233 residues)
  Global LDDT score: 0.8955
  (1194245 conserved distances out of 1333644 checked, over 4 thresholds)
   --> Computing lDDT between model chain B and reference chain B
  Coverage: 0.991379 (230 out of 232 residues)
  Global LDDT score: 0.8998
  (1200391 conserved distances out of 1334056 checked, over 4 thresholds)
   --> Computing oligomeric lDDT score
  Reference reference.pdb has: 2 chains
  Model model.pdb has: 2 chains
  Coverage: 0.991398 (461 out of 465 residues)
  Oligo lDDT score: 0.8977
   --> Computing weighted lDDT score
  Weighted lDDT score: 0.8976
  ################################################################################
  Saving output into output.json

This reads the model and reference file and calculates QS- and lDDT-scores
between them. In the example above the output file looks as follows (FASTA
alignments were cut in display here for readability):

..  code snippet to fix output.json generated above
  import json
  json_data = json.load(open("output.json"))
  mapping = json_data["result"]["model.pdb"]["reference.pdb"]["info"]["mapping"]
  new_alns = list()
  for aln in mapping["alignments"]:
    aln_lines = aln.splitlines()
    aln_lines[1] = aln_lines[1][:15] + "..."
    aln_lines[3] = aln_lines[3][:15] + "..."
    new_alns.append("\n".join(aln_lines))
  mapping["alignments"] = new_alns
  json_data["options"]["parameter_file"] = "Path to stage/share/openstructure/stereo_chemical_props.txt"
  json_data["options"]["compound_library"] = "Path to stage/share/openstructure/compounds.chemlib"
  json_data["options"]["cwd"] = "Path to current working directory"
  with open("output_fixed.json", "w") as outfile:
    json.dump(json_data, outfile, indent=2, sort_keys=True)

.. code-block:: json

  {
    "options": {
      "angle_tolerance": 15.0, 
      "bond_tolerance": 15.0, 
      "c_alpha_only": false, 
      "chain_mapping": null, 
      "clean_element_column": true, 
      "compound_library": "Path to stage/share/openstructure/compounds.chemlib", 
      "consistency_checks": true, 
      "cwd": "Path to current working directory", 
      "dump_structures": false, 
      "dump_suffix": ".compare.structures.pdb", 
      "fault_tolerant": false, 
      "inclusion_radius": 15.0, 
      "lddt": true, 
      "map_nonstandard_residues": true, 
      "model": "model.pdb", 
      "model_selection": "", 
      "molck": true, 
      "output": "output.json", 
      "parameter_file": "Path to stage/share/openstructure/stereo_chemical_props.txt", 
      "qs_max_mappings_extensive": 1000000, 
      "qs_rmsd": false, 
      "qs_score": true, 
      "reference": "reference.pdb", 
      "reference_selection": "", 
      "remove": [
        "oxt", 
        "hyd", 
        "unk", 
        "nonstd"
      ], 
      "residue_number_alignment": true, 
      "save_per_residue_scores": false, 
      "sequence_separation": 0, 
      "structural_checks": true, 
      "verbosity": 3
    }, 
    "result": {
      "model.pdb": {
        "reference.pdb": {
          "info": {
            "mapping": {
              "alignments": [
                ">reference:A\n-NAMKIGIVGAMAQE...\n>model:A\n---MKIGIVGAMAQE...", 
                ">reference:B\n-NAMKIGIVGAMAQE...\n>model:B\n---MKIGIVGAMAQE..."
              ], 
              "chain_mapping": {
                "A": "A", 
                "B": "B"
              }, 
              "chain_mapping_scheme": "strict"
            }, 
            "residue_names_consistent": true
          }, 
          "lddt": {
            "oligo_lddt": {
              "error": "", 
              "global_score": 0.8977285786061329, 
              "status": "SUCCESS"
            }, 
            "single_chain_lddt": [
              {
                "conserved_contacts": 1194245, 
                "error": "", 
                "global_score": 0.8954750895500183, 
                "model_chain": "A", 
                "reference_chain": "A", 
                "status": "SUCCESS", 
                "total_contacts": 1333644
              }, 
              {
                "conserved_contacts": 1200391, 
                "error": "", 
                "global_score": 0.8998055458068848, 
                "model_chain": "B", 
                "reference_chain": "B", 
                "status": "SUCCESS", 
                "total_contacts": 1334056
              }
            ], 
            "weighted_lddt": {
              "error": "", 
              "global_score": 0.8976406520766181, 
              "status": "SUCCESS"
            }
          }, 
          "qs_score": {
            "best_score": 0.9619749105661133, 
            "error": "", 
            "global_score": 0.9619749105661133, 
            "status": "SUCCESS"
          }
        }
      }
    }
  }

If all the structures are clean and have matching residue numbers, one can omit
all the checking steps and calculate scores directly as here:

.. code:: console

  $ $OST_ROOT/bin/ost compare-structures \
        --model model.pdb --reference reference.pdb --output output_qs.json \
        --qs-score --residue-number-alignment

  ################################################################################
  Reading input files (fault_tolerant=False)
   --> reading model from model.pdb
  imported 2 chains, 462 residues, 3400 atoms; with 0 helices and 0 strands
   --> reading reference from reference.pdb
  imported 3 chains, 471 residues, 3465 atoms; with 0 helices and 0 strands
  ################################################################################
  Comparing model.pdb to reference.pdb
  Chains removed from reference.pdb: _
  Chains in reference.pdb: AB
  Chains in model.pdb: AB
  Chemically equivalent chain-groups in reference.pdb: [['A', 'B']]
  Chemically equivalent chain-groups in model.pdb: [['A', 'B']]
  Chemical chain-groups mapping: {('A', 'B'): ('A', 'B')}
  Identifying Symmetry Groups...
  Symmetry threshold 0.1 used for angles of reference.pdb
  Symmetry threshold 0.1 used for axis of reference.pdb
  Symmetry threshold 0.1 used for angles of model.pdb
  Symmetry threshold 0.1 used for axis of model.pdb
  Selecting Symmetry Groups...
  Symmetry-groups used in reference.pdb: [('A',), ('B',)]
  Symmetry-groups used in model.pdb: [('A',), ('B',)]
  Closed Symmetry with strict parameters
  Mapping found: {'A': 'A', 'B': 'B'}
  --------------------------------------------------------------------------------
  Checking consistency between model.pdb and reference.pdb
  Consistency check: OK
  --------------------------------------------------------------------------------
  Computing QS-score
  QSscore reference.pdb, model.pdb: best: 0.96, global: 0.96
  ################################################################################
  Saving output into output_qs.json
