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

You can compare two structures from the command line with the
``ost compare-structures`` action.

.. warning::

  ``compare-structures`` underwent a complete rewrite in OpenStructure
  release 2.4.0. The old version is still available as
  ``compare-structures-legacy`` with documentation available
  :doc:`here <deprecated_actions>`.

Details on the usage (output of ``ost compare-structures --help``):

.. code-block:: console

  usage: ost compare-structures [-h] -m MODEL -r REFERENCE [-o OUTPUT]
                                [-mf {pdb,cif,mmcif}] [-rf {pdb,cif,mmcif}]
                                [-mb MODEL_BIOUNIT] [-rb REFERENCE_BIOUNIT]
                                [-rna] [-ec] [-d] [-ds DUMP_SUFFIX] [-ft]
                                [-c CHAIN_MAPPING [CHAIN_MAPPING ...]] [--lddt]
                                [--local-lddt] [--cad-score] [--local-cad-score]
                                [--cad-exec CAD_EXEC] [--qs-score]
                                [--rigid-scores] [--interface-scores]
                                [--patch-scores]
  
  Evaluate model against reference 
  
  Example: ost compare-structures -m model.pdb -r reference.cif
  
  Loads the structures and performs basic cleanup:
  
   * Assign elements according to the PDB compound dictionary
   * Map nonstandard residues to their parent residues as defined by the PDB
     compound dictionary, e.g. phospho-serine => serine
   * Remove hydrogens
   * Remove OXT atoms
   * Remove unknown atoms, i.e. atoms that are not expected according to the PDB
     compound dictionary
  
  The cleaned structures are optionally dumped using -d/--dump-structures
  
  Output is written in JSON format (default: out.json). In case of no additional
  options, this is a dictionary with five keys:
  
   * "chain_mapping": A dictionary with reference chain names as keys and the
     mapped model chain names as values.
   * "aln": Pairwise sequence alignment for each pair of mapped chains in fasta
     format.
   * "chem_groups": Groups of polypeptides/polynucleotides that are considered
     chemically equivalent. You can derive stoichiometry from this.
   * "inconsistent_residues": List of strings that represent name mismatches of
     aligned residues in form
     <trg_cname>.<trg_rname><trg_rnum>-<mdl_cname>.<mdl_rname><mdl_rnum>.
     Inconsistencies may lead to corrupt results but do not abort the program.
     Program abortion in these cases can be enforced with
     -ec/--enforce-consistency.
   * "status": SUCCESS if everything ran through. In case of failure, the only
     content of the JSON output will be "status" set to FAILURE and an
     additional key: "traceback".
  
  The pairwise sequence alignments are computed with Needleman-Wunsch using
  BLOSUM62 (NUC44 for nucleotides). Many benchmarking scenarios preprocess the
  structures to ensure matching residue numbers (CASP/CAMEO). In these cases,
  enabling -rna/--residue-number-alignment is recommended.
  
  Each score is opt-in and can be enabled with optional arguments.
  
  Example to compute local and per-residue lDDT values as well as QS-score:
  
  ost compare-structures -m model.pdb -r reference.cif --lddt --local-lddt --qs-score
  
  Example to inject custom chain mapping
  
  ost compare-structures -m model.pdb -r reference.cif -c A:B B:A

  optional arguments:
    -h, --help            show this help message and exit
    -m MODEL, --model MODEL
                          Path to model file.
    -r REFERENCE, --reference REFERENCE
                          Path to reference file.
    -o OUTPUT, --output OUTPUT
                          Output file name. The output will be saved as a JSON
                          file. default: out.json
    -mf {pdb,cif,mmcif}, --model-format {pdb,cif,mmcif}
                          Format of model file. pdb reads pdb but also pdb.gz,
                          same applies to cif/mmcif. Inferred from filepath if
                          not given.
    -rf {pdb,cif,mmcif}, --reference-format {pdb,cif,mmcif}
                          Format of reference file. pdb reads pdb but also
                          pdb.gz, same applies to cif/mmcif. Inferred from
                          filepath if not given.
    -mb MODEL_BIOUNIT, --model-biounit MODEL_BIOUNIT
                          Only has an effect if model is in mmcif format. By
                          default, the assymetric unit (AU) is used for scoring.
                          If there are biounits defined in the mmcif file, you
                          can specify the index of the one which should be used.
    -rb REFERENCE_BIOUNIT, --reference-biounit REFERENCE_BIOUNIT
                          Only has an effect if reference is in mmcif format. By
                          default, the assymetric unit (AU) is used for scoring.
                          If there are biounits defined in the mmcif file, you
                          can specify the index of the one which should be used.
    -rna, --residue-number-alignment
                          Make alignment based on residue number instead of
                          using a global BLOSUM62-based alignment (NUC44 for
                          nucleotides).
    -ec, --enforce-consistency
                          Enforce consistency. By default residue name
                          discrepancies between a model and reference are
                          reported but the program proceeds. If this flag is ON,
                          the program fails for these cases.
    -d, --dump-structures
                          Dump cleaned structures used to calculate all the
                          scores as PDB files using specified suffix. Files will
                          be dumped to the same location as original files.
    -ds DUMP_SUFFIX, --dump-suffix DUMP_SUFFIX
                          Use this suffix to dump structures. Defaults to
                          .compare.structures.pdb.
    -ft, --fault-tolerant
                          Fault tolerant parsing.
    -c CHAIN_MAPPING [CHAIN_MAPPING ...], --chain-mapping CHAIN_MAPPING [CHAIN_MAPPING ...]
                          Custom mapping of chains between the reference and the
                          model. Each separate mapping consist of key:value
                          pairs where key is the chain name in reference and
                          value is the chain name in model.
    --lddt                Compute global lDDT score with default
                          parameterization and store as key "lddt".
                          Stereochemical irregularities affecting lDDT are
                          reported as keys "model_clashes", "model_bad_bonds",
                          "model_bad_angles" and the respective reference
                          counterparts.
    --local-lddt          Compute per-residue lDDT scores with default
                          parameterization and store as key "local_lddt". Score
                          for model residue with number 42 in chain X can be
                          extracted with: data["local_lddt"]["X"]["42"]. If
                          there is an insertion code, lets say A, the last key
                          becomes "42A" Stereochemical irregularities affecting
                          lDDT are reported as keys "model_clashes",
                          "model_bad_bonds", "model_bad_angles" and the
                          respective reference counterparts.
    --cad-score           Compute global CAD's atom-atom (AA) score and store as
                          key "cad_score". --residue-number-alignment must be
                          enabled to compute this score. Requires
                          voronota_cadscore executable in PATH. Alternatively
                          you can set cad-exec.
    --local-cad-score     Compute local CAD's atom-atom (AA) scores and store as
                          key "local_cad_score". Score for model residue with
                          number 42 in chain X can be extracted with:
                          data["local_cad_score"]["X"]["42"]. --residue-number-
                          alignments must be enabled to compute this score.
                          Requires voronota_cadscore executable in PATH.
                          Alternatively you can set cad-exec.
    --cad-exec CAD_EXEC   Path to voronota-cadscore executable (installed from
                          https://github.com/kliment-olechnovic/voronota).
                          Searches PATH if not set.
    --qs-score            Compute QS-score, stored as key "qs_global", and the
                          QS-best variant, stored as key "qs_best".
    --rigid-scores        Computes rigid superposition based scores. They're
                          based on a Kabsch superposition of all mapped CA
                          positions (C3' for nucleotides). Makes the following
                          keys available: "oligo_gdtts": GDT with distance
                          thresholds [1.0, 2.0, 4.0, 8.0] given these positions
                          and transformation, "oligo_gdtha": same with
                          thresholds [0.5, 1.0, 2.0, 4.0], "rmsd": RMSD given
                          these positions and transformation, "transform": the
                          used 4x4 transformation matrix that superposes model
                          onto reference.
    --interface-scores    Per interface scores for each interface that has at
                          least one contact in the reference, i.e. at least one
                          pair of heavy atoms within 5A. The respective
                          interfaces are available from key "interfaces" which
                          is a list of tuples in form (ref_ch1, ref_ch2,
                          mdl_ch1, mdl_ch2). Per-interface scores are available
                          as lists referring to these interfaces and have the
                          following keys: "nnat" (number of contacts in
                          reference), "nmdl" (number of contacts in model),
                          "fnat" (fraction of reference contacts which are also
                          there in model), "fnonnat" (fraction of model contacts
                          which are not there in target), "irmsd" (interface
                          RMSD), "lrmsd" (ligand RMSD), "dockq_scores" (per-
                          interface score computed from "fnat", "irmsd" and
                          "lrmsd"), "interface_qs_global" and
                          "interface_qs_best" (per-interface versions of the two
                          QS-score variants). The DockQ score is strictly
                          designed to score each interface individually. We also
                          provide two averaged versions to get one full model
                          score: "dockq_ave", "dockq_wave". The first is simply
                          the average of "dockq_scores", the latter is a
                          weighted average with weights derived from "nnat".
                          These two scores only consider interfaces that are
                          present in both, the model and the reference.
                          "dockq_ave_full" and "dockq_wave_full" add zeros in
                          the average computation for each interface that is
                          only present in the reference but not in the model.
    --patch-scores        Local interface quality score used in CASP15. Scores
                          each model residue that is considered in the interface
                          (CB pos within 8A of any CB pos from another chain (CA
                          for GLY)). The local neighborhood gets represented by
                          "interface patches" which are scored with QS-score and
                          DockQ. Scores where not the full patches are
                          represented by the reference are set to None. Model
                          interface residues are available as key
                          "model_interface_residues", reference interface
                          residues as key "reference_interface_residues".
                          Residues are represented as string in form
                          <num><inscode>. The respective scores are available as
                          keys "patch_qs" and "patch_dockq"
