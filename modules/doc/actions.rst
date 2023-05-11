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
``ost compare-structures`` action. This can be considered a command line
interface to :class:`ost.mol.alg.scoring.Scorer`

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
                                [--local-lddt] [--bb-lddt] [--bb-local-lddt]
                                [--cad-score] [--local-cad-score]
                                [--cad-exec CAD_EXEC] [--qs-score]
                                [--rigid-scores] [--interface-scores]
                                [--patch-scores]
  
  Evaluate model against reference 
  
  Example: ost compare-structures -m model.pdb -r reference.cif
  
  Loads the structures and performs basic cleanup:
  
   * Assign elements according to the PDB Chemical Component Dictionary
   * Map nonstandard residues to their parent residues as defined by the PDB
     Chemical Component Dictionary, e.g. phospho-serine => serine
   * Remove hydrogens
   * Remove OXT atoms
   * Remove unknown atoms, i.e. atoms that are not expected according to the PDB
     Chemical Component Dictionary
   * Select for peptide/nucleotide residues
  
  The cleaned structures are optionally dumped using -d/--dump-structures
  
  Output is written in JSON format (default: out.json). In case of no additional
  options, this is a dictionary with 8 keys:

   * "reference_chains": Chain names of reference
   * "model_chains": Chain names of model
   * "chem_groups": Groups of polypeptides/polynucleotides from reference that
     are considered chemically equivalent. You can derive stoichiometry from this.
     Contains only chains that are considered in chain mapping, i.e. pass a
     size threshold (defaults: 10 for peptides, 4 for nucleotides).
   * "chem_mapping": List of same length as "chem_groups". Assigns model chains to
     the respective chem group. Again, only contains chains that are considered
     in chain mapping.
   * "chain_mapping": A dictionary with reference chain names as keys and the
     mapped model chain names as values. Missing chains are either not mapped
     (but present in "chem_groups", "chem_mapping") or were not considered in
     chain mapping (short peptides etc.)
   * "aln": Pairwise sequence alignment for each pair of mapped chains in fasta
     format.
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

  Example to compute global and per-residue lDDT values as well as QS-score:

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
                          default, the asymmetric unit (AU) is used for scoring.
                          If there are biounits defined in the mmcif file, you
                          can specify the (0-based) index of the one which
                          should be used.
    -rb REFERENCE_BIOUNIT, --reference-biounit REFERENCE_BIOUNIT
                          Only has an effect if reference is in mmcif format. By
                          default, the asymmetric unit (AU) is used for scoring.
                          If there are biounits defined in the mmcif file, you
                          can specify the (0-based) index of the one which
                          should be used.
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
    --bb-lddt             Compute global lDDT score with default
                          parameterization and store as key "bb_lddt". lDDT in
                          this case is only computed on backbone atoms: CA for
                          peptides and C3' for nucleotides
    --bb-local-lddt       Compute per-residue lDDT scores with default
                          parameterization and store as key "bb_local_lddt".
                          lDDT in this case is only computed on backbone atoms:
                          CA for peptides and C3' for nucleotides. Score for
                          model residue with number 42 in chain X can be
                          extracted with: data["local_lddt"]["X"]["42"]. If
                          there is an insertion code, lets say A, the last key
                          becomes "42A"
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
    --tm-score            Computes TM-score with the USalign tool. Also computes
                          a chain mapping in case of complexes that is stored
                          in the same format as the default mapping. TM-score
                          and the mapping are available as keys "tm_score" and
                          "usalign_mapping"



.. _ost compare ligand structures:

Comparing two structures with ligands
--------------------------------------------------------------------------------

You can compare two structures with non-polymer/small molecule ligands and
compute lDDT-PLI and ligand RMSD scores from the command line with the
``ost compare-ligand-structures`` action. This can be considered a command
line interface to :class:`ost.mol.alg.ligand_scoring.LigandScorer`.

Details on the usage (output of ``ost compare-ligand-structures --help``):

.. code-block:: console

    usage: ost compare-ligand-structures [-h] -m MODEL [-ml [MODEL_LIGANDS ...]]
                                     -r REFERENCE
                                     [-rl [REFERENCE_LIGANDS ...]]
                                     [-o OUTPUT] [-mf {pdb,mmcif,cif}]
                                     [-rf {pdb,mmcif,cif}] [-ft] [-rna] [-ec]
                                     [-sm] [--lddt-pli] [--rmsd]
                                     [--radius RADIUS]
                                     [--lddt-pli-radius LDDT_PLI_RADIUS]
                                     [--lddt-bs-radius LDDT_BS_RADIUS]
                                     [-v VERBOSITY]

    Evaluate model with non-polymer/small molecule ligands against reference.

    Example: ost compare-ligand-structures \
        -m model.pdb \
        -ml ligand.sdf \
        -r reference.cif \
        --lddt-pli --rmsd

    Structures of polymer entities (proteins and nucleotides) can be given in PDB
    or mmCIF format. If the structure is given in mmCIF format, only the asymmetric
    unit (AU) is used for scoring.

    Ligands can be given as path to SDF files containing the ligand for both model
    (--model-ligands/-ml) and reference (--reference-ligands/-rl). If omitted,
    ligands will be detected in the model and reference structures. For structures
    given in mmCIF format, this is based on the annotation as "non polymer entity"
    (i.e. ligands in the _pdbx_entity_nonpoly mmCIF category) and works reliably.
    For structures given in PDB format, this is based on the HET records and is
    normally not what you want. You should always give ligands as SDF for
    structures in PDB format.

    Polymer/oligomeric ligands (saccharides, peptides, nucleotides) are not
    supported.

    Only minimal cleanup steps are performed (remove hydrogens, and for structures
    of polymers only, remove unknown atoms and cleanup element column).

    Ligands in mmCIF and PDB files must comply with the PDB component dictionary
    definition, and have properly named residues and atoms, in order for
    ligand connectivity to be loaded correctly. Ligands loaded from SDF files
    are exempt from this restriction, meaning any arbitrary ligand can be assessed.

    Output is written in JSON format (default: out.json). In case of no additional
    options, this is a dictionary with three keys:

     * "model_ligands": A list of ligands in the model. If ligands were provided
       explicitly with --model-ligands, elements of the list will be the paths to
       the ligand SDF file(s). Otherwise, they will be the chain name and residue
       number of the ligand, separated by a dot.
     * "reference_ligands": A list of ligands in the reference. If ligands were
       provided explicitly with --reference-ligands, elements of the list will be
       the paths to the ligand SDF file(s). Otherwise, they will be the chain name
       and residue number of the ligand, separated by a dot.
     * "status": SUCCESS if everything ran through. In case of failure, the only
       content of the JSON output will be "status" set to FAILURE and an
       additional key: "traceback".

    Each score is opt-in and, be enabled with optional arguments and is added
    to the output. Keys correspond to the values in "model_ligands" above.
    Only assigned mapped ligands are reported.

    options:
      -h, --help            show this help message and exit
      -m MODEL, --mdl MODEL, --model MODEL
                            Path to model file.
      -ml [MODEL_LIGANDS ...], --mdl-ligands [MODEL_LIGANDS ...],
                            --model-ligands [MODEL_LIGANDS ...]
                            Path to model ligand files.
      -r REFERENCE, --ref REFERENCE, --reference REFERENCE
                            Path to reference file.
      -rl [REFERENCE_LIGANDS ...], --ref-ligands [REFERENCE_LIGANDS ...],
                            --reference-ligands [REFERENCE_LIGANDS ...]
                            Path to reference ligand files.
      -o OUTPUT, --out OUTPUT, --output OUTPUT
                            Output file name. The output will be saved as a JSON
                            file. default: out.json
      -mf {pdb,mmcif,cif}, --mdl-format {pdb,mmcif,cif},
                            --model-format {pdb,mmcif,cif}
                            Format of model file. Inferred from path if not
                            given.
      -rf {pdb,mmcif,cif}, --reference-format {pdb,mmcif,cif},
                            --ref-format {pdb,mmcif,cif}
                            Format of reference file. Inferred from path if not
                            given.
      -ft, --fault-tolerant
                            Fault tolerant parsing.
      -rna, --residue-number-alignment
                            Make alignment based on residue number instead of
                            using a global BLOSUM62-based alignment (NUC44 for
                            nucleotides).
      -ec, --enforce-consistency
                            Enforce consistency of residue names between the
                            reference binding site and the model. By default
                            residue name discrepancies are reported but the
                            program proceeds. If this is set to True, the program
                            will fail with an error message if the residues names
                            differ. Note: more binding site mappings may be
                            explored during scoring, but only inconsistencies in
                            the selected mapping are reported.
      -sm, --substructure-match
                            Allow incomplete target ligands.
      --lddt-pli            Compute lDDT-PLI score and store as key "lddt-pli".
      --rmsd                Compute RMSD score and store as key "rmsd".
      --radius RADIUS       Inclusion radius for the binding site. Any residue
                            with atoms within this distance of the ligand will be
                            included in the binding site.
      --lddt-pli-radius LDDT_PLI_RADIUS
                            lDDT inclusion radius for lDDT-PLI.
      --lddt-bs-radius LDDT_BS_RADIUS
                            lDDT inclusion radius for lDDT-BS.
      -v VERBOSITY, --verbosity VERBOSITY
                            Set verbosity level. Defaults to 3 (INFO).


Additional information about the scores and output values is available in
:meth:`rmsd_details <ost.mol.alg.ligand_scoring.LigandScorer.rmsd_details>` and
:meth:`lddt_pli_details <ost.mol.alg.ligand_scoring.LigandScorer.lddt_pli_details>`.