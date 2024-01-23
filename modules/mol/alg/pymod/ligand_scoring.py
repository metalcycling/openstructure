import warnings

import numpy as np
import numpy.ma as np_ma
import networkx

from ost import mol
from ost import geom
from ost import LogError, LogWarning, LogScript, LogInfo, LogVerbose, LogDebug
from ost.mol.alg import chain_mapping


class LigandScorer:
    """ Scorer class to compute various small molecule ligand (non polymer) scores.

    .. note ::
      Extra requirements:

      - Python modules `numpy` and `networkx` must be available
        (e.g. use ``pip install numpy networkx``)

    At the moment, two scores are available:

    * lDDT-PLI, that looks at the conservation of protein-ligand contacts
      with :class:`lDDT <ost.mol.alg.lddt.lDDTScorer>`.
    * Binding-site superposed, symmetry-corrected RMSD that assesses the
      accuracy of the ligand pose (BiSyRMSD, hereinafter referred to as RMSD).

    Both scores involve local chain mapping of the reference binding site
    onto the model, symmetry-correction, and finally assignment (mapping)
    of model and target ligands, as described in (Manuscript in preparation).

    The binding site is defined based on a radius around the target ligand
    in the reference structure only. It only contains protein and nucleic
    acid chains that pass the criteria for the
    :class:`chain mapping <ost.mol.alg.chain_mapping>`. This means ignoring
    other ligands, waters, short polymers as well as any incorrectly connected
    chains that may be in proximity.

    Results are available as matrices (`(lddt_pli|rmsd)_matrix`), where every
    target-model score is reported in a matrix; as `(lddt_pli|rmsd)` where
    a model-target assignment has been determined (see below) and reported in
    a dictionary; and as (`(lddt_pli|rmsd)_details`) methods, which report
    additional details about different aspects of the scoring such as chain
    mapping.

    The behavior of chain mapping and ligand assignment can be controlled
    with the `global_chain_mapping` and `rmsd_assignment` arguments.

    By default, chain mapping is performed locally, ie. only within the
    binding site. As a result, different ligand scores can correspond to
    different chain mappings. This tends to produce more favorable scores,
    especially in large, partially regular oligomeric complexes.
    Setting `global_chain_mapping=True` enforces a single global chain mapping,
    as per :meth:`ost.mol.alg.chain_mapping.ChainMapper.GetMapping`.
    Note that this global chain mapping currently ignores non polymer entities
    such as small ligands, and may result in overly pessimistic scores.

    By default, target-model ligand assignments are computed independently
    for the RMSD and lDDT-PLI scores. For RMSD, each model ligand is uniquely
    assigned to a target ligand, starting from the "best" possible mapping
    (lowest RMSD) and using each target and model ligand in a single
    assignment. Ties are resolved by best (highest) lDDT-PLI. Similarly,
    for lDDT-PLI, the assignment is based on the highest lDDT-PLI, and ties
    broken by lowest RMSD. Setting `rmsd_assignment=True` forces a single
    ligand assignment, based on RMSD only. Ties are broken arbitrarily.

    By default, only exact matches between target and model ligands are
    considered. This is a problem when the target only contains a subset
    of the expected atoms (for instance if atoms are missing in an
    experimental structure, which often happens in the PDB). With
    `substructure_match=True`, complete model ligands can be scored against
    partial target ligands. One problem with this approach is that it is
    very easy to find good matches to small, irrelevant ligands like EDO, CO2
    or GOL. To counter that, the assignment algorithm considers the coverage,
    expressed as the fraction of atoms of the model ligand atoms covered in the
    target. Higher coverage matches are prioritized, but a match with a better
    score will be preferred if it falls within a window of `coverage_delta`
    (by default 0.2) of a worse-scoring match. As a result, for instance,
    with a delta of 0.2, a low-score match with coverage 0.96 would be
    preferred to a high-score match with coverage 0.90.

    Assumptions:

    The class generally assumes that the
    :attr:`~ost.mol.ResidueHandle.is_ligand` property is properly set on all
    the ligand atoms, and only ligand atoms. This is typically the case for
    entities loaded from mmCIF (tested with mmCIF files from the PDB and
    SWISS-MODEL), but it will most likely not work for most entities loaded
    from PDB files.

    The class doesn't perform any cleanup of the provided structures.
    It is up to the caller to ensure that the data is clean and suitable for
    scoring. :ref:`Molck <molck>` should be used with extra
    care, as many of the options (such as `rm_non_std` or `map_nonstd_res`) can
    cause ligands to be removed from the structure. If cleanup with Molck is
    needed, ligands should be kept aside and passed separately. Non-ligand residues
    should be valid compounds with atom names following the naming conventions
    of the component dictionary. Non-standard residues are acceptable, and if
    the model contains a standard residue at that position, only atoms with
    matching names will be considered.

    Unlike most of OpenStructure, this class does not assume that the ligands
    (either in the model or the target) are part of the PDB component
    dictionary. They may have arbitrary residue names. Residue names do not
    have to match between the model and the target. Matching is based on
    the calculation of isomorphisms which depend on the atom element name and
    atom connectivity (bond order is ignored).
    It is up to the caller to ensure that the connectivity of atoms is properly
    set before passing any ligands to this class. Ligands with improper
    connectivity will lead to bogus results.

    Note, however, that atom names should be unique within a residue (ie two
    distinct atoms cannot have the same atom name).

    This only applies to the ligand. The rest of the model and target
    structures (protein, nucleic acids) must still follow the usual rules and
    contain only residues from the compound library.

    Although it isn't a requirement, hydrogen atoms should be removed from the
    structures. Here is an example code snippet that will perform a reasonable
    cleanup. Keep in mind that this is most likely not going to work as
    expected with entities loaded from PDB files, as the `is_ligand` flag is
    probably not set properly.

    Here is a snippet example of how to use this code::

        from ost.mol.alg.ligand_scoring import LigandScorer
        from ost.mol.alg import Molck, MolckSettings

        # Load data
        # Structure model in PDB format, containing the receptor only
        model = io.LoadPDB("path_to_model.pdb")
        # Ligand model as SDF file
        model_ligand = io.LoadEntity("path_to_ligand.sdf", format="sdf")
        # Target loaded from mmCIF, containing the ligand
        target = io.LoadMMCIF("path_to_target.cif")

        # Cleanup a copy of the structures
        cleaned_model = model.Copy()
        cleaned_target = target.Copy()
        molck_settings = MolckSettings(rm_unk_atoms=True,
                                       rm_non_std=False,
                                       rm_hyd_atoms=True,
                                       rm_oxt_atoms=False,
                                       rm_zero_occ_atoms=False,
                                       colored=False,
                                       map_nonstd_res=False,
                                       assign_elem=True)
        Molck(cleaned_model, conop.GetDefaultLib(), molck_settings)
        Molck(cleaned_target, conop.GetDefaultLib(), molck_settings)

        # Setup scorer object and compute lDDT-PLI
        model_ligands = [model_ligand.Select("ele != H")]
        ls = LigandScorer(model=cleaned_model, target=cleaned_target, model_ligands=model_ligands)
        print("lDDT-PLI:", ls.lddt_pli)
        print("RMSD:", ls.rmsd)

    :param model: Model structure - a deep copy is available as :attr:`model`.
                  No additional processing (ie. Molck), checks,
                  stereochemistry checks or sanitization is performed on the
                  input. Hydrogen atoms are kept.
    :type model: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param target: Target structure - a deep copy is available as :attr:`target`.
                  No additional processing (ie. Molck), checks or sanitization
                  is performed on the input. Hydrogen atoms are kept.
    :type target: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param model_ligands: Model ligands, as a list of
                  :class:`~ost.mol.ResidueHandle` belonging to the model
                  entity. Can be instantiated with either a :class:list of
                  :class:`~ost.mol.ResidueHandle`/:class:`ost.mol.ResidueView`
                  or of :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`.
                  If `None`, ligands will be extracted from the `model` entity,
                  from chains with :class:`~ost.mol.ChainType`
                  `CHAINTYPE_NON_POLY` (this is normally set properly in
                  entities loaded from mmCIF).
    :type model_ligands: :class:`list`
    :param target_ligands: Target ligands, as a list of
                  :class:`~ost.mol.ResidueHandle` belonging to the target
                  entity. Can be instantiated either a :class:list of
                  :class:`~ost.mol.ResidueHandle`/:class:`ost.mol.ResidueView`
                  or of :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
                  containing a single residue each. If `None`, ligands will be
                  extracted from the `target` entity, from chains with
                  :class:`~ost.mol.ChainType` `CHAINTYPE_NON_POLY` (this is
                  normally set properly in entities loaded from mmCIF).
    :type target_ligands: :class:`list`
    :param resnum_alignments: Whether alignments between chemically equivalent
                              chains in *model* and *target* can be computed
                              based on residue numbers. This can be assumed in
                              benchmarking setups such as CAMEO/CASP.
    :type resnum_alignments: :class:`bool`
    :param check_resnames:  On by default. Enforces residue name matches
                            between mapped model and target residues.
    :type check_resnames: :class:`bool`
    :param rename_ligand_chain: If a residue with the same chain name and
                                residue number than an explicitly passed model
                                or target ligand exits in the structure,
                                and `rename_ligand_chain` is False, a
                                RuntimeError will be raised. If
                                `rename_ligand_chain` is True, the ligand will
                                be moved to a new chain instead, and the move
                                will be logged to the console with SCRIPT
                                level.
    :type rename_ligand_chain: :class:`bool`
    :param chain_mapper: a chain mapper initialized for the target structure.
                         If None (default), a chain mapper will be initialized
                         lazily as required.
    :type chain_mapper:  :class:`ost.mol.alg.chain_mapping.ChainMapper`
    :param substructure_match: Set this to True to allow partial target ligand.
    :type substructure_match: :class:`bool`
    :param coverage_delta: the coverage delta for partial ligand assignment.
    :type coverage_delta: :class:`float`
    :param radius: Inclusion radius for the binding site. Residues with
                   atoms within this distance of the ligand will be considered
                   for inclusion in the binding site.
    :type radius: :class:`float`
    :param lddt_pli_radius: lDDT inclusion radius for lDDT-PLI.
    :type lddt_pli_radius: :class:`float`
    :param lddt_lp_radius: lDDT inclusion radius for lDDT-LP.
    :type lddt_lp_radius: :class:`float`
    :param binding_sites_topn: maximum number of target binding site
                               representations to assess, per target ligand.
                               Ignored if `global_chain_mapping` is True.
    :type binding_sites_topn: :class:`int`
    :param global_chain_mapping: set to True to use a global chain mapping for
                                 the polymer (protein, nucleotide) chains.
                                 Defaults to False, in which case only local
                                 chain mappings are allowed (where different
                                 ligand may be scored against different chain
                                 mappings).
    :type global_chain_mapping: :class:`bool`
    :param custom_mapping: Provide custom chain mapping between *model* and
                           *target* that is used as global chain mapping.
                           Dictionary with target chain names as key and model
                           chain names as value. Only has an effect if
                           *global_chain_mapping* is True.
    :type custom_mapping: :class:`dict`
    :param rmsd_assignment: assign ligands based on RMSD only. The default
                            (False) is to use a combination of lDDT-PLI and
                            RMSD for the assignment.
    :type rmsd_assignment: :class:`bool`
    :param n_max_naive: Parameter for global chain mapping. If *model* and
                        *target* have less or equal that number of chains,
                        the full
                        mapping solution space is enumerated to find the
                        the optimum. A heuristic is used otherwise.
    :type n_max_naive: :class:`int`
    :param max_symmetries: If more than that many isomorphisms exist for
                       a target-ligand pair, it will be ignored and reported
                       as unassigned.
    :type max_symmetries: :class:`int`
    :param unassigned: If True, unassigned model ligands are reported in
                       the output together with assigned ligands, with a score
                       of None, and reason for not being assigned in the
                       \\*_details matrix. Defaults to False.
    :type unassigned: :class:`bool`
    """
    def __init__(self, model, target, model_ligands=None, target_ligands=None,
                 resnum_alignments=False, check_resnames=True,
                 rename_ligand_chain=False,
                 chain_mapper=None, substructure_match=False,
                 coverage_delta=0.2,
                 radius=4.0, lddt_pli_radius=6.0, lddt_lp_radius=10.0,
                 binding_sites_topn=100000, global_chain_mapping=False,
                 rmsd_assignment=False, n_max_naive=12, max_symmetries=1e5,
                 custom_mapping=None, unassigned=False):

        if isinstance(model, mol.EntityView):
            self.model = mol.CreateEntityFromView(model, False)
        elif isinstance(model, mol.EntityHandle):
            self.model = model.Copy()
        else:
            raise RuntimeError("model must be of type EntityView/EntityHandle")

        if isinstance(target, mol.EntityView):
            self.target = mol.CreateEntityFromView(target, False)
        elif isinstance(target, mol.EntityHandle):
            self.target = target.Copy()
        else:
            raise RuntimeError("target must be of type EntityView/EntityHandle")

        # Extract ligands from target
        if target_ligands is None:
            self.target_ligands = self._extract_ligands(self.target)
        else:
            self.target_ligands = self._prepare_ligands(self.target, target,
                                                        target_ligands,
                                                        rename_ligand_chain)
        if len(self.target_ligands) == 0:
            LogWarning("No ligands in the target")

        # Extract ligands from model
        if model_ligands is None:
            self.model_ligands = self._extract_ligands(self.model)
        else:
            self.model_ligands = self._prepare_ligands(self.model, model,
                                                       model_ligands,
                                                       rename_ligand_chain)
        if len(self.model_ligands) == 0:
            LogWarning("No ligands in the model")
            if len(self.target_ligands) == 0:
                raise ValueError("No ligand in the model and in the target")

        self._chain_mapper = chain_mapper
        self.resnum_alignments = resnum_alignments
        self.check_resnames = check_resnames
        self.rename_ligand_chain = rename_ligand_chain
        self.substructure_match = substructure_match
        self.radius = radius
        self.lddt_pli_radius = lddt_pli_radius
        self.lddt_lp_radius = lddt_lp_radius
        self.binding_sites_topn = binding_sites_topn
        self.global_chain_mapping = global_chain_mapping
        self.rmsd_assignment = rmsd_assignment
        self.n_max_naive = n_max_naive
        self.max_symmetries = max_symmetries
        self.unassigned = unassigned
        self.coverage_delta = coverage_delta

        # scoring matrices
        self._rmsd_matrix = None
        self._rmsd_full_matrix = None
        self._lddt_pli_matrix = None
        self._lddt_pli_full_matrix = None

        # lazily computed scores
        self._rmsd = None
        self._rmsd_details = None
        self._lddt_pli = None
        self._lddt_pli_details = None

        # lazily precomputed variables
        self._binding_sites = {}
        self.__model_mapping = None

        # Bookkeeping of unassigned ligands
        self._unassigned_target_ligands = None
        self._unassigned_model_ligands = None
        self._unassigned_target_ligands_reason = {}
        self._unassigned_target_ligand_short = None
        self._unassigned_model_ligand_short = None
        self._unassigned_target_ligand_descriptions = None
        self._unassigned_model_ligand_descriptions = None
        # Keep track of symmetries/isomorphisms (regardless of scoring)
        # 0.0: no isomorphism
        # 1.0: isomorphic
        # np.nan: not assessed yet - that's why we can't use a boolean
        self._assignment_isomorphisms = None
        # Keep track of match coverage (only in case there was a score)
        self._assignment_match_coverage = None

        if custom_mapping is not None:
            self._set_custom_mapping(custom_mapping)

    @property
    def chain_mapper(self):
        """ Chain mapper object for the given :attr:`target`.

        :type: :class:`ost.mol.alg.chain_mapping.ChainMapper`
        """
        if self._chain_mapper is None:
            self._chain_mapper = chain_mapping.ChainMapper(self.target,
                                                           n_max_naive=1e9,
                                                           resnum_alignments=self.resnum_alignments)
        return self._chain_mapper

    @property
    def _model_mapping(self):
        """Get the global chain mapping for the model."""
        if self.__model_mapping is None:
            self.__model_mapping = self.chain_mapper.GetMapping(self.model,
                                                                n_max_naive=self.n_max_naive)
        return self.__model_mapping

    @staticmethod
    def _extract_ligands(entity):
        """Extract ligands from entity. Return a list of residues.

        Assumes that ligands are contained in one or more chain with chain type
        `mol.ChainType.CHAINTYPE_NON_POLY`. This is typically the case
        for entities loaded from mmCIF (tested with mmCIF files from the PDB
        and SWISS-MODEL), but it will most likely not work for most entities
        loaded from PDB files.

        As a deviation from the mmCIF semantics, we allow a chain, set as
        `CHAINTYPE_NON_POLY`, to contain more than one ligand. This function
        performs basic checks to ensure that the residues in this chain are
        not forming polymer bonds (ie peptide/nucleotide ligands) and will
        raise a RuntimeError if this assumption is broken.

        Note: This will not extract ligands based on the HET record in the old
        PDB style, as this is not a reliable indicator and depends on how the
        entity was loaded.

        :param entity: the entity to extract ligands from
        :type entity: :class:`~ost.mol.EntityHandle`
        :rtype: :class:`list` of :class:`~ost.mol.ResidueHandle`

        """
        extracted_ligands = []
        for chain in entity.chains:
            if chain.chain_type == mol.ChainType.CHAINTYPE_NON_POLY:
                for residue in chain.residues:
                    if mol.InSequence(residue, residue.next):
                        raise RuntimeError("Connected residues in non polymer "
                                           "chain %s" % (chain.name))
                    residue.SetIsLigand(True)  # just in case
                    extracted_ligands.append(residue)
                    LogVerbose("Detected residue %s as ligand" % residue)
        return extracted_ligands

    @staticmethod
    def _prepare_ligands(new_entity, old_entity, ligands, rename_chain):
        """Prepare the ligands given into a list of ResidueHandles which are
        part of the copied entity, suitable for the model_ligands and
        target_ligands properties.

        This function takes a list of ligands as (Entity|Residue)(Handle|View).
        Entities can contain multiple ligands, which will be considered as
        separate ligands.

        Ligands which are part of the entity are simply fetched in the new
        copied entity. Otherwise, they are copied over to the copied entity.
        """
        extracted_ligands = []

        next_chain_num = 1
        new_editor = None

        def _copy_residue(residue, rename_chain):
            """ Copy the residue into the new chain.
            Return the new residue handle."""
            nonlocal next_chain_num, new_editor

            # Instantiate the editor
            if new_editor is None:
                new_editor = new_entity.EditXCS()

            new_chain = new_entity.FindChain(residue.chain.name)
            if not new_chain.IsValid():
                new_chain = new_editor.InsertChain(residue.chain.name)
            else:
                # Does a residue with the same name already exist?
                already_exists = new_chain.FindResidue(residue.number).IsValid()
                if already_exists:
                    if rename_chain:
                        chain_ext = 2  # Extend the chain name by this
                        while True:
                            new_chain_name = residue.chain.name + "_" + str(chain_ext)
                            new_chain = new_entity.FindChain(new_chain_name)
                            if new_chain.IsValid():
                                chain_ext += 1
                                continue
                            else:
                                new_chain = new_editor.InsertChain(new_chain_name)
                                break
                        LogScript("Moved ligand residue %s to new chain %s" % (
                            residue.qualified_name, new_chain.name))
                    else:
                        msg = "A residue number %s already exists in chain %s" % (
                            residue.number, residue.chain.name)
                        raise RuntimeError(msg)

            # Add the residue with its original residue number
            new_res = new_editor.AppendResidue(new_chain, residue.name, residue.number)
            # Add atoms
            for old_atom in residue.atoms:
                new_editor.InsertAtom(new_res, old_atom.name, old_atom.pos, 
                    element=old_atom.element, occupancy=old_atom.occupancy,
                    b_factor=old_atom.b_factor, is_hetatm=old_atom.is_hetatom)
            # Add bonds
            for old_atom in residue.atoms:
                for old_bond in old_atom.bonds:
                    new_first = new_res.FindAtom(old_bond.first.name)
                    new_second = new_res.FindAtom(old_bond.second.name)
                    new_editor.Connect(new_first, new_second)
            return new_res

        def _process_ligand_residue(res, rename_chain):
            """Copy or fetch the residue. Return the residue handle."""
            new_res = None
            if res.entity.handle == old_entity.handle:
                # Residue is part of the old_entity handle.
                # However, it may not be in the copied one, for instance it may have been a view
                # We try to grab it first, otherwise we copy it
                new_res = new_entity.FindResidue(res.chain.name, res.number)
            if new_res and new_res.valid:
                LogVerbose("Ligand residue %s already in entity" % res.handle.qualified_name)
            else:
                # Residue is not part of the entity, need to copy it first
                new_res = _copy_residue(res, rename_chain)
                LogVerbose("Copied ligand residue %s" % res.handle.qualified_name)
            new_res.SetIsLigand(True)
            return new_res

        for ligand in ligands:
            if isinstance(ligand, mol.EntityHandle) or isinstance(ligand, mol.EntityView):
                for residue in ligand.residues:
                    new_residue = _process_ligand_residue(residue, rename_chain)
                    extracted_ligands.append(new_residue)
            elif isinstance(ligand, mol.ResidueHandle) or isinstance(ligand, mol.ResidueView):
                new_residue = _process_ligand_residue(ligand, rename_chain)
                extracted_ligands.append(new_residue)
            else:
                raise RuntimeError("Ligands should be given as Entity or Residue")

        if new_editor is not None:
            new_editor.UpdateICS()
        return extracted_ligands

    def _get_binding_sites(self, ligand):
        """Find representations of the binding site of *ligand* in the model.

        Only consider protein and nucleic acid chains that pass the criteria
        for the :class:`ost.mol.alg.chain_mapping`. This means ignoring other
        ligands, waters, short polymers as well as any incorrectly connected
        chain that may be in proximity.

        :param ligand: Defines the binding site to identify.
        :type ligand: :class:`~ost.mol.ResidueHandle`
        """
        if ligand.hash_code not in self._binding_sites:

            # create view of reference binding site
            ref_residues_hashes = set()  # helper to keep track of added residues
            ignored_residue_hashes = {ligand.hash_code}
            for ligand_at in ligand.atoms:
                close_atoms = self.target.FindWithin(ligand_at.GetPos(), self.radius)
                for close_at in close_atoms:
                    # Skip any residue not in the chain mapping target
                    ref_res = close_at.GetResidue()
                    h = ref_res.handle.GetHashCode()
                    if h not in ref_residues_hashes and \
                            h not in ignored_residue_hashes:
                        if self.chain_mapper.target.ViewForHandle(ref_res).IsValid():
                            h = ref_res.handle.GetHashCode()
                            ref_residues_hashes.add(h)
                        elif ref_res.is_ligand:
                            LogWarning("Ignoring ligand %s in binding site of %s" % (
                                ref_res.qualified_name, ligand.qualified_name))
                            ignored_residue_hashes.add(h)
                        elif ref_res.chem_type == mol.ChemType.WATERS:
                            pass # That's ok, no need to warn
                        else:
                            LogWarning("Ignoring residue %s in binding site of %s" % (
                                ref_res.qualified_name, ligand.qualified_name))
                            ignored_residue_hashes.add(h)

            if ref_residues_hashes:
                # reason for doing that separately is to guarantee same ordering of
                # residues as in underlying entity. (Reorder by ResNum seems only
                # available on ChainHandles)
                ref_bs = self.target.CreateEmptyView()
                for ch in self.target.chains:
                    for r in ch.residues:
                        if r.handle.GetHashCode() in ref_residues_hashes:
                            ref_bs.AddResidue(r, mol.ViewAddFlag.INCLUDE_ALL)
                if len(ref_bs.residues) == 0:
                    raise RuntimeError("Failed to add proximity residues to "
                                       "the reference binding site entity")

                # Find the representations
                if self.global_chain_mapping:
                    self._binding_sites[ligand.hash_code] = self.chain_mapper.GetRepr(
                        ref_bs, self.model, inclusion_radius=self.lddt_lp_radius,
                        global_mapping = self._model_mapping)
                else:
                    self._binding_sites[ligand.hash_code] = self.chain_mapper.GetRepr(
                        ref_bs, self.model, inclusion_radius=self.lddt_lp_radius,
                        topn=self.binding_sites_topn)

                # Flag empty representation
                if not self._binding_sites[ligand.hash_code]:
                    self._unassigned_target_ligands_reason[ligand] = (
                        "model_representation",
                        "No representation of the reference binding site was "
                        "found in the model")

            else:  # if ref_residues_hashes
                # Flag missing binding site
                self._unassigned_target_ligands_reason[ligand] = ("binding_site",
                    "No residue in proximity of the target ligand")
                self._binding_sites[ligand.hash_code] = []

        return self._binding_sites[ligand.hash_code]

    @staticmethod
    def _build_binding_site_entity(ligand, residues, extra_residues=[]):
        """ Build an entity with all the binding site residues in chain A
        and the ligand in chain _. Residues are renumbered consecutively from
        1. The ligand is assigned residue number 1 and residue name LIG.
        Residues in extra_residues not in `residues` in the model are added
        at the end of chain A.

        :param ligand: the Residue Handle of the ligand
        :type ligand: :class:`~ost.mol.ResidueHandle`
        :param residues: a list of binding site residues
        :type residues: :class:`list` of :class:`~ost.mol.ResidueHandle`
        :param extra_residues: an optional list with addition binding site
                               residues. Residues in this list which are not
                               in `residues` will be added at the end of chain
                               A. This allows for instance adding unmapped
                               residues missing from the model into the
                               reference binding site.
        :type extra_residues: :class:`list` of :class:`~ost.mol.ResidueHandle`
        :rtype: :class:`~ost.mol.EntityHandle`
        """
        bs_ent = mol.CreateEntity()
        ed = bs_ent.EditXCS()
        bs_chain = ed.InsertChain("A")
        seen_res_qn = []
        for resnum, old_res in enumerate(residues, 1):
            seen_res_qn.append(old_res.qualified_name)
            new_res = ed.AppendResidue(bs_chain, old_res.handle,
                                       deep=True)
            ed.SetResidueNumber(new_res, mol.ResNum(resnum))

        # Add extra residues at the end.
        for extra_res in extra_residues:
            if extra_res.qualified_name not in seen_res_qn:
                resnum += 1
                seen_res_qn.append(extra_res.qualified_name)
                new_res = ed.AppendResidue(bs_chain,
                                           extra_res.handle,
                                           deep=True)
                ed.SetResidueNumber(new_res, mol.ResNum(resnum))
        # Add the ligand in chain _
        ligand_chain = ed.InsertChain("_")
        ligand_res = ed.AppendResidue(ligand_chain, ligand,
                                      deep=True)
        ed.RenameResidue(ligand_res, "LIG")
        ed.SetResidueNumber(ligand_res, mol.ResNum(1))
        ed.UpdateICS()

        return bs_ent

    def _compute_scores(self):
        """
        Compute the RMSD and lDDT-PLI scores for every possible target-model
        ligand pair and store the result in internal matrices.
        """
        # Create the result matrices
        rmsd_full_matrix = np.empty(
            (len(self.target_ligands), len(self.model_ligands)), dtype=dict)
        lddt_pli_full_matrix = np.empty(
            (len(self.target_ligands), len(self.model_ligands)), dtype=dict)
        self._assignment_isomorphisms = np.full(
            (len(self.target_ligands), len(self.model_ligands)), fill_value=np.nan)
        self._assignment_match_coverage = np.zeros(
            (len(self.target_ligands), len(self.model_ligands)))

        for target_i, target_ligand in enumerate(self.target_ligands):
            LogVerbose("Analyzing target ligand %s" % target_ligand)

            for binding_site in self._get_binding_sites(target_ligand):
                LogVerbose("Found binding site with chain mapping %s" % (binding_site.GetFlatChainMapping()))

                ref_bs_ent = self._build_binding_site_entity(
                    target_ligand, binding_site.ref_residues,
                    binding_site.substructure.residues)
                ref_bs_ent_ligand = ref_bs_ent.FindResidue("_", 1)  # by definition

                custom_compounds = {
                    ref_bs_ent_ligand.name:
                        mol.alg.lddt.CustomCompound.FromResidue(
                            ref_bs_ent_ligand)}
                lddt_scorer = mol.alg.lddt.lDDTScorer(
                    ref_bs_ent,
                    custom_compounds=custom_compounds,
                    inclusion_radius=self.lddt_pli_radius)

                for model_i, model_ligand in enumerate(self.model_ligands):
                    try:
                        symmetries = _ComputeSymmetries(
                            model_ligand, target_ligand,
                            substructure_match=self.substructure_match,
                            by_atom_index=True,
                            max_symmetries=self.max_symmetries)
                        LogVerbose("Ligands %s and %s symmetry match" % (
                            str(model_ligand), str(target_ligand)))
                    except NoSymmetryError:
                        # Ligands are different - skip
                        LogVerbose("No symmetry between %s and %s" % (
                            str(model_ligand), str(target_ligand)))
                        self._assignment_isomorphisms[target_i, model_i] = 0.
                        continue
                    except TooManySymmetriesError:
                        # Ligands are too symmetrical - skip
                        LogVerbose("Too many symmetries between %s and %s" % (
                            str(model_ligand), str(target_ligand)))
                        self._assignment_isomorphisms[target_i, model_i] = -1.
                        continue
                    except DisconnectedGraphError:
                        # Disconnected graph is handled elsewhere
                        continue
                    substructure_match = len(symmetries[0][0]) != len(
                        model_ligand.atoms)
                    coverage = len(symmetries[0][0]) / len(model_ligand.atoms)
                    self._assignment_match_coverage[target_i, model_i] = coverage
                    self._assignment_isomorphisms[target_i, model_i] = 1.

                    rmsd = _SCRMSD_symmetries(symmetries, model_ligand, 
                        target_ligand, transformation=binding_site.transform)
                    LogDebug("RMSD: %.4f" % rmsd)

                    # Save results?
                    if not rmsd_full_matrix[target_i, model_i] or \
                            rmsd_full_matrix[target_i, model_i]["rmsd"] > rmsd:
                        rmsd_full_matrix[target_i, model_i] = {
                            "rmsd": rmsd,
                            "lddt_lp": binding_site.lDDT,
                            "bs_ref_res": binding_site.substructure.residues,
                            "bs_ref_res_mapped": binding_site.ref_residues,
                            "bs_mdl_res_mapped": binding_site.mdl_residues,
                            "bb_rmsd": binding_site.bb_rmsd,
                            "target_ligand": target_ligand,
                            "model_ligand": model_ligand,
                            "chain_mapping": binding_site.GetFlatChainMapping(),
                            "transform": binding_site.transform,
                            "substructure_match": substructure_match,
                            "coverage": coverage,
                            "inconsistent_residues": binding_site.inconsistent_residues,
                        }
                        if self.unassigned:
                            rmsd_full_matrix[target_i, model_i][
                                "unassigned"] = False
                        LogDebug("Saved RMSD")

                    mdl_bs_ent = self._build_binding_site_entity(
                        model_ligand, binding_site.mdl_residues, [])
                    mdl_bs_ent_ligand = mdl_bs_ent.FindResidue("_", 1)  # by definition

                    # Now for each symmetry, loop and rename atoms according
                    # to ref.
                    mdl_editor = mdl_bs_ent.EditXCS()
                    for i, (trg_sym, mdl_sym) in enumerate(symmetries):
                        # Prepare Entities for RMSD
                        for mdl_anum, trg_anum in zip(mdl_sym, trg_sym):
                            # Rename model atoms according to symmetry
                            trg_atom = ref_bs_ent_ligand.atoms[trg_anum]
                            mdl_atom = mdl_bs_ent_ligand.atoms[mdl_anum]
                            mdl_editor.RenameAtom(mdl_atom, trg_atom.name)
                        mdl_editor.UpdateICS()

                        global_lddt, local_lddt, lddt_tot, lddt_cons, n_res, \
                            n_cont, n_cons = lddt_scorer.lDDT(
                                mdl_bs_ent, chain_mapping={"A": "A", "_": "_"},
                                no_intrachain=True,
                                return_dist_test=True,
                                check_resnames=self.check_resnames)
                        LogDebug("lDDT-PLI for symmetry %d: %.4f" % (i, global_lddt))

                        # Save results?
                        if not lddt_pli_full_matrix[target_i, model_i]:
                            # First iteration
                            save_lddt = True
                        else:
                            last_best_lddt = lddt_pli_full_matrix[
                                target_i, model_i]["lddt_pli"]
                            last_best_rmsd = lddt_pli_full_matrix[
                                target_i, model_i]["rmsd"]
                            if global_lddt > last_best_lddt:
                                # Better lDDT-PLI
                                save_lddt = True
                            elif global_lddt == last_best_lddt and \
                                    rmsd < last_best_rmsd:
                                # Same lDDT-PLI, better RMSD
                                save_lddt = True
                            else:
                                save_lddt = False
                        if save_lddt:
                            lddt_pli_full_matrix[target_i, model_i] = {
                                "lddt_pli": global_lddt,
                                "rmsd": rmsd,
                                "lddt_lp": binding_site.lDDT,
                                "lddt_pli_n_contacts": lddt_tot,
                                "bs_ref_res": binding_site.substructure.residues,
                                "bs_ref_res_mapped": binding_site.ref_residues,
                                "bs_mdl_res_mapped": binding_site.mdl_residues,
                                "bb_rmsd": binding_site.bb_rmsd,
                                "target_ligand": target_ligand,
                                "model_ligand": model_ligand,
                                "chain_mapping": binding_site.GetFlatChainMapping(),
                                "transform": binding_site.transform,
                                "substructure_match": substructure_match,
                                "coverage": coverage,
                                "inconsistent_residues": binding_site.inconsistent_residues,
                            }
                            if self.unassigned:
                                lddt_pli_full_matrix[target_i, model_i][
                                    "unassigned"] = False
                            LogDebug("Saved lDDT-PLI")

        self._rmsd_full_matrix = rmsd_full_matrix
        self._lddt_pli_full_matrix = lddt_pli_full_matrix

    @staticmethod
    def _find_ligand_assignment(mat1, mat2=None, coverage=None, coverage_delta=None):
        """ Find the ligand assignment based on mat1. If mat2 is provided, it
        will be used to break ties in mat1. If mat2 is not provided, ties will
        be resolved by taking the first match arbitrarily.

        Both mat1 and mat2 should "look" like RMSD - ie be between inf (bad)
        and 0 (good).
        """
        # We will modify mat1 and mat2, so make copies of it first
        mat1 = np.copy(mat1)
        if mat2 is None:
            mat2 = np.copy(mat1)
            mat2[~np.isnan(mat2)] = np.inf
        else:
            mat2 = np.copy(mat2)
        if coverage is None:
            coverage = np.copy(mat1)
            coverage[:] = 1  # Assume full coverage by default
        else:
            coverage = np.copy(coverage)

        assignments = []
        if 0 in mat1.shape:
            # No model or target ligand
            LogDebug("No model or target ligand, returning no assignment.")
            return assignments

        def _get_best_match(mat1_val, coverage_val):
            """ Extract the row/column indices of the prediction matching the
                given values."""
            mat1_match_idx = np.argwhere((mat1 == mat1_val) & (coverage >= coverage_val))
            # Multiple "best" - use mat2 to disambiguate
            if len(mat1_match_idx) > 1:
                # Get the values of mat2 at these positions
                best_mat2_match = [mat2[tuple(x)] for x in mat1_match_idx]
                # Find the index of the best mat2
                # Note: argmin returns the first value which is min.
                best_mat2_idx = np.array(best_mat2_match).argmin()
                # Now get the original indices
                return mat1_match_idx[best_mat2_idx]
            else:
                return mat1_match_idx[0]

        # First only consider top coverage matches.
        min_coverage = np.max(coverage)
        while min_coverage > 0:
            LogVerbose("Looking for matches with coverage >= %s" % min_coverage)
            min_mat1 = LigandScorer._nanmin_nowarn(mat1, coverage < min_coverage)
            while not np.isnan(min_mat1):
                max_i_trg, max_i_mdl = _get_best_match(min_mat1, min_coverage)

                # Would we have a match for this model ligand with higher score
                # but lower coverage?
                alternative_matches = (mat1[:, max_i_mdl] < min_mat1) & (
                        coverage[:, max_i_mdl] > (min_coverage - coverage_delta))
                if np.any(alternative_matches):
                    # Get the scores of these matches
                    LogVerbose("Found match with lower coverage but better score")
                    min_mat1 = np.nanmin(mat1[alternative_matches])
                    max_i_trg, max_i_mdl = _get_best_match(min_mat1, min_coverage - coverage_delta)

                # Disable row and column
                mat1[max_i_trg, :] = np.nan
                mat1[:, max_i_mdl] = np.nan
                mat2[max_i_trg, :] = np.nan
                mat2[:, max_i_mdl] = np.nan
                coverage[max_i_trg, :] = -np.inf
                coverage[:, max_i_mdl] = -np.inf

                # Save
                assignments.append((max_i_trg, max_i_mdl))

                # Recompute min
                min_mat1 = LigandScorer._nanmin_nowarn(mat1, coverage < min_coverage)
            # Recompute min_coverage
            min_coverage = np.max(coverage)
        return assignments

    @staticmethod
    def _nanmin_nowarn(array, mask):
        """Compute np.nanmin but ignore the RuntimeWarning."""
        masked_array = np_ma.masked_array(array, mask=mask)
        with warnings.catch_warnings():  # RuntimeWarning: All-NaN slice encountered
            warnings.simplefilter("ignore")
            min = np.nanmin(masked_array, )
            if np_ma.is_masked(min):
                return np.nan  # Everything was masked
            else:
                return min

    @staticmethod
    def _reverse_lddt(lddt):
        """Reverse lDDT means turning it from a number between 0 and 1 to a
        number between infinity and 0 (0 being better).

        In practice, this is 1/lDDT. If lDDT is 0, the result is infinity.
        """
        with warnings.catch_warnings():  # RuntimeWarning: divide by zero
            warnings.simplefilter("ignore")
            return np.float64(1) / lddt

    def _assign_ligands_rmsd(self):
        """Assign (map) ligands between model and target.

        Sets self._rmsd and self._rmsd_details.
        """
        mat2 = self._reverse_lddt(self.lddt_pli_matrix)

        mat_tuple = self._assign_matrices(self.rmsd_matrix,
                                          mat2,
                                          self._rmsd_full_matrix,
                                          "rmsd")
        self._rmsd = mat_tuple[0]
        self._rmsd_details = mat_tuple[1]
        # Ignore unassigned ligands - they are dealt with in lddt_pli.
        # So the following lines should stay commented out:
        # self._unassigned_target_ligands = mat_tuple[2]
        # self._unassigned_model_ligands = mat_tuple[3]

    def _assign_matrices(self, mat1, mat2, data, main_key):
        """
        Perform the ligand assignment, ie find the mapping between model and
        target ligands.

        The algorithm starts by assigning the "best" mapping, and then discards
        the target and model ligands (row, column) so that every model ligand
        can be assigned to a single target ligand, and every target ligand
        is only assigned to a single model ligand. Repeat until there is
        nothing left to assign.

        In case of a tie in values in `mat1`, it uses `mat2` to break the tie.

        This algorithm doesn't guarantee a globally optimal assignment.

        Both `mat1` and `mat2` should contain values between 0 and infinity,
        with lower values representing better scores. Use the
        :meth:`_reverse_lddt` method to convert lDDT values to such a score.

        :param mat1: the main ligand assignment criteria (RMSD or lDDT-PLI)
        :param mat2: the secondary ligand assignment criteria (lDDT-PLI or RMSD)
        :param data: the data (either self._rmsd_full_matrix or self._lddt_pli_matrix)
        :param main_key: the key of data (dictionnaries within `data`) to
               assign into out_main.
        :return: a tuple with 2 dictionaries of matrices containing the main
                 data, and details, respectively.
        """
        assignments = self._find_ligand_assignment(mat1, mat2,
                                                   self._assignment_match_coverage,
                                                   self.coverage_delta)
        out_main = {}
        out_details = {}
        assigned_trg = [False] * len(self.target_ligands)
        assigned_mdl = [False] * len(self.model_ligands)
        for assignment in assignments:
            trg_idx, mdl_idx = assignment
            assigned_mdl[mdl_idx] = True
            assigned_trg[trg_idx] = True
            mdl_lig = self.model_ligands[mdl_idx]
            mdl_cname = mdl_lig.chain.name
            mdl_resnum = mdl_lig.number
            if mdl_cname not in out_main:
                out_main[mdl_cname] = {}
                out_details[mdl_cname] = {}
            out_main[mdl_cname][mdl_resnum] = data[
                trg_idx, mdl_idx][main_key]
            out_details[mdl_cname][mdl_resnum] = data[
                trg_idx, mdl_idx]


        unassigned_trg, unassigned_mdl = self._assign_unassigned(
            assigned_trg, assigned_mdl, [out_main], [out_details], [main_key])
        return out_main, out_details, unassigned_trg, unassigned_mdl

    def _assign_unassigned(self, assigned_trg, assigned_mdl,
                           out_main, out_details, main_key):
        unassigned_trg = {}
        unassigned_mdl = {}

        unassigned_trg_idx = [i for i, x in enumerate(assigned_trg) if not x]
        unassigned_mdl_idx = [i for i, x in enumerate(assigned_mdl) if not x]

        for mdl_idx in unassigned_mdl_idx:
            mdl_lig = self.model_ligands[mdl_idx]
            reason = self._find_unassigned_model_ligand_reason(mdl_lig, check=False)
            mdl_cname = mdl_lig.chain.name
            mdl_resnum = mdl_lig.number
            if mdl_cname not in unassigned_mdl:
                unassigned_mdl[mdl_cname] = {}
            unassigned_mdl[mdl_cname][mdl_resnum] = reason
            if self.unassigned:
                for i, _ in enumerate(out_main):
                    if mdl_cname not in out_main[i]:
                        out_main[i][mdl_cname] = {}
                        out_details[i][mdl_cname] = {}
                    out_main[i][mdl_cname][mdl_resnum] = None
                    out_details[i][mdl_cname][mdl_resnum] = {
                        "unassigned": True,
                        "reason_short": reason[0],
                        "reason_long": reason[1],
                        main_key[i]: None,
                    }
                    LogInfo("Model ligand %s is unassigned: %s" % (
                        mdl_lig.qualified_name, reason[1]))

        for trg_idx in unassigned_trg_idx:
            trg_lig = self.target_ligands[trg_idx]
            reason = self._find_unassigned_target_ligand_reason(trg_lig, check=False)
            trg_cname = trg_lig.chain.name
            trg_resnum = trg_lig.number
            if trg_cname not in unassigned_trg:
                unassigned_trg[trg_cname] = {}
            unassigned_trg[trg_cname][trg_resnum] = reason
            LogInfo("Target ligand %s is unassigned: %s" % (
                trg_lig.qualified_name, reason[1]))

        return unassigned_trg, unassigned_mdl


    def _assign_matrix(self, mat, data1, main_key1, data2, main_key2):
        """
        Perform the ligand assignment, ie find the mapping between model and
        target ligands, based on a single matrix

        The algorithm starts by assigning the "best" mapping, and then discards
        the target and model ligands (row, column) so that every model ligand
        can be assigned to a single target ligand, and every target ligand
        is only assigned to a single model ligand. Repeat until there is
        nothing left to assign.

        This algorithm doesn't guarantee a globally optimal assignment.

        `mat` should contain values between 0 and infinity,
        with lower values representing better scores. Use the
        :meth:`_reverse_lddt` method to convert lDDT values to such a score.

        :param mat: the ligand assignment criteria (RMSD or lDDT-PLI)
        :param data1: the first data (either self._rmsd_full_matrix or self._lddt_pli_matrix)
        :param main_key1: the first key of data (dictionnaries within `data`) to
               assign into out_main.
        :param data2: the second data (either self._rmsd_full_matrix or self._lddt_pli_matrix)
        :param main_key2: the second key of data (dictionnaries within `data`) to
               assign into out_main.
        :return: a tuple with 4 dictionaries of matrices containing the main
                 data1, details1, main data2 and details2, respectively.
        """
        assignments = self._find_ligand_assignment(mat,
                                                   coverage=self._assignment_match_coverage,
                                                   coverage_delta=self.coverage_delta)
        out_main1 = {}
        out_details1 = {}
        out_main2 = {}
        out_details2 = {}
        assigned_trg = [False] * len(self.target_ligands)
        assigned_mdl = [False] * len(self.model_ligands)
        for assignment in assignments:
            trg_idx, mdl_idx = assignment
            assigned_mdl[mdl_idx] = True
            assigned_trg[trg_idx] = True
            mdl_lig = self.model_ligands[mdl_idx]
            mdl_cname = mdl_lig.chain.name
            mdl_resnum = mdl_lig.number
            # Data 1
            if mdl_cname not in out_main1:
                out_main1[mdl_cname] = {}
                out_details1[mdl_cname] = {}
            out_main1[mdl_cname][mdl_resnum] = data1[
                trg_idx, mdl_idx][main_key1]
            out_details1[mdl_cname][mdl_resnum] = data1[
                trg_idx, mdl_idx]
            # Data2
            if mdl_cname not in out_main2:
                out_main2[mdl_cname] = {}
                out_details2[mdl_cname] = {}
            out_main2[mdl_cname][mdl_resnum] = data2[
                trg_idx, mdl_idx][main_key2]
            out_details2[mdl_cname][mdl_resnum] = data2[
                trg_idx, mdl_idx]

        unassigned_trg, unassigned_mdl = self._assign_unassigned(
            assigned_trg, assigned_mdl,
            [out_main1, out_main2], [out_details1, out_details2],
            [main_key1, main_key2])

        return out_main1, out_details1, out_main2, out_details2, \
            unassigned_trg, unassigned_mdl

    def _assign_ligands_lddt_pli(self):
        """ Assign ligands based on lDDT-PLI.

        Sets self._lddt_pli and self._lddt_pli_details.
        """
        mat1 = self._reverse_lddt(self.lddt_pli_matrix)

        mat_tuple = self._assign_matrices(mat1,
                                          self.rmsd_matrix,
                                          self._lddt_pli_full_matrix,
                                          "lddt_pli")
        self._lddt_pli = mat_tuple[0]
        self._lddt_pli_details = mat_tuple[1]
        self._unassigned_target_ligands = mat_tuple[2]
        self._unassigned_model_ligands = mat_tuple[3]

    def _assign_ligands_rmsd_only(self):
        """Assign (map) ligands between model and target based on RMSD only.

        Sets self._rmsd, self._rmsd_details, self._lddt_pli and
        self._lddt_pli_details.
        """
        mat_tuple = self._assign_matrix(self.rmsd_matrix,
                                        self._rmsd_full_matrix,
                                        "rmsd",
                                        self._lddt_pli_full_matrix,
                                        "lddt_pli")
        self._rmsd = mat_tuple[0]
        self._rmsd_details = mat_tuple[1]
        self._lddt_pli = mat_tuple[2]
        self._lddt_pli_details = mat_tuple[3]
        self._unassigned_target_ligands = mat_tuple[4]
        self._unassigned_model_ligands = mat_tuple[5]

    @property
    def rmsd_matrix(self):
        """ Get the matrix of RMSD values.

        Target ligands are in rows, model ligands in columns.

        NaN values indicate that no RMSD could be computed (i.e. different
        ligands).

        :rtype: :class:`~numpy.ndarray`
        """
        if self._rmsd_full_matrix is None:
            self._compute_scores()
        if self._rmsd_matrix is None:
            # convert
            shape = self._rmsd_full_matrix.shape
            self._rmsd_matrix = np.full(shape, np.nan)
            for i, j in np.ndindex(shape):
                if self._rmsd_full_matrix[i, j] is not None:
                    self._rmsd_matrix[i, j] = self._rmsd_full_matrix[
                        i, j]["rmsd"]
        return self._rmsd_matrix

    @property
    def lddt_pli_matrix(self):
        """ Get the matrix of lDDT-PLI values.

        Target ligands are in rows, model ligands in columns.

        NaN values indicate that no lDDT-PLI could be computed (i.e. different
        ligands).

        :rtype: :class:`~numpy.ndarray`
        """
        if self._lddt_pli_full_matrix is None:
            self._compute_scores()
        if self._lddt_pli_matrix is None:
            # convert
            shape = self._lddt_pli_full_matrix.shape
            self._lddt_pli_matrix = np.full(shape, np.nan)
            for i, j in np.ndindex(shape):
                if self._lddt_pli_full_matrix[i, j] is not None:
                    self._lddt_pli_matrix[i, j] = self._lddt_pli_full_matrix[
                        i, j]["lddt_pli"]
        return self._lddt_pli_matrix

    @property
    def coverage_matrix(self):
        """ Get the matrix of model ligand atom coverage in the target.

        Target ligands are in rows, model ligands in columns.

        A value of 0 indicates that there was no isomorphism between the model
        and target ligands. If `substructure_match=False`, only full match
        isomorphisms are considered, and therefore only values of 1.0 and 0.0
        are reported.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._assignment_match_coverage is None:
            self._compute_scores()
        return self._assignment_match_coverage

    @property
    def rmsd(self):
        """Get a dictionary of RMSD score values, keyed by model ligand
        (chain name, :class:`~ost.mol.ResNum`).

        If the scoring object was instantiated with `unassigned=True`, some
        scores may be `None`.

        :rtype: :class:`dict`
        """
        if self._rmsd is None:
            if self.rmsd_assignment:
                self._assign_ligands_rmsd_only()
            else:
                self._assign_ligands_rmsd()
        return self._rmsd

    @property
    def rmsd_details(self):
        """Get a dictionary of RMSD score details (dictionaries), keyed by
        model ligand (chain name, :class:`~ost.mol.ResNum`).

        The value is a dictionary. For ligands that were assigned (mapped) to
        the target, the dictionary contain the following information:

        * `rmsd`: the RMSD score value.
        * `lddt_lp`: the lDDT score of the ligand pocket (lDDT-LP).
        * `bs_ref_res`: a list of residues (:class:`~ost.mol.ResidueHandle`)
          that define the binding site in the reference.
        * `bs_ref_res_mapped`: a list of residues
          (:class:`~ost.mol.ResidueHandle`) in the reference binding site
          that could be mapped to the model.
        * `bs_mdl_res_mapped`: a list of residues
          (:class:`~ost.mol.ResidueHandle`) in the model that were mapped to
          the reference binding site. The residues are in the same order as
          `bs_ref_res_mapped`.
        * `bb_rmsd`: the RMSD of the binding site backbone after superposition
        * `target_ligand`: residue handle of the target ligand.
        * `model_ligand`: residue handle of the model ligand.
        * `chain_mapping`: local chain mapping as a dictionary, with target
          chain name as key and model chain name as value.
        * `transform`: transformation to superpose the model onto the target.
        * `substructure_match`: whether the score is the result of a partial
          (substructure) match. A value of `True` indicates that the target
          ligand covers only part of the model, while `False` indicates a
          perfect match.
        * `coverage`: the fraction of model atoms covered by the assigned
          target ligand, in the interval (0, 1]. If `substructure_match`
          is `False`, this will always be 1.
        * `inconsistent_residues`: a list of tuples of mapped residues views
          (:class:`~ost.mol.ResidueView`) with residue names that differ
          between the reference and the model, respectively.
          The list is empty if all residue names match, which is guaranteed
          if `check_resnames=True`.
          Note: more binding site mappings may be explored during scoring,
          but only inconsistencies in the selected mapping are reported.
        * `unassigned`: only if the scorer was instantiated with
          `unassigned=True`: `False`

        If the scoring object was instantiated with `unassigned=True`, in
        addition the unassigned ligands will be reported with a score of `None`
        and the following information:

        * `unassigned`: `True`,
        * `reason_short`: a short token of the reason, see
          :attr:`unassigned_model_ligands` for details and meaning.
        * `reason_long`: a human-readable text of the reason, see
          :attr:`unassigned_model_ligands` for details and meaning.
        * `rmsd`: `None`

        :rtype: :class:`dict`
        """
        if self._rmsd_details is None:
            if self.rmsd_assignment:
                self._assign_ligands_rmsd_only()
            else:
                self._assign_ligands_rmsd()
        return self._rmsd_details

    @property
    def lddt_pli(self):
        """Get a dictionary of lDDT-PLI score values, keyed by model ligand
        (chain name, :class:`~ost.mol.ResNum`).

        If the scoring object was instantiated with `unassigned=True`, some
        scores may be `None`.

        :rtype: :class:`dict`
        """
        if self._lddt_pli is None:
            if self.rmsd_assignment:
                self._assign_ligands_rmsd_only()
            else:
                self._assign_ligands_lddt_pli()
        return self._lddt_pli

    @property
    def lddt_pli_details(self):
        """Get a dictionary of lDDT-PLI score details (dictionaries), keyed by
        model ligand (chain name, :class:`~ost.mol.ResNum`).

        Each sub-dictionary contains the following information:

        * `lddt_pli`: the lDDT-PLI score value.
        * `rmsd`: the RMSD score value corresponding to the lDDT-PLI
          chain mapping and assignment. This may differ from the RMSD-based
          assignment. Note that a different isomorphism than `lddt_pli` may
          be used.
        * `lddt_lp`: the lDDT score of the ligand pocket (lDDT-LP).
        * `lddt_pli_n_contacts`: number of total contacts used in lDDT-PLI,
          summed over all thresholds. Can be divided by 8 to obtain the number
          of atomic contacts.
        * `bs_ref_res`: a list of residues (:class:`~ost.mol.ResidueHandle`)
          that define the binding site in the reference.
        * `bs_ref_res_mapped`: a list of residues
          (:class:`~ost.mol.ResidueHandle`) in the reference binding site
          that could be mapped to the model.
        * `bs_mdl_res_mapped`: a list of residues
          (:class:`~ost.mol.ResidueHandle`) in the model that were mapped to
          the reference binding site. The residues are in the same order as
          `bs_ref_res_mapped`.
        * `bb_rmsd`: the RMSD of the binding site backbone after superposition.
          Note: not used for lDDT-PLI computation.
        * `target_ligand`: residue handle of the target ligand.
        * `model_ligand`: residue handle of the model ligand.
        * `chain_mapping`: local chain mapping as a dictionary, with target
          chain name as key and model chain name as value.
        * `transform`: transformation to superpose the model onto the target
          (for RMSD only).
        * `substructure_match`: whether the score is the result of a partial
          (substructure) match. A value of `True` indicates that the target
          ligand covers only part of the model, while `False` indicates a
          perfect match.
        * `inconsistent_residues`: a list of tuples of mapped residues views
          (:class:`~ost.mol.ResidueView`) with residue names that differ
          between the reference and the model, respectively.
          The list is empty if all residue names match, which is guaranteed
          if `check_resnames=True`.
          Note: more binding site mappings may be explored during scoring,
          but only inconsistencies in the selected mapping are reported.
        * `unassigned`: only if the scorer was instantiated with
          `unassigned=True`: `False`

        If the scoring object was instantiated with `unassigned=True`, in
        addition the unmapped ligands will be reported with a score of `None`
        and the following information:

        * `unassigned`: `True`,
        * `reason_short`: a short token of the reason, see
          :attr:`unassigned_model_ligands` for details and meaning.
        * `reason_long`: a human-readable text of the reason, see
          :attr:`unassigned_model_ligands` for details and meaning.
        * `lddt_pli`: `None`

        :rtype: :class:`dict`
        """
        if self._lddt_pli_details is None:
            if self.rmsd_assignment:
                self._assign_ligands_rmsd_only()
            else:
                self._assign_ligands_lddt_pli()
        return self._lddt_pli_details

    @property
    def unassigned_target_ligands(self):
        """Get a dictionary of target ligands not assigned to any model ligand,
        keyed by target ligand (chain name, :class:`~ost.mol.ResNum`).

        The assignment for the lDDT-PLI score is used (and is controlled
        by the `rmsd_assignment` argument).

        Each item contains a string from a controlled dictionary
        about the reason for the absence of assignment.
        A human-readable description can be obtained from the
        :attr:`unassigned_target_ligand_descriptions` property.

        Currently, the following reasons are reported:

        * `no_ligand`: there was no ligand in the model.
        * `disconnected`: the ligand graph was disconnected.
        * `binding_site`: no residues were in proximity of the ligand.
        * `model_representation`: no representation of the reference binding
          site was found in the model. (I.e. the binding site was not modeled.
          Remember: the binding site is defined in the target structure,
          the position of the model ligand itself is ignored at this point.)
        * `identity`: the ligand was not found in the model (by graph
          isomorphism). Check your ligand connectivity, and enable the
          `substructure_match` option if the target ligand is incomplete.
        * `stoichiometry`: there was a possible assignment in the model, but
          the model ligand was already assigned to a different target ligand.
          This indicates different stoichiometries.
        * `symmetries`: too many symmetries were found (by graph isomorphisms).
          Increase `max_symmetries`.

        Some of these reasons can be overlapping, but a single reason will be
        reported.

        :rtype: :class:`dict`
        """
        if self._unassigned_target_ligand_short is None:
            if self.rmsd_assignment:
                self._assign_ligands_rmsd_only()
            else:
                self._assign_ligands_lddt_pli()
            self._unassigned_target_ligand_short = {}
            self._unassigned_target_ligand_descriptions = {}
            for cname, res in self._unassigned_target_ligands.items():
                self._unassigned_target_ligand_short[cname] = {}
                for resnum, val in res.items():
                    self._unassigned_target_ligand_short[cname][resnum] = val[0]
                    self._unassigned_target_ligand_descriptions[val[0]] = val[1]
        return self._unassigned_target_ligand_short

    @property
    def unassigned_target_ligand_descriptions(self):
        """Get a human-readable description of why target ligands were
        unassigned, as a dictionary keyed by the controlled dictionary
        from :attr:`unassigned_target_ligands`.
        """
        if self._unassigned_target_ligand_descriptions is None:
            _ = self.unassigned_target_ligands  # assigned there
        return self._unassigned_target_ligand_descriptions

    @property
    def unassigned_model_ligands(self):
        """Get a dictionary of model ligands not assigned to any target ligand,
        keyed by model ligand (chain name, :class:`~ost.mol.ResNum`).

        The assignment for the lDDT-PLI score is used (and is controlled
        by the `rmsd_assignment` argument).

        Each item contains a string from a controlled dictionary
        about the reason for the absence of assignment.
        A human-readable description can be obtained from the
        :attr:`unassigned_model_ligand_descriptions` property.
        Currently, the following reasons are reported:

        * `no_ligand`: there was no ligand in the target.
        * `disconnected`: the ligand graph is disconnected.
        * `binding_site`: a potential assignment was found in the target, but
          there were no polymer residues in proximity of the ligand in the
          target.
        * `model_representation`: a potential assignment was found in the target,
          but no representation of the binding site was found in the model.
          (I.e. the binding site was not modeled. Remember: the binding site
          is defined in the target structure, the position of the model ligand
          itself is ignored at this point.)
        * `identity`: the ligand was not found in the target (by graph
          isomorphism). Check your ligand connectivity, and enable the
          `substructure_match` option if the target ligand is incomplete.
        * `stoichiometry`: there was a possible assignment in the target, but
          the model target was already assigned to a different model ligand.
          This indicates different stoichiometries.
        * `symmetries`: too many symmetries were found (by graph isomorphisms).
          Increase `max_symmetries`.

        Some of these reasons can be overlapping, but a single reason will be
        reported.

        :rtype: :class:`dict`
        """
        if self._unassigned_model_ligand_short is None:
            if self.rmsd_assignment:
                self._assign_ligands_rmsd_only()
            else:
                self._assign_ligands_lddt_pli()
            self._unassigned_model_ligand_short = {}
            self._unassigned_model_ligand_descriptions = {}
            for cname, res in self._unassigned_model_ligands.items():
                self._unassigned_model_ligand_short[cname] = {}
                for resnum, val in res.items():
                    self._unassigned_model_ligand_short[cname][resnum] = val[0]
                    self._unassigned_model_ligand_descriptions[val[0]] = val[1]
        return self._unassigned_model_ligand_short

    @property
    def unassigned_model_ligand_descriptions(self):
        """Get a human-readable description of why model ligands were
        unassigned, as a dictionary keyed by the controlled dictionary
        from :attr:`unassigned_model_ligands`.
        """
        if self._unassigned_model_ligand_descriptions is None:
            _ = self.unassigned_model_ligands  # assigned there
        return self._unassigned_model_ligand_descriptions


    def _set_custom_mapping(self, mapping):
        """ sets self.__model_mapping with a full blown MappingResult object

        :param mapping: mapping with trg chains as key and mdl ch as values
        :type mapping: :class:`dict`
        """
        chain_mapper = self.chain_mapper
        chem_mapping, chem_group_alns, mdl = \
        chain_mapper.GetChemMapping(self.model)

        # now that we have a chem mapping, lets do consistency checks
        # - check whether chain names are unique and available in structures
        # - check whether the mapped chains actually map to the same chem groups
        if len(mapping) != len(set(mapping.keys())):
            raise RuntimeError(f"Expect unique trg chain names in mapping. Got "
                               f"{mapping.keys()}")
        if len(mapping) != len(set(mapping.values())):
            raise RuntimeError(f"Expect unique mdl chain names in mapping. Got "
                               f"{mapping.values()}")

        trg_chains = set([ch.GetName() for ch in chain_mapper.target.chains])
        mdl_chains = set([ch.GetName() for ch in mdl.chains])
        for k,v in mapping.items():
            if k not in trg_chains:
                raise RuntimeError(f"Target chain \"{k}\" is not available "
                                   f"in target processed for chain mapping "
                                   f"({trg_chains})")
            if v not in mdl_chains:
                raise RuntimeError(f"Model chain \"{v}\" is not available "
                                   f"in model processed for chain mapping "
                                   f"({mdl_chains})")

        for trg_ch, mdl_ch in mapping.items():
            trg_group_idx = None
            mdl_group_idx = None
            for idx, group in enumerate(chain_mapper.chem_groups):
                if trg_ch in group:
                    trg_group_idx = idx
                    break
            for idx, group in enumerate(chem_mapping):
                if mdl_ch in group:
                    mdl_group_idx = idx
                    break
            if trg_group_idx is None or mdl_group_idx is None:
                raise RuntimeError("Could not establish a valid chem grouping "
                                   "of chain names provided in custom mapping.")
            
            if trg_group_idx != mdl_group_idx:
                raise RuntimeError(f"Chem group mismatch in custom mapping: "
                                   f"target chain \"{trg_ch}\" groups with the "
                                   f"following chemically equivalent target "
                                   f"chains: "
                                   f"{chain_mapper.chem_groups[trg_group_idx]} "
                                   f"but model chain \"{mdl_ch}\" maps to the "
                                   f"following target chains: "
                                   f"{chain_mapper.chem_groups[mdl_group_idx]}")

        pairs = set([(trg_ch, mdl_ch) for trg_ch, mdl_ch in mapping.items()])
        ref_mdl_alns =  \
        chain_mapping._GetRefMdlAlns(chain_mapper.chem_groups,
                                     chain_mapper.chem_group_alignments,
                                     chem_mapping,
                                     chem_group_alns,
                                     pairs = pairs)

        # translate mapping format
        final_mapping = list()
        for ref_chains in chain_mapper.chem_groups:
            mapped_mdl_chains = list()
            for ref_ch in ref_chains:
                if ref_ch in mapping:
                    mapped_mdl_chains.append(mapping[ref_ch])
                else:
                    mapped_mdl_chains.append(None)
            final_mapping.append(mapped_mdl_chains)

        alns = dict()
        for ref_group, mdl_group in zip(chain_mapper.chem_groups,
                                        final_mapping):
            for ref_ch, mdl_ch in zip(ref_group, mdl_group):
                if ref_ch is not None and mdl_ch is not None:
                    aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                    trg_view = chain_mapper.target.Select(f"cname={ref_ch}")
                    mdl_view = mdl.Select(f"cname={mdl_ch}")
                    aln.AttachView(0, trg_view)
                    aln.AttachView(1, mdl_view)
                    alns[(ref_ch, mdl_ch)] = aln

        self.__model_mapping = chain_mapping.MappingResult(chain_mapper.target, mdl,
                                                           chain_mapper.chem_groups,
                                                           chem_mapping,
                                                           final_mapping, alns)

    def _find_unassigned_model_ligand_reason(self, ligand, assignment="lddt_pli", check=True):
        # Is this a model ligand?
        try:
            ligand_idx = self.model_ligands.index(ligand)
        except ValueError:
            # Raise with a better error message
            raise ValueError("Ligand %s is not in self.model_ligands" % ligand)

        # Ensure we are unassigned
        if check:
            details = getattr(self, assignment + "_details")
            if ligand.chain.name in details and ligand.number in details[ligand.chain.name]:
                ligand_details = details[ligand.chain.name][ligand.number]
                if not ("unassigned" in ligand_details and ligand_details["unassigned"]):
                    raise RuntimeError("Ligand %s is mapped to %s" % (ligand, ligand_details["target_ligand"]))

        # Were there any ligands in the target?
        if len(self.target_ligands) == 0:
            return ("no_ligand", "No ligand in the target")

        # Is the ligand disconnected?
        graph = _ResidueToGraph(ligand)
        if not networkx.is_connected(graph):
            return ("disconnected", "Ligand graph is disconnected")

        # Do we have isomorphisms with the target?
        for trg_lig_idx, assigned in enumerate(self._assignment_isomorphisms[:, ligand_idx]):
            if np.isnan(assigned):
                try:
                    _ComputeSymmetries(
                        self.model_ligands[ligand_idx],
                        self.target_ligands[trg_lig_idx],
                        substructure_match=self.substructure_match,
                        by_atom_index=True,
                        return_symmetries=False)
                except (NoSymmetryError, DisconnectedGraphError):
                    assigned = 0.
                except TooManySymmetriesError:
                    assigned = -1.
                else:
                    assigned = 1.
                self._assignment_isomorphisms[trg_lig_idx,ligand_idx] = assigned
            if assigned == 1.:
                # Could have been assigned
                # So what's up with this target ligand?
                assignment_matrix = getattr(self, assignment + "_matrix")
                all_nan = np.all(np.isnan(assignment_matrix[:, ligand_idx]))
                if all_nan:
                    # The assignment matrix is all nans so we have a problem
                    # with the binding site or the representation
                    trg_ligand = self.target_ligands[trg_lig_idx]
                    return self._unassigned_target_ligands_reason[trg_ligand]
                else:
                    # Ligand was already assigned
                    return ("stoichiometry",
                            "Ligand was already assigned to an other "
                            "model ligand (different stoichiometry)")
            elif assigned == -1:
                # Symmetries / isomorphisms exceeded limit
                return ("symmetries",
                        "Too many symmetries were found.")

        # Could not be assigned to any ligand - must be different
        if self.substructure_match:
            iso = "subgraph isomorphism"
        else:
            iso = "full graph isomorphism"
        return ("identity", "Ligand was not found in the target (by %s)" % iso)

    def _find_unassigned_target_ligand_reason(self, ligand, assignment="lddt_pli", check=True):
        # Is this a target ligand?
        try:
            ligand_idx = self.target_ligands.index(ligand)
        except ValueError:
            # Raise with a better error message
            raise ValueError("Ligand %s is not in self.target_ligands" % ligand)

        # Ensure we are unassigned
        if check:
            details = getattr(self, assignment + "_details")
            for cname, chain_ligands in details.items():
                for rnum, details in chain_ligands.items():
                    if "unassigned" in details and details["unassigned"]:
                        continue
                    if details['target_ligand'] == ligand:
                        raise RuntimeError("Ligand %s is mapped to %s.%s" % (
                            ligand, cname, rnum))

        # Were there any ligands in the model?
        if len(self.model_ligands) == 0:
            return ("no_ligand", "No ligand in the model")

        # Is the ligand disconnected?
        graph = _ResidueToGraph(ligand)
        if not networkx.is_connected(graph):
            return ("disconnected", "Ligand graph is disconnected")

        # Is it because there was no valid binding site or no representation?
        if ligand in self._unassigned_target_ligands_reason:
            return self._unassigned_target_ligands_reason[ligand]
        # Or because no symmetry?
        for model_lig_idx, assigned in enumerate(
                self._assignment_isomorphisms[ligand_idx, :]):
            if np.isnan(assigned):
                try:
                    _ComputeSymmetries(
                        self.model_ligands[model_lig_idx],
                        self.target_ligands[ligand_idx],
                        substructure_match=self.substructure_match,
                        by_atom_index=True,
                        return_symmetries=False)
                except (NoSymmetryError, DisconnectedGraphError):
                    assigned = 0.
                except TooManySymmetriesError:
                    assigned = -1.
                else:
                    assigned = 1.
                self._assignment_isomorphisms[ligand_idx,model_lig_idx] = assigned
            if assigned == 1:
                # Could have been assigned but was assigned to a different ligand
                return ("stoichiometry",
                        "Ligand was already assigned to an other "
                        "target ligand (different stoichiometry)")
            elif assigned == -1:
                # Symmetries / isomorphisms exceeded limit
                return ("symmetries",
                        "Too many symmetries were found.")

        # Could not be assigned to any ligand - must be different
        if self.substructure_match:
            iso = "subgraph isomorphism"
        else:
            iso = "full graph isomorphism"
        return ("identity", "Ligand was not found in the model (by %s)" % iso)


def _ResidueToGraph(residue, by_atom_index=False):
    """Return a NetworkX graph representation of the residue.

    :param residue: the residue from which to derive the graph
    :type residue: :class:`ost.mol.ResidueHandle` or
                   :class:`ost.mol.ResidueView`
    :param by_atom_index: Set this parameter to True if you need the nodes to
                          be labeled by atom index (within the residue).
                          Otherwise, if False, the nodes will be labeled by
                          atom names.
    :type by_atom_index: :class:`bool`
    :rtype: :class:`~networkx.classes.graph.Graph`

    Nodes are labeled with the Atom's uppercase :attr:`~ost.mol.AtomHandle.element`.
    """
    nxg = networkx.Graph()

    for atom in residue.atoms:
        nxg.add_node(atom.name, element=atom.element.upper())

    # This will list all edges twice - once for every atom of the pair.
    # But as of NetworkX 3.0 adding the same edge twice has no effect, so we're good.
    nxg.add_edges_from([(
        b.first.name,
        b.second.name) for a in residue.atoms for b in a.GetBondList()])

    if by_atom_index:
        nxg = networkx.relabel_nodes(nxg,
                                     {a: b for a, b in zip(
                                         [a.name for a in residue.atoms],
                                         range(len(residue.atoms)))},
                                     True)
    return nxg


def SCRMSD(model_ligand, target_ligand, transformation=geom.Mat4(),
           substructure_match=False, max_symmetries=1e6):
    """Calculate symmetry-corrected RMSD.

    Binding site superposition must be computed separately and passed as
    `transformation`.

    :param model_ligand: The model ligand
    :type model_ligand: :class:`ost.mol.ResidueHandle` or
                        :class:`ost.mol.ResidueView`
    :param target_ligand: The target ligand
    :type target_ligand: :class:`ost.mol.ResidueHandle` or
                         :class:`ost.mol.ResidueView`
    :param transformation: Optional transformation to apply on each atom
                           position of model_ligand.
    :type transformation: :class:`ost.geom.Mat4`
    :param substructure_match: Set this to True to allow partial target
                               ligand.
    :type substructure_match: :class:`bool`
    :param max_symmetries: If more than that many isomorphisms exist, raise
      a :class:`TooManySymmetriesError`. This can only be assessed by
      generating at least that many isomorphisms and can take some time.
    :type max_symmetries: :class:`int`
    :rtype: :class:`float`
    :raises: :class:`NoSymmetryError` when no symmetry can be found,
             :class:`DisconnectedGraphError` when ligand graph is disconnected,
             :class:`TooManySymmetriesError` when more than `max_symmetries`
             isomorphisms are found.
    """

    symmetries = _ComputeSymmetries(model_ligand, target_ligand,
                                    substructure_match=substructure_match,
                                    by_atom_index=True,
                                    max_symmetries=max_symmetries)
    return _SCRMSD_symmetries(symmetries, model_ligand, target_ligand,
                              transformation)


def _SCRMSD_symmetries(symmetries, model_ligand, target_ligand, 
                       transformation):
    """Compute SCRMSD with pre-computed symmetries. Internal. """

    best_rmsd = np.inf
    for i, (trg_sym, mdl_sym) in enumerate(symmetries):
        # Prepare Entities for RMSD
        trg_lig_rmsd_ent = mol.CreateEntity()
        trg_lig_rmsd_editor = trg_lig_rmsd_ent.EditXCS()
        trg_lig_rmsd_chain = trg_lig_rmsd_editor.InsertChain("_")
        trg_lig_rmsd_res = trg_lig_rmsd_editor.AppendResidue(trg_lig_rmsd_chain, "LIG")

        mdl_lig_rmsd_ent = mol.CreateEntity()
        mdl_lig_rmsd_editor = mdl_lig_rmsd_ent.EditXCS()
        mdl_lig_rmsd_chain = mdl_lig_rmsd_editor.InsertChain("_")
        mdl_lig_rmsd_res = mdl_lig_rmsd_editor.AppendResidue(mdl_lig_rmsd_chain, "LIG")

        for mdl_anum, trg_anum in zip(mdl_sym, trg_sym):
            # Rename model atoms according to symmetry
            trg_atom = target_ligand.atoms[trg_anum]
            mdl_atom = model_ligand.atoms[mdl_anum]
            # Add atoms in the correct order to the RMSD entities
            trg_lig_rmsd_editor.InsertAtom(trg_lig_rmsd_res, trg_atom.name, trg_atom.pos)
            mdl_lig_rmsd_editor.InsertAtom(mdl_lig_rmsd_res, mdl_atom.name, mdl_atom.pos)

        trg_lig_rmsd_editor.UpdateICS()
        mdl_lig_rmsd_editor.UpdateICS()

        rmsd = mol.alg.CalculateRMSD(mdl_lig_rmsd_ent.CreateFullView(),
                                     trg_lig_rmsd_ent.CreateFullView(),
                                     transformation)
        if rmsd < best_rmsd:
            best_rmsd = rmsd

    return best_rmsd


def _ComputeSymmetries(model_ligand, target_ligand, substructure_match=False,
                       by_atom_index=False, return_symmetries=True,
                       max_symmetries=1e6):
    """Return a list of symmetries (isomorphisms) of the model onto the target
    residues.

    :param model_ligand: The model ligand
    :type model_ligand: :class:`ost.mol.ResidueHandle` or
                        :class:`ost.mol.ResidueView`
    :param target_ligand: The target ligand
    :type target_ligand: :class:`ost.mol.ResidueHandle` or
                         :class:`ost.mol.ResidueView`
    :param substructure_match: Set this to True to allow partial ligands
                               in the reference.
    :type substructure_match: :class:`bool`
    :param by_atom_index: Set this parameter to True if you need the symmetries
                          to refer to atom index (within the residue).
                          Otherwise, if False, the symmetries refer to atom
                          names.
    :type by_atom_index: :class:`bool`
    :type return_symmetries: If Truthy, return the mappings, otherwise simply
                             return True if a mapping is found (and raise if
                             no mapping is found). This is useful to quickly
                             find out if a mapping exist without the expensive
                             step to find all the mappings.
    :type return_symmetries: :class:`bool`
    :param max_symmetries: If more than that many isomorphisms exist, raise
      a :class:`TooManySymmetriesError`. This can only be assessed by
      generating at least that many isomorphisms and can take some time.
    :type max_symmetries: :class:`int`
    :raises: :class:`NoSymmetryError` when no symmetry can be found;
             :class:`TooManySymmetriesError` when more than `max_symmetries`
             isomorphisms are found.

    """

    # Get the Graphs of the ligands
    model_graph = _ResidueToGraph(model_ligand, by_atom_index=by_atom_index)
    target_graph = _ResidueToGraph(target_ligand, by_atom_index=by_atom_index)

    if not networkx.is_connected(model_graph):
        raise DisconnectedGraphError("Disconnected graph for model ligand %s" % model_ligand)
    if not networkx.is_connected(target_graph):
        raise DisconnectedGraphError("Disconnected graph for target ligand %s" % target_ligand)

    # Note the argument order (model, target) which differs from spyrmsd.
    # This is because a subgraph of model is isomorphic to target - but not the opposite
    # as we only consider partial ligands in the reference.
    # Make sure to generate the symmetries correctly in the end
    gm = networkx.algorithms.isomorphism.GraphMatcher(
        model_graph, target_graph, node_match=lambda x, y:
        x["element"] == y["element"])
    if gm.is_isomorphic():
        if not return_symmetries:
            return True
        symmetries = []
        for i, isomorphism in enumerate(gm.isomorphisms_iter()):
            if i >= max_symmetries:
                raise TooManySymmetriesError(
                    "Too many symmetries between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
            symmetries.append((list(isomorphism.values()), list(isomorphism.keys())))
        assert len(symmetries) > 0
        LogDebug("Found %s isomorphic mappings (symmetries)" % len(symmetries))
    elif gm.subgraph_is_isomorphic() and substructure_match:
        if not return_symmetries:
            return True
        symmetries = []
        for i, isomorphism in enumerate(gm.subgraph_isomorphisms_iter()):
            if i >= max_symmetries:
                raise TooManySymmetriesError(
                    "Too many symmetries between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
            symmetries.append((list(isomorphism.values()), list(isomorphism.keys())))
        assert len(symmetries) > 0
        # Assert that all the atoms in the target are part of the substructure
        assert len(symmetries[0][0]) == len(target_ligand.atoms)
        LogDebug("Found %s subgraph isomorphisms (symmetries)" % len(symmetries))
    elif gm.subgraph_is_isomorphic():
        LogDebug("Found subgraph isomorphisms (symmetries), but"
                 " ignoring because substructure_match=False")
        raise NoSymmetryError("No symmetry between %s and %s" % (
            str(model_ligand), str(target_ligand)))
    else:
        LogDebug("Found no isomorphic mappings (symmetries)")
        raise NoSymmetryError("No symmetry between %s and %s" % (
            str(model_ligand), str(target_ligand)))

    return symmetries


class NoSymmetryError(ValueError):
    """Exception raised when no symmetry can be found.
    """
    pass


class TooManySymmetriesError(ValueError):
    """Exception raised when too many symmetries are found.
    """
    pass

class DisconnectedGraphError(Exception):
    """Exception raised when the ligand graph is disconnected.
    """
    pass


__all__ = ["LigandScorer", "SCRMSD", "NoSymmetryError", 
           "TooManySymmetriesError", "DisconnectedGraphError"]
