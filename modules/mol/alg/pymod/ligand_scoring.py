import warnings

import numpy as np
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

    * lDDT-PLI
    * Symmetry-corrected RMSD

    The class takes care to perform chain mapping and assignment (mapping) of
    model and target ligands. This assignment may differ between scores.

    It mostly expects cleaned up structures (you can use the
    :class:`~ost.mol.alg.scoring.Scorer` outputs for that). In addition,
    you probably want to remove hydrogen atoms from the structures before
    calling this function. You can do this easily with a selection::

        target_noH = target.Select("ele != H")
        model_noH = model.Select("ele != H")
        LigandScorer(model_noH, target_noH, ...)

    Make sure to remove hydrogen atoms from the ligands too.

    The class generally assumes that the
    :attr:`~ost.mol.ResidueHandle.is_ligand` property is properly set on all
    the ligand atoms, and only ligand atoms. This is typically the case for
    entities loaded from mmCIF (tested with mmCIF files from the PDB and
    SWISS-MODEL), but it will most likely not work for most entities loaded
    from PDB files.

    Unlike most of OpenStructure, this class does not assume that the ligands
    (either in the model or the target) are part of the PDB component
    dictionary. They may have arbitrary residue names. Residue names do not
    have to match between the model and the target.
    It is up to the caller to ensure that the connectivity of atoms is properly
    set before passing any ligands to this class. Ligands with improper
    connectivity will lead to bogus results.

    Note, however, that atom names should be unique within a residue (ie two
    distinct atoms cannot have the same atom name).

    This only applies to the ligand. The rest of the model and target
    structures (protein, nucleic acids) must still follow the usual rules and
    contain only residues from the compound library.


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
    :param chain_mapper: a chain mapper initialized for the target structure.
                         If None (default), a chain mapper will be initialized
                         lazily as required.
    :type chain_mapper:  :class:`ost.mol.alg.chain_mapping.ChainMapper`
    :param substructure_match: Set this to True to allow partial target ligand.
    :type substructure_match: :class:`bool`
    :param radius: Inclusion radius for the binding site. Any residue with
                   atoms within this distance of the ligand will be included
                   in the binding site.
    :type radius: :class:`float`
    :param lddt_pli_radius: lDDT inclusion radius for lDDT-PLI.
    :type lddt_pli_radius: :class:`float`
    :param lddt_bs_radius: lDDT inclusion radius for lDDT-BS.
    :type lddt_bs_radius: :class:`float`
    """
    def __init__(self, model, target, model_ligands=None, target_ligands=None,
                 resnum_alignments=False, check_resnames=True,
                 chain_mapper=None, substructure_match=False,
                 radius=4.0, lddt_pli_radius=6.0, lddt_bs_radius=10.0):

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
            raise RuntimeError("model must be of type EntityView/EntityHandle")

        # Extract ligands from target
        if target_ligands is None:
            self.target_ligands = self._extract_ligands(self.target)
        else:
            self.target_ligands = self._prepare_ligands(self.target, target, target_ligands)

        # Extract ligands from model
        if model_ligands is None:
            self.model_ligands = self._extract_ligands(self.model)
        else:
            self.model_ligands = self._prepare_ligands(self.model, model, model_ligands)

        self._chain_mapper = chain_mapper
        self.resnum_alignments = resnum_alignments
        self.check_resnames = check_resnames
        self.substructure_match = substructure_match
        self.radius = radius
        self.lddt_pli_radius = lddt_pli_radius
        self.lddt_bs_radius = lddt_bs_radius

        # scoring matrices
        self._rmsd_matrix = None
        self._rmsd_full_matrix = None
        self._lddt_pli_matrix = None
        self._lddt_pli_full_matrix = None

        # lazily computed scores
        self._rmsd = None
        self._rmsd_assignment = None
        self._rmsd_details = None
        self._lddt_pli = None
        self._lddt_pli_assignment = None
        self._lddt_pli_details = None

        # lazily precomputed variables
        self._binding_sites = {}

    @property
    def chain_mapper(self):
        """ Chain mapper object for given :attr:`target`

        :type: :class:`ost.mol.alg.chain_mapping.ChainMapper`
        """
        if self._chain_mapper is None:
            self._chain_mapper = chain_mapping.ChainMapper(self.target,
                                                           n_max_naive=1e9,
                                                           resnum_alignments=self.resnum_alignments)
        return self._chain_mapper

    @staticmethod
    def _extract_ligands(entity):
        """Extracts ligands from entity. Returns a list of residues.

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
    def _prepare_ligands(new_entity, old_entity, ligands):
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

        def _copy_residue(handle):
            """ Copy the residue handle into the new chain.
            Return the new residue handle."""
            nonlocal next_chain_num, new_editor

            # Does a residue with the same name already exist?
            already_exists = new_entity.FindResidue(handle.chain.name,
                                                    handle.number).IsValid()
            if already_exists:
                msg = "A residue number %s already exists in chain %s" % (
                    handle.number, handle.chain.name)
                raise RuntimeError(msg)

            # Instantiate the editor
            if new_editor is None:
                new_editor = new_entity.EditXCS()

            # Get or create the chain
            new_chain = new_entity.FindChain(handle.chain.name)
            if not new_chain.IsValid():
                new_chain = new_editor.InsertChain(handle.chain.name)
            # Add the residue with its original residue number
            new_res = new_editor.AppendResidue(new_chain, handle, deep=True)
            return new_res

        def _process_ligand_residue(res):
            """Copy or fetch the residue. Return the residue handle."""
            if res.entity.handle == old_entity.handle:
                # Residue is already in copied entity. We only need to grab it
                new_res = new_entity.FindResidue(res.chain.name, res.number)
                LogVerbose("Ligand residue %s already in entity" % res.handle.qualified_name)
            else:
                # Residue is not part of the entity, need to copy it first
                new_res = _copy_residue(res.handle)
                LogVerbose("Copied ligand residue %s" % res.handle.qualified_name)
            new_res.SetIsLigand(True)
            return new_res

        for ligand in ligands:
            if isinstance(ligand, mol.EntityHandle) or isinstance(ligand, mol.EntityView):
                for residue in ligand.residues:
                    new_residue = _process_ligand_residue(residue)
                    extracted_ligands.append(new_residue)
            elif isinstance(ligand, mol.ResidueHandle) or isinstance(ligand, mol.ResidueView):
                new_residue = _process_ligand_residue(ligand)
                extracted_ligands.append(new_residue)
            else:
                raise RuntimeError("Ligands should be given as Entity or Residue")

        if new_editor is not None:
            new_editor.UpdateICS()
        return extracted_ligands

    def _get_binding_sites(self, ligand, topn=100000):
        """Find representations of the binding site of *ligand* in the model.

        Ignore other ligands and waters that may be in proximity.

        :param ligand: Defines the binding site to identify.
        :type ligand: :class:`~ost.mol.ResidueHandle`
        """
        if ligand.hash_code not in self._binding_sites:

            # create view of reference binding site
            ref_residues_hashes = set()  # helper to keep track of added residues
            for ligand_at in ligand.atoms:
                close_atoms = self.target.FindWithin(ligand_at.GetPos(), self.radius)
                for close_at in close_atoms:
                    # Skip other ligands and waters.
                    # This assumes that .IsLigand() is properly set on the entity's
                    # residues.
                    ref_res = close_at.GetResidue()
                    if not (ref_res.is_ligand or
                            ref_res.chem_type == mol.ChemType.WATERS):
                        h = ref_res.handle.GetHashCode()
                        if h not in ref_residues_hashes:
                            ref_residues_hashes.add(h)

            # reason for doing that separately is to guarantee same ordering of
            # residues as in underlying entity. (Reorder by ResNum seems only
            # available on ChainHandles)
            ref_bs = self.target.CreateEmptyView()
            for ch in self.target.chains:
                for r in ch.residues:
                    if r.handle.GetHashCode() in ref_residues_hashes:
                        ref_bs.AddResidue(r, mol.ViewAddFlag.INCLUDE_ALL)

            # Find the representations
            self._binding_sites[ligand.hash_code] = self.chain_mapper.GetRepr(
                ref_bs, self.model, inclusion_radius=self.lddt_bs_radius,
                topn=topn)
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
        """"""
        # Create the matrix
        self._rmsd_full_matrix = np.empty(
            (len(self.target_ligands), len(self.model_ligands)), dtype=dict)
        self._lddt_pli_full_matrix = np.empty(
            (len(self.target_ligands), len(self.model_ligands)), dtype=dict)
        for target_i, target_ligand in enumerate(self.target_ligands):
            LogVerbose("Analyzing target ligand %s" % target_ligand)

            for binding_site in self._get_binding_sites(target_ligand):
                if len(binding_site.substructure.residues) == 0:
                    LogWarning("No residue in proximity of target ligand "
                               "%s" % str(target_ligand))
                    continue  # next binding site
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
                            by_atom_index=True)
                        LogVerbose("Ligands %s and %s symmetry match" % (
                            str(model_ligand), str(target_ligand)))
                    except NoSymmetryError:
                        # Ligands are different - skip
                        LogVerbose("No symmetry between %s and %s" % (
                            str(model_ligand), str(target_ligand)))
                        continue

                    #LogVerbose("Compute RMSD for model ligand %s" % model_ligand)
                    rmsd = SCRMSD(model_ligand, target_ligand,
                                  transformation=binding_site.transform,
                                  substructure_match=self.substructure_match)
                    self._rmsd_full_matrix[target_i, model_i] = {
                        "rmsd": rmsd,
                        "chain_mapping": binding_site.GetFlatChainMapping(),
                        "lddt_bs": binding_site.lDDT,
                        "bb_rmsd": binding_site.bb_rmsd,
                        "bs_num_res": len(binding_site.substructure.residues),
                        "bs_num_overlap_res": len(binding_site.ref_residues),
                        "target_ligand": target_ligand,
                        "model_ligand": model_ligand
                    }

                    mdl_bs_ent = self._build_binding_site_entity(
                        model_ligand, binding_site.mdl_residues, [])
                    mdl_bs_ent_ligand = mdl_bs_ent.FindResidue("_", 1)  # by definition

                    # Prepare to save the data for this target/model mapping
                    # TODO: figure out if this try/except is still needed
                    # try:
                    #     bb_rmsd = binding_site.bb_rmsd
                    # except Exception as err:
                    #     # TODO: switch to whole backbone superposition - and drop try/except
                    #     LogWarning("Can't calculate backbone RMSD: %s"
                    #                " - setting to Infinity" % str(err))
                    #     bb_rmsd = float("inf")
                    self._lddt_pli_full_matrix[target_i, model_i] = {
                        "lddt_pli": 0,
                        "lddt_pli_n_contacts": None,
                        "rmsd": rmsd,
                        # "symmetry_number": i,
                        "chain_mapping": binding_site.GetFlatChainMapping(),
                        "lddt_bs": binding_site.lDDT,
                        "bb_rmsd": binding_site.bb_rmsd,
                        "bs_num_res": len(binding_site.substructure.residues),
                        "bs_num_overlap_res": len(binding_site.ref_residues),
                        "target_ligand": target_ligand,
                        "model_ligand": model_ligand
                    }

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

                        # Save results?
                        best_lddt = self._lddt_pli_full_matrix[
                            target_i, model_i]["lddt_pli"]
                        if global_lddt > best_lddt:
                            self._lddt_pli_full_matrix[target_i, model_i].update({
                                "lddt_pli": global_lddt,
                                "lddt_pli_n_contacts": lddt_tot,
                            })

    @staticmethod
    def _find_ligand_assignment(mat1, mat2):
        """ Find the ligand assignment based on mat1

        Both mat1 and mat2 should "look" like RMSD - ie be between inf (bad)
        and 0 (good).
        """
        # We will modify mat1 and mat2, so make copies of it first
        mat1 = np.copy(mat1)
        mat2 = np.copy(mat2)
        assignments = []
        min_mat1 = mat1.min()
        while min_mat1 < np.inf:
            best_mat1 = np.argwhere(mat1 == min_mat1)
            # Multiple "best" - use mat2 to disambiguate
            if len(best_mat1) > 1:
                # Get the values of mat2 at these positions
                best_mat2_match = [mat2[tuple(x)] for x in best_mat1]
                # Find the index of the best mat2
                best_mat2_idx = np.array(best_mat2_match).argmin()
                # Now get the original indices
                max_i_trg, max_i_mdl = best_mat1[best_mat2_idx]
            else:
                max_i_trg, max_i_mdl = best_mat1[0]

            # Disable row and column
            mat1[max_i_trg, :] = np.inf
            mat1[:, max_i_mdl] = np.inf
            mat2[max_i_trg, :] = np.inf
            mat2[:, max_i_mdl] = np.inf

            # Save
            assignments.append((max_i_trg, max_i_mdl))

            # Recompute min
            min_mat1 = mat1.min()
        return assignments

    def _assign_ligands_rmsd(self):
        """Assign (map) ligands between model and target
        """
        # Transform lddt_pli to be on the scale of RMSD
        with warnings.catch_warnings():  # RuntimeWarning: divide by zero
            warnings.simplefilter("ignore")
            mat2 = np.float64(1) / self.lddt_pli_matrix

        assignments = self._find_ligand_assignment(self.rmsd_matrix, mat2)
        self._rmsd = {}
        self._rmsd_assignment = {}
        self._rmsd_details = {}
        for assignment in assignments:
            trg_idx, mdl_idx = assignment
            mdl_lig_qname = self.model_ligands[mdl_idx].qualified_name
            self._rmsd[mdl_lig_qname] = self._rmsd_full_matrix[
                trg_idx, mdl_idx]["rmsd"]
            self._rmsd_assignment[mdl_lig_qname] = self._rmsd_full_matrix[
                trg_idx, mdl_idx]["target_ligand"].qualified_name
            self._rmsd_details[mdl_lig_qname] = self._rmsd_full_matrix[
                trg_idx, mdl_idx]

    def _assign_ligands_lddt_pli(self):
        """ Assign ligands based on lDDT-PLI.
        """
        # Transform lddt_pli to be on the scale of RMSD
        with warnings.catch_warnings():  # RuntimeWarning: divide by zero
            warnings.simplefilter("ignore")
            mat1 = np.float64(1) / self.lddt_pli_matrix

        assignments = self._find_ligand_assignment(mat1, self.rmsd_matrix)
        self._lddt_pli = {}
        self._lddt_pli_assignment = {}
        self._lddt_pli_details = {}
        for assignment in assignments:
            trg_idx, mdl_idx = assignment
            mdl_lig_qname = self.model_ligands[mdl_idx].qualified_name
            self._lddt_pli[mdl_lig_qname] = self._lddt_pli_full_matrix[
                trg_idx, mdl_idx]["lddt_pli"]
            self._lddt_pli_assignment[mdl_lig_qname] = self._lddt_pli_full_matrix[
                trg_idx, mdl_idx]["target_ligand"].qualified_name
            self._lddt_pli_details[mdl_lig_qname] = self._lddt_pli_full_matrix[
                trg_idx, mdl_idx]

    @property
    def rmsd_matrix(self):
        """ Get the matrix of RMSD values.

        Target ligands are in rows, model ligands in columns.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._rmsd_full_matrix is None:
            self._compute_scores()
        if self._rmsd_matrix is None:
            # convert
            shape = self._rmsd_full_matrix.shape
            self._rmsd_matrix = np.full(shape, np.inf)
            for i, j in np.ndindex(shape):
                if self._rmsd_full_matrix[i, j] is not None:
                    self._rmsd_matrix[i, j] = self._rmsd_full_matrix[
                        i, j]["rmsd"]
        return self._rmsd_matrix

    @property
    def lddt_pli_matrix(self):
        """ Get the matrix of lDDT-PLI values.

        Target ligands are in rows, model ligands in columns.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._lddt_pli_full_matrix is None:
            self._compute_scores()
        if self._lddt_pli_matrix is None:
            # convert
            shape = self._lddt_pli_full_matrix.shape
            self._lddt_pli_matrix = np.zeros(shape)
            for i, j in np.ndindex(shape):
                if self._lddt_pli_full_matrix[i, j] is not None:
                    self._lddt_pli_matrix[i, j] = self._lddt_pli_full_matrix[
                        i, j]["lddt_pli"]
        return self._lddt_pli_matrix

    @property
    def rmsd(self):
        """Get a dictionary of RMSD score values, keyed by model ligand
        qualified names.

        :rtype: :class:`dict`
        """
        if self._rmsd is None:
            self._assign_ligands_rmsd()
        return self._rmsd

    @property
    def rmsd_assignment(self):
        """Get a dictionary of RMSD-based ligand assignment, keyed by model
        ligand qualified names. Values are the qualified names of the
        corresponding target ligand.

        :rtype: :class:`dict`
        """
        if self._rmsd_assignment is None:
            self._assign_ligands_rmsd()
        return self._rmsd_assignment

    @property
    def rmsd_details(self):
        """Get a dictionary of RMSD score details (dictionaries), keyed by
        model ligand qualified names.

        Each sub-dictionary contains the following information:

        * `rmsd`: the RMSD score value
        * `lddt_bs`: the lDDT-BS score of the binding site
        * `bs_num_res`: number of residues in the target binding site
        * `bs_num_overlap_res`: number of residues in the model overlapping
          with the target binding site.
        * `bb_rmsd`: the RMSD of the binding site backbone after superposition
        * `target_ligand`: residue handle of the target ligand
        * `model_ligand`: residue handle of the model ligand
        * `chain_mapping`: local chain mapping as a dictionary, with target
          chain name as key and model chain name as value.

        :rtype: :class:`dict`
        """
        if self._rmsd_details is None:
            self._assign_ligands_rmsd()
        return self._rmsd_details

    @property
    def lddt_pli(self):
        """Get a dictionary of lDDT-PLI score values, keyed by model ligand
        qualified names.

        :rtype: :class:`dict`
        """
        if self._lddt_pli is None:
            self._assign_ligands_lddt_pli()
        return self._lddt_pli

    @property
    def lddt_pli_assignment(self):
        """Get a dictionary of lDDT-PLI-based ligand assignment, keyed by model
        ligand qualified names. Values are the qualified names of the
        corresponding target ligand.

        :rtype: :class:`dict`
        """
        if self._lddt_pli_assignment is None:
            self._assign_ligands_lddt_pli()
        return self._lddt_pli_assignment

    @property
    def lddt_pli_details(self):
        """Get a dictionary of lDDT-PLI score details (dictionaries), keyed by
        model ligand qualified names.

        Each sub-dictionary contains the following information:

        * `lddt_pli`: the lDDT-PLI score value
        * `lddt_pli_n_contacts`: number of total contacts used in lDDT-PLI,
          summed over all thresholds. Can be divided by 8 to obtain the number
          of atomic contacts.
        * `rmsd`: the RMSD score value corresponding to the lDDT-PLI
          assignment. This may differ from the RMSD-based assignment.
        * `lddt_bs`: the lDDT-BS score of the binding site
        * `bs_num_res`: number of residues in the target binding site
        * `bs_num_overlap_res`: number of residues in the model overlapping
          with the target binding site.
        * `bb_rmsd`: the RMSD of the binding site backbone after superposition.
          Note: not used for lDDT-PLI computation.
        * `target_ligand`: residue handle of the target ligand
        * `model_ligand`: residue handle of the model ligand
        * `chain_mapping`: local chain mapping as a dictionary, with target
          chain name as key and model chain name as value.

        :rtype: :class:`dict`
        """
        if self._lddt_pli_details is None:
            self._assign_ligands_lddt_pli()
        return self._lddt_pli_details


def ResidueToGraph(residue, by_atom_index=False):
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

    Nodes are labeled with the Atom's :attr:`~ost.mol.AtomHandle.element`.
    """
    nxg = networkx.Graph()

    for atom in residue.atoms:
        nxg.add_node(atom.name, element=atom.element)

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
           substructure_match=False):
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
    :rtype: :class:`float`
    :raises: :class:`NoSymmetryError` when no symmetry can be found.
    """

    symmetries = _ComputeSymmetries(model_ligand, target_ligand,
                                    substructure_match=substructure_match,
                                    by_atom_index=True)

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
                       by_atom_index=False):
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
    :raises: :class:`NoSymmetryError` when no symmetry can be found.

    """

    # Get the Graphs of the ligands
    model_graph = ResidueToGraph(model_ligand, by_atom_index=by_atom_index)
    target_graph = ResidueToGraph(target_ligand, by_atom_index=by_atom_index)

    if not networkx.is_connected(model_graph):
        raise RuntimeError("Disconnected graph for model ligand %s" % model_ligand)
    if not networkx.is_connected(target_graph):
        raise RuntimeError("Disconnected graph for target ligand %s" % target_ligand)

    # Note the argument order (model, target) which differs from spyrmsd.
    # This is because a subgraph of model is isomorphic to target - but not the opposite
    # as we only consider partial ligands in the reference.
    # Make sure to generate the symmetries correctly in the end
    gm = networkx.algorithms.isomorphism.GraphMatcher(
        model_graph, target_graph, node_match=lambda x, y:
        x["element"] == y["element"])
    if gm.is_isomorphic():
        symmetries = [
            (list(isomorphism.values()), list(isomorphism.keys()))
            for isomorphism in gm.isomorphisms_iter()]
        assert len(symmetries) > 0
        LogDebug("Found %s isomorphic mappings (symmetries)" % len(symmetries))
    elif gm.subgraph_is_isomorphic() and substructure_match:
        symmetries = [(list(isomorphism.values()), list(isomorphism.keys())) for isomorphism in
                      gm.subgraph_isomorphisms_iter()]
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


class NoSymmetryError(Exception):
    """Exception to be raised when no symmetry can be found.
    """
    pass


__all__ = ["LigandScorer", "ResidueToGraph", "SCRMSD", "NoSymmetryError"]
