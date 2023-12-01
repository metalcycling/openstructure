import os
from ost import mol
from ost import seq
from ost import io
from ost import conop
from ost import settings
from ost import geom
from ost.mol.alg import lddt
from ost.mol.alg import qsscore
from ost.mol.alg import chain_mapping
from ost.mol.alg import stereochemistry
from ost.mol.alg import dockq
from ost.mol.alg.lddt import lDDTScorer
from ost.mol.alg.qsscore import QSScorer
from ost.mol.alg.contact_score import ContactScorer
from ost.mol.alg import Molck, MolckSettings
from ost import bindings
from ost.bindings import cadscore
from ost.bindings import tmtools
import numpy as np

class lDDTBSScorer:
    """Scorer specific for a reference/model pair

    Finds best possible binding site representation of reference in model given
    lDDT score. Uses :class:`ost.mol.alg.chain_mapping.ChainMapper` to deal with
    chain mapping.

    :param reference: Reference structure
    :type reference: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param model: Model structure
    :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param residue_number_alignment: Passed to ChainMapper constructor
    :type residue_number_alignment: :class:`bool`
    """
    def __init__(self, reference, model,
                 residue_number_alignment=False):
        self.chain_mapper = chain_mapping.ChainMapper(reference,
            resnum_alignments=residue_number_alignment)
        self.ref = self.chain_mapper.target
        self.mdl = model

    def ScoreBS(self, ligand, radius = 4.0, lddt_radius=10.0):
        """Computes binding site lDDT score given *ligand*. Best possible
        binding site representation is selected by lDDT but other scores such as
        CA based RMSD and GDT are computed too and returned.

        :param ligand: Defines the scored binding site, i.e. provides positions
                       to perform proximity search
        :type ligand: r'((Residue)|(Chain)|(Entity))((View)|(Handle))'
        :param radius: Reference residues with any atom position within *radius*
                       of *ligand* consitute the scored binding site
        :type radius: :class:`float`
        :param lddt_radius: Passed as *inclusion_radius* to
                            :class:`ost.mol.alg.lddt.lDDTScorer`
        :type lddt_radius: :class:`float`
        :returns: Object of type :class:`ost.mol.alg.chain_mapping.ReprResult`
                  containing all atom lDDT score and mapping information.
                  None if no representation could be found.
        """

        # create view of reference binding site
        ref_residues_hashes = set() # helper to keep track of added residues
        for ligand_at in ligand.atoms:
            close_atoms = self.ref.FindWithin(ligand_at.GetPos(), radius)
            for close_at in close_atoms:
                ref_res = close_at.GetResidue()
                h = ref_res.handle.GetHashCode()
                if h not in ref_residues_hashes:
                    ref_residues_hashes.add(h)

        # reason for doing that separately is to guarantee same ordering of
        # residues as in underlying entity. (Reorder by ResNum seems only
        # available on ChainHandles)
        ref_bs = self.ref.CreateEmptyView()
        for ch in self.ref.chains:
            for r in ch.residues:
                if r.handle.GetHashCode() in ref_residues_hashes:
                    ref_bs.AddResidue(r, mol.ViewAddFlag.INCLUDE_ALL)

        # gogogo
        bs_repr = self.chain_mapper.GetRepr(ref_bs, self.mdl,
                                            inclusion_radius = lddt_radius)
        if len(bs_repr) >= 1:
            return bs_repr[0]
        else:
            return None


class Scorer:
    """ Helper class to access the various scores available from ost.mol.alg

    Deals with structure cleanup, chain mapping, interface identification etc.
    Intermediate results are available as attributes.

    :param model: Model structure - a deep copy is available as :attr:`model`.
                  Additionally, :func:`ost.mol.alg.Molck` using *molck_settings*
                  is applied.
    :type model: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param target: Target structure - a deep copy is available as :attr:`target`.
                  Additionally, :func:`ost.mol.alg.Molck` using *molck_settings*
                  is applied.
    :type target: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param resnum_alignments: Whether alignments between chemically equivalent
                              chains in *model* and *target* can be computed
                              based on residue numbers. This can be assumed in
                              benchmarking setups such as CAMEO/CASP.
    :type resnum_alignments: :class:`bool`
    :param molck_settings: Settings used for Molck on *model* and *target*, if
                           set to None, a default object is constructed by
                           setting everything except rm_zero_occ_atoms and
                           colored to True in
                           :class:`ost.mol.alg.MolckSettings` constructor.
    :type molck_settings: :class:`ost.mol.alg.MolckSettings`
    :param cad_score_exec: Explicit path to voronota-cadscore executable from
                           voronota installation from 
                           https://github.com/kliment-olechnovic/voronota. If
                           not given, voronota-cadscore must be in PATH if any
                           of the CAD score related attributes is requested.
    :type cad_score_exec: :class:`str`
    :param custom_mapping: Provide custom chain mapping between *model* and
                           *target*. Dictionary with target chain names as key
                           and model chain names as value.
    :type custom_mapping: :class:`dict`
    :param usalign_exec: Explicit path to USalign executable used to compute
                         TM-score. If not given, TM-score will be computed
                         with OpenStructure internal copy of USalign code.
    :type usalign_exec: :class:`str`
    :param lddt_no_stereochecks: Whether to compute lDDT without stereochemistry
                                checks
    :type lddt_no_stereochecks: :class:`bool`
    :param n_max_naive: Parameter for chain mapping. If the number of possible
                        mappings is <= *n_max_naive*, the full
                        mapping solution space is enumerated to find the
                        the optimum. A heuristic is used otherwise. The default
                        of 40320 corresponds to an octamer (8! = 40320).
                        A structure with stoichiometry A6B2 would be
                        6!*2! = 1440 etc.
    :type n_max_naive: :class:`int`
    :param oum: Override USalign Mapping. Inject mapping of :class:`Scorer`
                object into USalign to compute TM-score. Experimental feature
                with limitations.
    :type oum: :class:`bool`
    """
    def __init__(self, model, target, resnum_alignments=False,
                 molck_settings = None, cad_score_exec = None,
                 custom_mapping=None, usalign_exec = None,
                 lddt_no_stereochecks=False, n_max_naive=40320,
                 oum=False):

        self._target_orig = target
        self._model_orig = model

        if isinstance(self._model_orig, mol.EntityView):
            self._model = mol.CreateEntityFromView(self._model_orig, False)
        else:
            self._model = self._model_orig.Copy()

        if isinstance(self._target_orig, mol.EntityView):
            self._target = mol.CreateEntityFromView(self._target_orig, False)
        else:
            self._target = self._target_orig.Copy()

        if molck_settings is None:
            molck_settings = MolckSettings(rm_unk_atoms=True,
                                           rm_non_std=False,
                                           rm_hyd_atoms=True,
                                           rm_oxt_atoms=True,
                                           rm_zero_occ_atoms=False,
                                           colored=False,
                                           map_nonstd_res=True,
                                           assign_elem=True)
        Molck(self._model, conop.GetDefaultLib(), molck_settings)
        Molck(self._target, conop.GetDefaultLib(), molck_settings)
        self._model = model.Select("peptide=True or nucleotide=True")
        self._target = target.Select("peptide=True or nucleotide=True")

        # catch models which have empty chain names
        for ch in self._model.chains:
            if ch.GetName().strip() == "":
                raise RuntimeError("Model chains must have valid chain names")
        
        # catch targets which have empty chain names
        for ch in self._target.chains:
            if ch.GetName().strip() == "":
                raise RuntimeError("Target chains must have valid chain names")

        if resnum_alignments:
            # In case of resnum_alignments, we have some requirements on 
            # residue numbers in the chain mapping: 1) no ins codes 2) strictly
            # increasing residue numbers.
            for ch in self._model.chains:
                ins_codes = [r.GetNumber().GetInsCode() for r in ch.residues]
                if len(set(ins_codes)) != 1 or ins_codes[0] != '\0':
                    raise RuntimeError("Residue numbers in each model chain "
                                       "must not contain insertion codes if "
                                       "resnum_alignments are enabled")
                nums = [r.GetNumber().GetNum() for r in ch.residues]
                if not all(i < j for i, j in zip(nums, nums[1:])):
                    raise RuntimeError("Residue numbers in each model chain "
                                       "must be strictly increasing if "
                                       "resnum_alignments are enabled")

            for ch in self._target.chains:
                ins_codes = [r.GetNumber().GetInsCode() for r in ch.residues]
                if len(set(ins_codes)) != 1 or ins_codes[0] != '\0':
                    raise RuntimeError("Residue numbers in each target chain "
                                       "must not contain insertion codes if "
                                       "resnum_alignments are enabled")
                nums = [r.GetNumber().GetNum() for r in ch.residues]
                if not all(i < j for i, j in zip(nums, nums[1:])):
                    raise RuntimeError("Residue numbers in each target chain "
                                       "must be strictly increasing if "
                                       "resnum_alignments are enabled")

        if usalign_exec is not None:
            if not os.path.exists(usalign_exec):
                raise RuntimeError(f"USalign exec ({usalign_exec}) "
                                   f"not found")
            if not os.access(usalign_exec, os.X_OK):
                raise RuntimeError(f"USalign exec ({usalign_exec}) "
                                   f"is not executable")

        self.resnum_alignments = resnum_alignments
        self.cad_score_exec = cad_score_exec
        self.usalign_exec = usalign_exec
        self.lddt_no_stereochecks = lddt_no_stereochecks
        self.n_max_naive = n_max_naive
        self.oum = oum

        # lazily evaluated attributes
        self._stereochecked_model = None
        self._stereochecked_target = None
        self._model_clashes = None
        self._model_bad_bonds = None
        self._model_bad_angles = None
        self._target_clashes = None
        self._target_bad_bonds = None
        self._target_bad_angles = None
        self._chain_mapper = None
        self._mapping = None
        self._model_interface_residues = None
        self._target_interface_residues = None
        self._aln = None
        self._stereochecked_aln = None
        self._pepnuc_aln = None

        # lazily constructed scorer objects
        self._lddt_scorer = None
        self._bb_lddt_scorer = None
        self._qs_scorer = None
        self._contact_scorer = None

        # lazily computed scores
        self._lddt = None
        self._local_lddt = None
        self._bb_lddt = None
        self._bb_local_lddt = None

        self._qs_global = None
        self._qs_best = None
        self._qs_target_interfaces = None
        self._qs_model_interfaces = None
        self._qs_interfaces = None
        self._per_interface_qs_global = None
        self._per_interface_qs_best = None

        self._contact_target_interfaces = None
        self._contact_model_interfaces = None
        self._native_contacts = None
        self._model_contacts = None
        self._ics_precision = None
        self._ics_recall = None
        self._ics = None
        self._per_interface_ics_precision = None
        self._per_interface_ics_recall = None
        self._per_interface_ics = None
        self._ips_precision = None
        self._ips_recall = None
        self._ips = None
        self._per_interface_ics_precision = None
        self._per_interface_ics_recall = None
        self._per_interface_ics = None

        self._dockq_target_interfaces = None
        self._dockq_interfaces = None
        self._fnat = None
        self._fnonnat = None
        self._irmsd = None
        self._lrmsd = None
        self._nnat = None
        self._nmdl = None
        self._dockq_scores = None
        self._dockq_ave = None
        self._dockq_wave = None
        self._dockq_ave_full = None
        self._dockq_wave_full = None

        self._mapped_target_pos = None
        self._mapped_model_pos = None
        self._transformed_mapped_model_pos = None
        self._n_target_not_mapped = None
        self._transform = None
        self._gdtts = None
        self._gdtha = None
        self._rmsd = None

        self._cad_score = None
        self._local_cad_score = None

        self._patch_qs = None
        self._patch_dockq = None

        self._tm_score = None
        self._usalign_mapping = None

        if custom_mapping is not None:
            self._set_custom_mapping(custom_mapping)

    @property
    def model(self):
        """ Model with Molck cleanup

        :type: :class:`ost.mol.EntityHandle`
        """
        return self._model

    @property
    def model_orig(self):
        """ The original model passed at object construction

        :type: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
        """
        return self._model_orig

    @property
    def target(self):
        """ Target with Molck cleanup

        :type: :class:`ost.mol.EntityHandle`
        """
        return self._target

    @property
    def target_orig(self):
        """ The original target passed at object construction

        :type: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
        """
        return self._target_orig

    @property
    def aln(self):
        """ Alignments of :attr:`model`/:attr:`target` chains

        Alignments for each pair of chains mapped in :attr:`mapping`.
        First sequence is target sequence, second sequence the model sequence.

        :type: :class:`list` of :class:`ost.seq.AlignmentHandle`
        """
        if self._aln is None:
            self._compute_aln()
        return self._aln

    @property
    def stereochecked_aln(self):
        """ Stereochecked equivalent of :attr:`aln`

        The alignments may differ, as stereochecks potentially remove residues

        :type: :class:`list` of :class:`ost.seq.AlignmentHandle`
        """
        if self._stereochecked_aln is None:
            self._compute_stereochecked_aln()
        return self._stereochecked_aln

    @property
    def pepnuc_aln(self):
        """ Alignments of :attr:`model_orig`/:attr:`target_orig` chains

        Selects for peptide and nucleotide residues before sequence
        extraction. Includes residues that would be removed by molck in
        structure preprocessing.

        :type: :class:`list` of :class:`ost.seq.AlignmentHandle`
        """
        if self._pepnuc_aln is None:
            self._compute_pepnuc_aln()
        return self._pepnuc_aln

    @property
    def stereochecked_model(self):
        """ View of :attr:`~model` that has stereochemistry checks applied

        First, a selection for peptide/nucleotide residues is performed,
        secondly peptide sidechains with stereochemical irregularities are
        removed (full residue if backbone atoms are involved). Irregularities
        are clashes or bond lengths/angles more than 12 standard deviations
        from expected values.

        :type: :class:`ost.mol.EntityView`
        """
        if self._stereochecked_model is None:
            self._do_stereochecks()
        return self._stereochecked_model

    @property
    def model_clashes(self):
        """ Clashing model atoms

        :type: :class:`list` of :class:`ost.mol.alg.stereochemistry.ClashInfo`
        """
        if self._model_clashes is None:
            self._do_stereochecks()
        return self._model_clashes

    @property
    def model_bad_bonds(self):
        """ Model bonds with unexpected stereochemistry

        :type: :class:`list` of
               :class:`ost.mol.alg.stereochemistry.BondViolationInfo`
        """
        if self._model_bad_bonds is None:
            self._do_stereochecks()
        return self._model_bad_bonds

    @property
    def model_bad_angles(self):
        """ Model angles with unexpected stereochemistry

        :type: :class:`list` of
               :class:`ost.mol.alg.stereochemistry.AngleViolationInfo`
        """
        if self._model_bad_angles is None:
            self._do_stereochecks()
        return self._model_bad_angles

    @property
    def stereochecked_target(self):
        """ Same as :attr:`~stereochecked_model` for :attr:`~target`

        :type: :class:`ost.mol.EntityView`
        """
        if self._stereochecked_target is None:
            self._do_stereochecks()
        return self._stereochecked_target

    @property
    def target_clashes(self):
        """ Clashing target atoms

        :type: :class:`list` of :class:`ost.mol.alg.stereochemistry.ClashInfo`
        """
        if self._target_clashes is None:
            self._do_stereochecks()
        return self._target_clashes

    @property
    def target_bad_bonds(self):
        """ Target bonds with unexpected stereochemistry

        :type: :class:`list` of
               :class:`ost.mol.alg.stereochemistry.BondViolationInfo`
        """
        if self._target_bad_bonds is None:
            self._do_stereochecks()
        return self._target_bad_bonds

    @property
    def target_bad_angles(self):
        """ Target angles with unexpected stereochemistry

        :type: :class:`list` of
               :class:`ost.mol.alg.stereochemistry.AngleViolationInfo`
        """
        if self._target_bad_angles is None:
            self._do_stereochecks()
        return self._target_bad_angles

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

    @property
    def mapping(self):
        """ Full chain mapping result for :attr:`target`/:attr:`model`

        :type: :class:`ost.mol.alg.chain_mapping.MappingResult` 
        """
        if self._mapping is None:
            self._mapping = \
            self.chain_mapper.GetMapping(self.model,
                                         n_max_naive = self.n_max_naive)
        return self._mapping

    @property
    def model_interface_residues(self):
        """ Interface residues in :attr:`~model`

        Thats all residues having a contact with at least one residue from
        another chain (CB-CB distance <= 8A, CA in case of Glycine)

        :type: :class:`dict` with chain names as key and and :class:`list`
                with residue numbers of the respective interface residues.
        """
        if self._model_interface_residues is None:
            self._model_interface_residues = \
            self._get_interface_residues(self.model)
        return self._model_interface_residues

    @property
    def target_interface_residues(self):
        """ Same as :attr:`~model_interface_residues` for :attr:`~target`

        :type: :class:`dict` with chain names as key and and :class:`list`
                with residue numbers of the respective interface residues.
        """
        if self._target_interface_residues is None:
            self._target_interface_residues = \
            self._get_interface_residues(self.target)
        return self._target_interface_residues

    @property
    def lddt_scorer(self):
        """ lDDT scorer for :attr:`~stereochecked_target` (default parameters)

        :type: :class:`ost.mol.alg.lddt.lDDTScorer`
        """
        if self._lddt_scorer is None:
            if self.lddt_no_stereochecks:
                self._lddt_scorer = lDDTScorer(self.target)
            else:
                self._lddt_scorer = lDDTScorer(self.stereochecked_target)
        return self._lddt_scorer

    @property
    def bb_lddt_scorer(self):
        """ Backbone only lDDT scorer for :attr:`~target`

        No stereochecks applied for bb only lDDT which considers CA atoms
        for peptides and C3' atoms for nucleotides.

        :type: :class:`ost.mol.alg.lddt.lDDTScorer`
        """
        if self._bb_lddt_scorer is None:
            self._bb_lddt_scorer = lDDTScorer(self.target, bb_only=True)
        return self._bb_lddt_scorer

    @property
    def qs_scorer(self):
        """ QS scorer constructed from :attr:`~mapping`

        The scorer object is constructed with default parameters and relates to
        :attr:`~model` and :attr:`~target` (no stereochecks).

        :type: :class:`ost.mol.alg.qsscore.QSScorer`
        """
        if self._qs_scorer is None:
            self._qs_scorer = QSScorer.FromMappingResult(self.mapping)
        return self._qs_scorer

    @property
    def contact_scorer(self):
        if self._contact_scorer is None:
            self._contact_scorer = ContactScorer.FromMappingResult(self.mapping)
        return self._contact_scorer
    

    @property
    def lddt(self):
        """ Global lDDT score in range [0.0, 1.0]

        Computed based on :attr:`~stereochecked_model`. In case of oligomers,
        :attr:`~mapping` is used.

        :type: :class:`float`
        """
        if self._lddt is None:
            self._compute_lddt()
        return self._lddt
    
    @property
    def local_lddt(self):
        """ Per residue lDDT scores in range [0.0, 1.0]

        Computed based on :attr:`~stereochecked_model` but scores for all 
        residues in :attr:`~model` are reported. If a residue has been removed
        by stereochemistry checks, the respective score is set to 0.0. If a
        residue is not covered by the target or is in a chain skipped by the
        chain mapping procedure (happens for super short chains), the respective
        score is set to None. In case of oligomers, :attr:`~mapping` is used.

        :type: :class:`dict`
        """
        if self._local_lddt is None:
            self._compute_lddt()
        return self._local_lddt

    @property
    def bb_lddt(self):
        """ Backbone only global lDDT score in range [0.0, 1.0]

        Computed based on :attr:`~model` on backbone atoms only. This is CA for
        peptides and C3' for nucleotides. No stereochecks are performed. In case
        of oligomers, :attr:`~mapping` is used.

        :type: :class:`float`
        """
        if self._bb_lddt is None:
            self._compute_bb_lddt()
        return self._bb_lddt
    
    @property
    def bb_local_lddt(self):
        """ Backbone only per residue lDDT scores in range [0.0, 1.0]

        Computed based on :attr:`~model` on backbone atoms only. This is CA for
        peptides and C3' for nucleotides. No stereochecks are performed. If a
        residue is not covered by the target or is in a chain skipped by the
        chain mapping procedure (happens for super short chains), the respective
        score is set to None. In case of oligomers, :attr:`~mapping` is used.

        :type: :class:`dict`
        """
        if self._bb_local_lddt is None:
            self._compute_bb_lddt()
        return self._bb_local_lddt

    @property
    def qs_global(self):
        """  Global QS-score

        Computed based on :attr:`model` using :attr:`mapping`

        :type: :class:`float`
        """
        if self._qs_global is None:
            self._compute_qs()
        return self._qs_global

    @property
    def qs_best(self):
        """  Global QS-score - only computed on aligned residues

        Computed based on :attr:`model` using :attr:`mapping`. The QS-score
        computation only considers contacts between residues with a mapping
        between target and model. As a result, the score won't be lowered in
        case of additional chains/residues in any of the structures.

        :type: :class:`float`
        """
        if self._qs_best is None:
            self._compute_qs()
        return self._qs_best

    @property
    def qs_target_interfaces(self):
        """ Interfaces in :attr:`~target` with non-zero contribution to
        :attr:`~qs_global`/:attr:`~qs_best`

        Chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each:
               (trg_ch1, trg_ch2)
        """
        if self._qs_target_interfaces is None:
            self._qs_target_interfaces = self.qs_scorer.qsent1.interacting_chains
            self._qs_target_interfaces = \
            [(min(x[0],x[1]), max(x[0],x[1])) for x in self._qs_target_interfaces]
        return self._qs_target_interfaces

    @property
    def qs_model_interfaces(self):
        """ Interfaces in :attr:`~model` with non-zero contribution to
        :attr:`~qs_global`/:attr:`~qs_best`

        Chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each:
               (mdl_ch1, mdl_ch2)
        """
        if self._qs_model_interfaces is None:
            self._qs_model_interfaces = self.qs_scorer.qsent2.interacting_chains
            self._qs_model_interfaces = \
            [(min(x[0],x[1]), max(x[0],x[1])) for x in self._qs_model_interfaces]

        return self._qs_model_interfaces

    @property
    def qs_interfaces(self):
        """ Interfaces in :attr:`~qs_target_interfaces` that can be mapped
        to :attr:`~model`.

        Target chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 4 elements each:
               (trg_ch1, trg_ch2, mdl_ch1, mdl_ch2)
        """
        if self._qs_interfaces is None:
            self._qs_interfaces = list()
            flat_mapping = self.mapping.GetFlatMapping()
            for i in self.qs_target_interfaces:
                if i[0] in flat_mapping and i[1] in flat_mapping:
                    self._qs_interfaces.append((i[0], i[1],
                                                flat_mapping[i[0]],
                                                flat_mapping[i[1]]))
        return self._qs_interfaces
    
    @property
    def per_interface_qs_global(self):
        """ QS-score for each interface in :attr:`qs_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_qs_global is None:
            self._compute_per_interface_qs_scores()
        return self._per_interface_qs_global
    
    @property
    def per_interface_qs_best(self):
        """ QS-score for each interface in :attr:`qs_interfaces`

        Only computed on aligned residues

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_qs_best is None:
            self._compute_per_interface_qs_scores()
        return self._per_interface_qs_best
    
    @property
    def native_contacts(self):
        """ Native contacts

        A contact is a pair or residues from distinct chains that have
        a minimal heavy atom distance < 5A. Contacts are specified as
        :class:`tuple` with two strings in format:
        <cname>.<rnum>.<ins_code>

        :type: :class:`list` of :class:`tuple`
        """
        if self._native_contacts is None:
            self._native_contacts = self.contact_scorer.cent1.hr_contacts
        return self._native_contacts

    @property
    def model_contacts(self):
        """ Same for model
        """
        if self._model_contacts is None:
            self._model_contacts = self.contact_scorer.cent2.hr_contacts
        return self._model_contacts

    @property
    def contact_target_interfaces(self):
        """ Interfaces in :class:`target` which have at least one contact

        Contact as defined in :attr:`~native_contacts`,
        chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each
               (trg_ch1, trg_ch2)
        """
        if self._contact_target_interfaces is None:
            tmp = self.contact_scorer.cent1.interacting_chains
            tmp = [(min(x[0],x[1]), max(x[0],x[1])) for x in tmp]
            self._contact_target_interfaces = tmp
        return self._contact_target_interfaces

    @property
    def contact_model_interfaces(self):
        """ Interfaces in :class:`model` which have at least one contact

        Contact as defined in :attr:`~native_contacts`,
        chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each
               (mdl_ch1, mdl_ch2)
        """
        if self._contact_model_interfaces is None:
            tmp = self.contact_scorer.cent2.interacting_chains
            tmp = [(min(x[0],x[1]), max(x[0],x[1])) for x in tmp]
            self._contact_model_interfaces = tmp
        return self._contact_model_interfaces

    @property
    def ics_precision(self):
        """ Fraction of model contacts that are also present in target

        :type: :class:`float`
        """
        if self._ics_precision is None:
            self._compute_ics_scores()
        return self._ics_precision
    
    @property
    def ics_recall(self):
        """ Fraction of target contacts that are correctly reproduced in model

        :type: :class:`float`
        """
        if self._ics_recall is None:
            self._compute_ics_scores()
        return self._ics_recall

    @property
    def ics(self):
        """ ICS (Interface Contact Similarity) score

        Combination of :attr:`~ics_precision` and :attr:`~ics_recall`
        using the F1-measure

        :type: :class:`float`
        """
        if self._ics is None:
            self._compute_ics_scores()
        return self._ics

    @property
    def per_interface_ics_precision(self):
        """ Per-interface ICS precision

        :attr:`~ics_precision` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ics_precision is None:
            self._compute_ics_scores()
        return self._per_interface_ics_precision


    @property
    def per_interface_ics_recall(self):
        """ Per-interface ICS recall

        :attr:`~ics_recall` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ics_recall is None:
            self._compute_ics_scores()
        return self._per_interface_ics_recall

    @property
    def per_interface_ics(self):
        """ Per-interface ICS (Interface Contact Similarity) score

        :attr:`~ics` for each interface in 
        :attr:`~contact_target_interfaces`

        :type: :class:`float`
        """

        if self._per_interface_ics is None:
            self._compute_ics_scores()
        return self._per_interface_ics
    

    @property
    def ips_precision(self):
        """ Fraction of model interface residues that are also interface
        residues in target

        :type: :class:`float`
        """
        if self._ips_precision is None:
            self._compute_ips_scores()
        return self._ips_precision
    
    @property
    def ips_recall(self):
        """ Fraction of target interface residues that are also interface
        residues in model

        :type: :class:`float`
        """
        if self._ips_recall is None:
            self._compute_ips_scores()
        return self._ips_recall

    @property
    def ips(self):
        """ IPS (Interface Patch Similarity) score

        Jaccard coefficient of interface residues in target and their mapped
        counterparts in model

        :type: :class:`float`
        """
        if self._ips is None:
            self._compute_ips_scores()
        return self._ips

    @property
    def per_interface_ips_precision(self):
        """ Per-interface IPS precision

        :attr:`~ips_precision` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ips_precision is None:
            self._compute_ips_scores()
        return self._per_interface_ips_precision


    @property
    def per_interface_ips_recall(self):
        """ Per-interface IPS recall

        :attr:`~ips_recall` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ics_recall is None:
            self._compute_ips_scores()
        return self._per_interface_ips_recall

    @property
    def per_interface_ips(self):
        """ Per-interface IPS (Interface Patch Similarity) score

        :attr:`~ips` for each interface in 
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """

        if self._per_interface_ips is None:
            self._compute_ips_scores()
        return self._per_interface_ips

    @property
    def dockq_target_interfaces(self):
        """ Interfaces in :attr:`target` that are relevant for DockQ

        In principle a subset of :attr:`~contact_target_interfaces` that only
        contains peptide sequences. Chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each:
               (trg_ch1, trg_ch2)
        """
        if self._dockq_target_interfaces is None:
            self._dockq_target_interfaces = list()
            pep_seqs = set([s.GetName() for s in self.chain_mapper.polypep_seqs])
            for interface in self.contact_target_interfaces:
                if interface[0] in pep_seqs and interface[1] in pep_seqs:
                    self._dockq_target_interfaces.append(interface)
        return self._dockq_target_interfaces

    @property
    def dockq_interfaces(self):
        """ Interfaces in :attr:`dockq_target_interfaces` that can be mapped
        to model

        Target chain names are lexicographically sorted

        :type: :class:`list` of :class:`tuple` with 4 elements each:
               (trg_ch1, trg_ch2, mdl_ch1, mdl_ch2)
        """
        if self._dockq_interfaces is None:
            self._dockq_interfaces = list()
            flat_mapping = self.mapping.GetFlatMapping()
            for i in self.dockq_target_interfaces:
                if i[0] in flat_mapping and i[1] in flat_mapping:
                    self._dockq_interfaces.append((i[0], i[1],
                                                   flat_mapping[i[0]],
                                                   flat_mapping[i[1]]))
        return self._dockq_interfaces
    
    @property
    def dockq_scores(self):
        """ DockQ scores for interfaces in :attr:`~dockq_interfaces` 

        :class:`list` of :class:`float`
        """
        if self._dockq_scores is None:
            self._compute_dockq_scores()
        return self._dockq_scores

    @property
    def fnat(self):
        """ fnat scores for interfaces in :attr:`~dockq_interfaces` 

        fnat: Fraction of native contacts that are also present in model

        :class:`list` of :class:`float`
        """
        if self._fnat is None:
            self._compute_dockq_scores()
        return self._fnat

    @property
    def nnat(self):
        """ N native contacts for interfaces in :attr:`~dockq_interfaces` 

        :class:`list` of :class:`int`
        """
        if self._nnat is None:
            self._compute_dockq_scores()
        return self._nnat

    @property
    def nmdl(self):
        """ N model contacts for interfaces in :attr:`~dockq_interfaces` 

        :class:`list` of :class:`int`
        """
        if self._nmdl is None:
            self._compute_dockq_scores()
        return self._nmdl

    @property
    def fnonnat(self):
        """ fnonnat scores for interfaces in :attr:`~dockq_interfaces` 

        fnat: Fraction of model contacts that are not present in target

        :class:`list` of :class:`float`
        """
        if self._fnonnat is None:
            self._compute_dockq_scores()
        return self._fnonnat

    @property
    def irmsd(self):
        """ irmsd scores for interfaces in :attr:`~dockq_interfaces` 

        irmsd: RMSD of interface (RMSD computed on N, CA, C, O atoms) which
        consists of each residue that has at least one heavy atom within 10A of
        other chain.

        :class:`list` of :class:`float`
        """
        if self._irmsd is None:
            self._compute_dockq_scores()
        return self._irmsd

    @property
    def lrmsd(self):
        """ lrmsd scores for interfaces in :attr:`~dockq_interfaces` 

        lrmsd: The interfaces are superposed based on the receptor (rigid
        min RMSD superposition) and RMSD for the ligand is reported.
        Superposition and RMSD are based on N, CA, C and O positions,
        receptor is the chain contributing to the interface with more
        residues in total.

        :class:`list` of :class:`float`
        """
        if self._lrmsd is None:
            self._compute_dockq_scores()
        return self._lrmsd
        
    @property
    def dockq_ave(self):
        """ Average of DockQ scores in :attr:`dockq_scores`

        In its original implementation, DockQ only operates on single
        interfaces. Thus the requirement to combine scores for higher order
        oligomers.

        :type: :class:`float`
        """
        if self._dockq_ave is None:
            self._compute_dockq_scores()
        return self._dockq_ave
    
    @property
    def dockq_wave(self):
        """ Same as :attr:`dockq_ave`, weighted by native contacts

        :type: :class:`float`
        """
        if self._dockq_wave is None:
            self._compute_dockq_scores()
        return self._dockq_wave
        
    @property
    def dockq_ave_full(self):
        """ Same as :attr:`~dockq_ave` but penalizing for missing interfaces

        Interfaces that are not covered in model are added as 0.0
        in average computation.

        :type: :class:`float`
        """
        if self._dockq_ave_full is None:
            self._compute_dockq_scores()
        return self._dockq_ave_full
    
    @property
    def dockq_wave_full(self):
        """ Same as :attr:`~dockq_ave_full`, but weighted

        Interfaces that are not covered in model are added as 0.0 in
        average computations and the respective weights are derived from
        number of contacts in respective target interface. 
        """
        if self._dockq_wave_full is None:
            self._compute_dockq_scores()
        return self._dockq_wave_full

    @property
    def mapped_target_pos(self):
        """ Mapped representative positions in target

        Thats CA positions for peptide residues and C3' positions for
        nucleotides. Has same length as :attr:`~mapped_model_pos` and mapping
        is based on :attr:`~mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._mapped_target_pos is None:
            self._extract_mapped_pos()
        return self._mapped_target_pos

    @property
    def mapped_model_pos(self):
        """ Mapped representative positions in model

        Thats CA positions for peptide residues and C3' positions for
        nucleotides. Has same length as :attr:`~mapped_target_pos` and mapping
        is based on :attr:`~mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._mapped_model_pos is None:
            self._extract_mapped_pos()
        return self._mapped_model_pos

    @property
    def transformed_mapped_model_pos(self):
        """ :attr:`~mapped_model_pos` with :attr:`~transform` applied

        :type: :class:`ost.geom.Vec3List`
        """
        if self._transformed_mapped_model_pos is None:
            self._transformed_mapped_model_pos = \
            geom.Vec3List(self.mapped_model_pos)
            self._transformed_mapped_model_pos.ApplyTransform(self.transform)
        return self._transformed_mapped_model_pos

    @property
    def n_target_not_mapped(self):
        """ Number of target residues which have no mapping to model

        :type: :class:`int`
        """
        if self._n_target_not_mapped is None:
            self._extract_mapped_pos()
        return self._n_target_not_mapped

    @property
    def transform(self):
        """ Transform: :attr:`~mapped_model_pos` onto :attr:`~mapped_target_pos`

        Computed using Kabsch minimal rmsd algorithm

        :type: :class:`ost.geom.Mat4`
        """
        if self._transform is None:
            try:
                res = mol.alg.SuperposeSVD(self.mapped_model_pos,
                                           self.mapped_target_pos)
                self._transform = res.transformation
            except:
                self._transform = geom.Mat4()
        return self._transform

    @property
    def gdtts(self):
        """ GDT with thresholds: 8.0A, 4.0A, 2.0A and 1.0A

        Computed on :attr:`~transformed_mapped_model_pos` and
        :attr:`mapped_target_pos`

        :type: :class:`float`
        """
        if self._gdtts is None:
            n = \
            self.mapped_target_pos.GetGDTTS(self.transformed_mapped_model_pos,
                                            norm=False)
            n_full = 4*len(self.mapped_target_pos) + 4*self.n_target_not_mapped
            if n_full > 0:
                self._gdtts = float(n) / n_full
            else:
                self._gdtts = 0.0
        return self._gdtts

    @property
    def gdtha(self):
        """ GDT with thresholds: 4.0A, 2.0A, 1.0A and 0.5A

        Computed on :attr:`~transformed_mapped_model_pos` and
        :attr:`mapped_target_pos`

        :type: :class:`float`
        """
        if self._gdtha is None:
            n = \
            self.mapped_target_pos.GetGDTHA(self.transformed_mapped_model_pos,
                                            norm=False)
            n_full = 4*len(self.mapped_target_pos) + 4*self.n_target_not_mapped
            if n_full > 0:
                self._gdtha = float(n) / n_full
            else:
                self._gdtha = 0.0
        return self._gdtha

    @property
    def rmsd(self):
        """ RMSD

        Computed on :attr:`~transformed_mapped_model_pos` and
        :attr:`mapped_target_pos`

        :type: :class:`float`
        """
        if self._rmsd is None:
            self._rmsd = \
            self.mapped_target_pos.GetRMSD(self.transformed_mapped_model_pos)
        return self._rmsd

    @property
    def cad_score(self):
        """ The global CAD atom-atom (AA) score

        Computed based on :attr:`~model`. In case of oligomers, :attr:`~mapping`
        is used.

        :type: :class:`float`
        """
        if self._cad_score is None:
            self._compute_cad_score()
        return self._cad_score

    @property
    def local_cad_score(self):
        """ The per-residue CAD atom-atom (AA) scores

        Computed based on :attr:`~model`. In case of oligomers, :attr:`~mapping`
        is used.

        :type: :class:`dict`
        """
        if self._local_cad_score is None:
            self._compute_cad_score()
        return self._local_cad_score

    @property
    def patch_qs(self):
        """ Patch QS-scores for each residue in :attr:`model_interface_residues`

        Representative patches for each residue r in chain c are computed as
        follows:
    
        * mdl_patch_one: All residues in c with CB (CA for GLY) positions within
          8A of r and within 12A of residues from any other chain.
        * mdl_patch_two: Closest residue x to r in any other chain gets
          identified. Patch is then constructed by selecting all residues from
          any other chain within 8A of x and within 12A from any residue in c.
        * trg_patch_one: Chain name and residue number based mapping from
          mdl_patch_one
        * trg_patch_two: Chain name and residue number based mapping from
          mdl_patch_two

        Results are stored in the same manner as
        :attr:`model_interface_residues`, with corresponding scores instead of
        residue numbers. Scores for residues which are not
        :class:`mol.ChemType.AMINOACIDS` are set to None. Additionally,
        interface patches are derived from :attr:`model`. If they contain
        residues which are not covered by :attr:`target`, the score is set to
        None too.

        :type: :class:`dict` with chain names as key and and :class:`list`
                with scores of the respective interface residues.
        """
        if self._patch_qs is None:
            self._compute_patchqs_scores()
        return self._patch_qs

    @property
    def patch_dockq(self):
        """ Same as :attr:`patch_qs` but for DockQ scores
        """
        if self._patch_dockq is None:
            self._compute_patchdockq_scores()
        return self._patch_dockq

    @property
    def tm_score(self):
        """ TM-score computed with USalign

        USalign executable can be specified with usalign_exec kwarg at Scorer
        construction, an OpenStructure internal copy of the USalign code is
        used otherwise.

        :type: :class:`float`
        """
        if self._tm_score is None:
            self._compute_tmscore()
        return self._tm_score

    @property
    def usalign_mapping(self):
        """ Mapping computed with USalign

        Dictionary with target chain names as key and model chain names as
        values. No guarantee that all chains are mapped. USalign executable
        can be specified with usalign_exec kwarg at Scorer construction, an
        OpenStructure internal copy of the USalign code is used otherwise.

        :type: :class:`dict`
        """
        if self._usalign_mapping is None:
            self._compute_tmscore()
        return self._usalign_mapping

    def _aln_helper(self, target, model):
        # perform required alignments - cannot take the alignments from the
        # mapping results as we potentially remove stuff there as compared
        # to self.model and self.target
        trg_seqs = dict()
        for ch in target.chains:
            cname = ch.GetName()
            s = ''.join([r.one_letter_code for r in ch.residues])
            s = seq.CreateSequence(ch.GetName(), s)
            s.AttachView(target.Select(f"cname='{cname}'"))
            trg_seqs[ch.GetName()] = s
        mdl_seqs = dict()
        for ch in model.chains:
            cname = ch.GetName()
            s = ''.join([r.one_letter_code for r in ch.residues])
            s = seq.CreateSequence(cname, s)
            s.AttachView(model.Select(f"cname='{cname}'"))
            mdl_seqs[ch.GetName()] = s

        alns = list()
        trg_pep_chains = [s.GetName() for s in self.chain_mapper.polypep_seqs]
        trg_nuc_chains = [s.GetName() for s in self.chain_mapper.polynuc_seqs]
        trg_pep_chains = set(trg_pep_chains)
        trg_nuc_chains = set(trg_nuc_chains)
        for trg_ch, mdl_ch in self.mapping.GetFlatMapping().items():
            if mdl_ch in mdl_seqs and trg_ch in trg_seqs:
                if trg_ch in trg_pep_chains:
                    stype = mol.ChemType.AMINOACIDS
                elif trg_ch in trg_nuc_chains:
                    stype = mol.ChemType.NUCLEOTIDES
                else:
                    raise RuntimeError("Chain name inconsistency... ask "
                                       "Gabriel")
                alns.append(self.chain_mapper.Align(trg_seqs[trg_ch],
                                                    mdl_seqs[mdl_ch],
                                                    stype))
                alns[-1].AttachView(0, trg_seqs[trg_ch].GetAttachedView())
                alns[-1].AttachView(1, mdl_seqs[mdl_ch].GetAttachedView())
        return alns

    def _compute_pepnuc_aln(self):
        query = "peptide=true or nucleotide=true"
        pep_nuc_target = self.target_orig.Select(query)
        pep_nuc_model = self.model_orig.Select(query)
        self._pepnuc_aln = self._aln_helper(pep_nuc_target, pep_nuc_model)

    def _compute_aln(self):
        self._aln = self._aln_helper(self.target, self.model)

    def _compute_stereochecked_aln(self):
        self._stereochecked_aln = self._aln_helper(self.stereochecked_target,
                                                   self.stereochecked_model)

    def _compute_lddt(self):
        # lDDT requires a flat mapping with mdl_ch as key and trg_ch as value
        flat_mapping = self.mapping.GetFlatMapping(mdl_as_key=True)

        # make alignments accessible by mdl seq name
        stereochecked_alns = dict()
        for aln in self.stereochecked_aln:
            mdl_seq = aln.GetSequence(1)
            stereochecked_alns[mdl_seq.name] = aln
        alns = dict()
        for aln in self.aln:
            mdl_seq = aln.GetSequence(1)
            alns[mdl_seq.name] = aln

        # score variables to be set
        lddt_score = None
        local_lddt = None

        if self.lddt_no_stereochecks:
            lddt_chain_mapping = dict()
            for mdl_ch, trg_ch in flat_mapping.items():
                if mdl_ch in alns:
                    lddt_chain_mapping[mdl_ch] = trg_ch
            lddt_score = self.lddt_scorer.lDDT(self.model,
                                               chain_mapping = lddt_chain_mapping,
                                               residue_mapping = alns,
                                               check_resnames=False,
                                               local_lddt_prop="lddt")[0]
            local_lddt = dict()
            for r in self.model.residues:
                cname = r.GetChain().GetName()
                if cname not in local_lddt:
                    local_lddt[cname] = dict()
                if r.HasProp("lddt"):
                    score = round(r.GetFloatProp("lddt"), 3)
                    local_lddt[cname][r.GetNumber()] = score
                else:
                    # not covered by trg or skipped in chain mapping procedure
                    # the latter happens if its part of a super short chain
                    local_lddt[cname][r.GetNumber()] = None

        else:
            lddt_chain_mapping = dict()
            for mdl_ch, trg_ch in flat_mapping.items():
                if mdl_ch in stereochecked_alns:
                    lddt_chain_mapping[mdl_ch] = trg_ch
            lddt_score = self.lddt_scorer.lDDT(self.stereochecked_model,
                                               chain_mapping = lddt_chain_mapping,
                                               residue_mapping = stereochecked_alns,
                                               check_resnames=False,
                                               local_lddt_prop="lddt")[0]
            local_lddt = dict()
            for r in self.model.residues:
                cname = r.GetChain().GetName()
                if cname not in local_lddt:
                    local_lddt[cname] = dict()
                if r.HasProp("lddt"):
                    score = round(r.GetFloatProp("lddt"), 3)
                    local_lddt[cname][r.GetNumber()] = score
                else:
                    # rsc => residue stereo checked...
                    mdl_res = self.stereochecked_model.FindResidue(cname, r.GetNumber())
                    if mdl_res.IsValid():
                        # not covered by trg or skipped in chain mapping procedure
                        # the latter happens if its part of a super short chain
                        local_lddt[cname][r.GetNumber()] = None
                    else:
                        # opt 1: removed by stereochecks => assign 0.0
                        # opt 2: removed by stereochecks AND not covered by ref
                        #        => assign None
                        # fetch trg residue from non-stereochecked aln
                        trg_r = None
                        if cname in flat_mapping:
                            for col in alns[cname]:
                                if col[0] != '-' and col[1] != '-':
                                    if col.GetResidue(1).number == r.number:
                                        trg_r = col.GetResidue(0)
                                        break
                        if trg_r is None:
                            local_lddt[cname][r.GetNumber()] = None
                        else:
                            local_lddt[cname][r.GetNumber()] = 0.0

        self._lddt = lddt_score
        self._local_lddt = local_lddt

    def _compute_bb_lddt(self):

        # make alignments accessible by mdl seq name
        alns = dict()
        for aln in self.aln:
            mdl_seq = aln.GetSequence(1)
            alns[mdl_seq.name] = aln

        # lDDT requires a flat mapping with mdl_ch as key and trg_ch as value
        flat_mapping = self.mapping.GetFlatMapping(mdl_as_key=True)
        lddt_chain_mapping = dict()
        for mdl_ch, trg_ch in flat_mapping.items():
            if mdl_ch in alns:
                lddt_chain_mapping[mdl_ch] = trg_ch

        lddt_score = self.bb_lddt_scorer.lDDT(self.model,
                                              chain_mapping = lddt_chain_mapping,
                                              residue_mapping = alns,
                                              check_resnames=False,
                                              local_lddt_prop="bb_lddt")[0]
        local_lddt = dict()
        for r in self.model.residues:
            cname = r.GetChain().GetName()
            if cname not in local_lddt:
                local_lddt[cname] = dict()
            if r.HasProp("bb_lddt"):
                score = round(r.GetFloatProp("bb_lddt"), 3)
                local_lddt[cname][r.GetNumber()] = score
            else:
                # not covered by trg or skipped in chain mapping procedure
                # the latter happens if its part of a super short chain
                local_lddt[cname][r.GetNumber()] = None

        self._bb_lddt = lddt_score
        self._bb_local_lddt = local_lddt

    def _compute_qs(self):
        qs_score_result = self.qs_scorer.Score(self.mapping.mapping)
        self._qs_global = qs_score_result.QS_global
        self._qs_best = qs_score_result.QS_best

    def _compute_per_interface_qs_scores(self):
        self._per_interface_qs_global = list()
        self._per_interface_qs_best = list()

        for interface in self.qs_interfaces:
            trg_ch1 = interface[0]
            trg_ch2 = interface[1]
            mdl_ch1 = interface[2]
            mdl_ch2 = interface[3]
            qs_res = self.qs_scorer.ScoreInterface(trg_ch1, trg_ch2,
                                                   mdl_ch1, mdl_ch2)
            self._per_interface_qs_best.append(qs_res.QS_best)
            self._per_interface_qs_global.append(qs_res.QS_global)

    def _compute_ics_scores(self):
        contact_scorer_res = self.contact_scorer.ScoreICS(self.mapping.mapping)
        self._ics_precision = contact_scorer_res.precision
        self._ics_recall = contact_scorer_res.recall
        self._ics = contact_scorer_res.ics
        self._per_interface_ics_precision = list()
        self._per_interface_ics_recall = list()
        self._per_interface_ics = list()
        flat_mapping = self.mapping.GetFlatMapping()
        for trg_int in self.contact_target_interfaces:
            trg_ch1 = trg_int[0]
            trg_ch2 = trg_int[1]
            if trg_ch1 in flat_mapping and trg_ch2 in flat_mapping:
                mdl_ch1 = flat_mapping[trg_ch1]
                mdl_ch2 = flat_mapping[trg_ch2]
                res = self.contact_scorer.ScoreICSInterface(trg_ch1, trg_ch2,
                                                            mdl_ch1, mdl_ch2)
                self._per_interface_ics_precision.append(res.precision)
                self._per_interface_ics_recall.append(res.recall)
                self._per_interface_ics.append(res.ics)
            else:
                self._per_interface_ics_precision.append(None)
                self._per_interface_ics_recall.append(None)
                self._per_interface_ics.append(None)

    def _compute_ips_scores(self):
        contact_scorer_res = self.contact_scorer.ScoreIPS(self.mapping.mapping)
        self._ips_precision = contact_scorer_res.precision
        self._ips_recall = contact_scorer_res.recall
        self._ips = contact_scorer_res.ips

        self._per_interface_ips_precision = list()
        self._per_interface_ips_recall = list()
        self._per_interface_ips = list()
        flat_mapping = self.mapping.GetFlatMapping()
        for trg_int in self.contact_target_interfaces:
            trg_ch1 = trg_int[0]
            trg_ch2 = trg_int[1]
            if trg_ch1 in flat_mapping and trg_ch2 in flat_mapping:
                mdl_ch1 = flat_mapping[trg_ch1]
                mdl_ch2 = flat_mapping[trg_ch2]
                res = self.contact_scorer.ScoreIPSInterface(trg_ch1, trg_ch2,
                                                            mdl_ch1, mdl_ch2)
                self._per_interface_ips_precision.append(res.precision)
                self._per_interface_ips_recall.append(res.recall)
                self._per_interface_ips.append(res.ips)
            else:
                self._per_interface_ips_precision.append(None)
                self._per_interface_ips_recall.append(None)
                self._per_interface_ips.append(None)

    def _compute_dockq_scores(self):
        # lists with values in contact_target_interfaces
        self._dockq_scores = list()
        self._fnat = list()
        self._fnonnat = list()
        self._irmsd = list()
        self._lrmsd = list()
        self._nnat = list()
        self._nmdl = list()

        dockq_alns = dict()
        for aln in self.aln:
            trg_s = aln.GetSequence(0)
            mdl_s = aln.GetSequence(1)
            dockq_alns[(trg_s.GetName(), mdl_s.GetName())] = aln

        for interface in self.dockq_interfaces:
            trg_ch1 = interface[0]
            trg_ch2 = interface[1]
            mdl_ch1 = interface[2]
            mdl_ch2 = interface[3]
            aln1 = dockq_alns[(trg_ch1, mdl_ch1)]
            aln2 = dockq_alns[(trg_ch2, mdl_ch2)]
            res = dockq.DockQ(self.model, self.target, mdl_ch1, mdl_ch2,
                              trg_ch1, trg_ch2, ch1_aln=aln1,
                              ch2_aln=aln2)
            self._fnat.append(res["fnat"])
            self._fnonnat.append(res["fnonnat"])
            self._irmsd.append(res["irmsd"])
            self._lrmsd.append(res["lrmsd"])
            self._dockq_scores.append(res["DockQ"])
            self._nnat.append(res["nnat"])
            self._nmdl.append(res["nmdl"])

        # keep track of native counts in target interfaces which are
        # not covered in model in order to compute
        # dockq_ave_full/dockq_wave_full in the end
        not_covered_counts = list()
        proc_trg_interfaces = set([(x[0], x[1]) for x in self.dockq_interfaces])
        for interface in self.dockq_target_interfaces:
            if interface not in proc_trg_interfaces:
                # let's run DockQ with trg as trg/mdl in order to get the native
                # contacts out - no need to pass alns as the residue numbers
                # match for sure
                trg_ch1 = interface[0]
                trg_ch2 = interface[1]
                res = dockq.DockQ(self.target, self.target,
                                  trg_ch1, trg_ch2, trg_ch1, trg_ch2)
                not_covered_counts.apend(res["nnat"])
  
        # there are 4 types of combined scores
        # - simple average
        # - average weighted by native_contacts
        # - the two above including nonmapped_contact_interfaces => set DockQ to 0.0
        scores = np.array([self._dockq_scores])
        weights = np.array([self._nnat])
        if len(scores) > 0:
            self._dockq_ave = np.mean(scores)
        else:
            self._dockq_ave = 0.0
        self._dockq_wave = np.sum(np.multiply(weights/np.sum(weights), scores))
        scores = np.append(scores, [0.0]*len(not_covered_counts))
        weights = np.append(weights, not_covered_counts)
        if len(scores) > 0:
            self._dockq_ave_full = np.mean(scores)
        else:
            self._dockq_ave_full = 0.0
        self._dockq_wave_full = np.sum(np.multiply(weights/np.sum(weights),
                                                   scores))

    def _extract_mapped_pos(self):
        self._mapped_target_pos = geom.Vec3List()
        self._mapped_model_pos = geom.Vec3List()
        self._n_target_not_mapped = 0
        processed_trg_chains = set()
        for trg_ch, mdl_ch in self.mapping.GetFlatMapping().items():
            processed_trg_chains.add(trg_ch)
            aln = self.mapping.alns[(trg_ch, mdl_ch)]
            for col in aln:
                if col[0] != '-' and col[1] != '-':
                    trg_res = col.GetResidue(0)
                    mdl_res = col.GetResidue(1)
                    trg_at = trg_res.FindAtom("CA")
                    mdl_at = mdl_res.FindAtom("CA")
                    if not trg_at.IsValid():
                        trg_at = trg_res.FindAtom("C3'")
                    if not mdl_at.IsValid():
                        mdl_at = mdl_res.FindAtom("C3'")
                    self._mapped_target_pos.append(trg_at.GetPos())
                    self._mapped_model_pos.append(mdl_at.GetPos())
                elif col[0] != '-':
                    self._n_target_not_mapped += 1
        # count number of trg residues from non-mapped chains
        for ch in self.mapping.target.chains:
            if ch.GetName() not in processed_trg_chains:
                self._n_target_not_mapped += len(ch.residues)

    def _compute_cad_score(self):
        if not self.resnum_alignments:
            raise RuntimeError("CAD score computations rely on residue numbers "
                               "that are consistent between target and model "
                               "chains, i.e. only work if resnum_alignments "
                               "is True at Scorer construction.")
        try:
            cad_score_exec = \
            settings.Locate("voronota-cadscore",
                            explicit_file_name=self.cad_score_exec)
        except Exception as e:
            raise RuntimeError("voronota-cadscore must be in PATH for CAD "
                               "score scoring") from e
        cad_bin_dir = os.path.dirname(cad_score_exec)
        m = self.mapping.GetFlatMapping(mdl_as_key=True)
        cad_result = cadscore.CADScore(self.model, self.target,
                                       mode = "voronota",
                                       label="localcad",
                                       old_regime=False,
                                       cad_bin_path=cad_bin_dir,
                                       chain_mapping=m)

        local_cad = dict()
        for r in self.model.residues:
            cname = r.GetChain().GetName()
            if cname not in local_cad:
                local_cad[cname] = dict()
            if r.HasProp("localcad"):
                score = round(r.GetFloatProp("localcad"), 3)
                local_cad[cname][r.GetNumber()] = score
            else:
                local_cad[cname][r.GetNumber()] = None

        self._cad_score = cad_result.globalAA
        self._local_cad_score = local_cad

    def _get_repr_view(self, ent):
        """ Returns view with representative peptide atoms => CB, CA for GLY
    
        Ensures that each residue has exactly one atom with assertions
    
        :param ent: Entity for which you want the representative view
        :param ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
        :returns: :class:`ost.mol.EntityView` derived from ent
        """
        repr_ent = ent.Select("(aname=\"CB\" or (rname=\"GLY\" and aname=\"CA\"))")
        for r in repr_ent.residues:
            assert(len(r.atoms) == 1)
        return repr_ent

    def _get_interface_residues(self, ent):
        """ Get interface residues
    
        Thats all residues having a contact with at least one residue from another
        chain (CB-CB distance <= 8A, CA in case of Glycine)
    
        :param ent: Model for which to extract interface residues
        :type ent: :class:`ost.mol.EntityView`
        :returns: :class:`dict` with chain names as key and and :class:`list`
                  with residue numbers of the respective interface residues.
        """
        # select for representative positions => CB, CA for GLY 
        repr_ent = self._get_repr_view(ent)
        result = {ch.GetName(): list() for ch in ent.chains}
        for ch in ent.chains:
            cname = ch.GetName()
            sel = repr_ent.Select(f"(cname='{cname}' and 8 <> [cname!='{cname}'])")
            result[cname] = [r.GetNumber() for r in sel.residues]
        return result

    def _do_stereochecks(self):
        """ Perform stereochemistry checks on model and target
        """
        data = stereochemistry.GetDefaultStereoData()
        l_data = stereochemistry.GetDefaultStereoLinkData()

        a, b, c, d = stereochemistry.StereoCheck(self.model, stereo_data = data,
                                                 stereo_link_data = l_data)
        self._stereochecked_model = a
        self._model_clashes = b
        self._model_bad_bonds = c
        self._model_bad_angles = d

        a, b, c, d = stereochemistry.StereoCheck(self.target, stereo_data = data,
                                                 stereo_link_data = l_data)
        self._stereochecked_target = a
        self._target_clashes = b
        self._target_bad_bonds = c
        self._target_bad_angles = d


    def _get_interface_patches(self, mdl_ch, mdl_rnum):
        """ Select interface patches representative for specified residue
    
        The patches for specified residue r in chain c are selected as follows:
    
        * mdl_patch_one: All residues in c with CB (CA for GLY) positions within 8A
          of r and within 12A of residues from any other chain.
        * mdl_patch_two: Closest residue x to r in any other chain gets identified.
          Patch is then constructed by selecting all residues from any other chain
          within 8A of x and within 12A from any residue in c.
        * trg_patch_one: Chain name and residue number based mapping from
          mdl_patch_one
        * trg_patch_two: Chain name and residue number based mapping from
          mdl_patch_two
    
        :param mdl_ch: Name of chain in *self.model* of residue of interest
        :type mdl_ch: :class:`str`
        :param mdl_rnum: Residue number of residue of interest
        :type mdl_rnum: :class:`ost.mol.ResNum`
        :returns: Tuple with 5 elements: 1) :class:`bool` flag whether all residues
                  in *mdl* patches are covered in *trg* 2) mtl_patch_one
                  3) mdl_patch_two 4) trg_patch_one 5) trg_patch_two
        """
        # select for representative positions => CB, CA for GLY 
        repr_mdl = self._get_repr_view(self.model.Select("peptide=true"))
    
        # get position for specified residue
        r = self.model.FindResidue(mdl_ch, mdl_rnum)
        if not r.IsValid():
            raise RuntimeError(f"Cannot find residue {mdl_rnum} in chain {mdl_ch}")
        if r.GetName() == "GLY":
            at = r.FindAtom("CA")
        else:
            at = r.FindAtom("CB")
        if not at.IsValid():
            raise RuntimeError("Cannot find interface views for res without CB/CA")
        r_pos = at.GetPos()
    
        # mdl_patch_one contains residues from the same chain as r
        # => all residues within 8A of r and within 12A of any other chain
    
        # q1 selects for everything in same chain and within 8A of r_pos
        q1 = f"(cname='{mdl_ch}' and 8 <> {{{r_pos[0]},{r_pos[1]},{r_pos[2]}}})"
        # q2 selects for everything within 12A of any other chain
        q2 = f"(12 <> [cname!={mdl_ch}])"
        mdl_patch_one = self.model.CreateEmptyView()
        sel = repr_mdl.Select(" and ".join([q1, q2]))
        for r in sel.residues:
            mdl_r = self.model.FindResidue(r.GetChain().GetName(), r.GetNumber())
            mdl_patch_one.AddResidue(mdl_r, mol.ViewAddFlag.INCLUDE_ALL)
    
        # mdl_patch_two contains residues from all other chains. In detail:
        # the closest residue to r is identified in any other chain, and the
        # patch is filled with residues that are within 8A of that residue and
        # within 12A of chain from r
        sel = repr_mdl.Select(f"(cname!='{mdl_ch}')")
        close_stuff = sel.FindWithin(r_pos, 8)
        min_pos = None
        min_dist = 42.0
        for close_at in close_stuff:
            dist = geom.Distance(r_pos, close_at.GetPos())
            if dist < min_dist:
                min_pos = close_at.GetPos()
                min_dist = dist
    
        # q1 selects for everything not in mdl_ch but within 8A of min_pos
        q1 = f"(cname!='{mdl_ch}' and 8 <> {{{min_pos[0]},{min_pos[1]},{min_pos[2]}}})"
        # q2 selects for everything within 12A of mdl_ch
        q2 = f"(12 <> [cname='{mdl_ch}'])"
        mdl_patch_two = self.model.CreateEmptyView()
        sel = repr_mdl.Select(" and ".join([q1, q2]))
        for r in sel.residues:
            mdl_r = self.model.FindResidue(r.GetChain().GetName(), r.GetNumber())
            mdl_patch_two.AddResidue(mdl_r, mol.ViewAddFlag.INCLUDE_ALL)
    
        # transfer mdl residues to trg
        flat_mapping = self.mapping.GetFlatMapping(mdl_as_key=True)
        full_trg_coverage = True
        trg_patch_one = self.target.CreateEmptyView()
        for r in mdl_patch_one.residues:
            trg_r = None
            mdl_cname = r.GetChain().GetName()
            if mdl_cname in flat_mapping:
                aln = self.mapping.alns[(flat_mapping[mdl_cname], mdl_cname)]
                for col in aln:
                    if col[0] != '-' and col[1] != '-':
                        if col.GetResidue(1).GetNumber() == r.GetNumber():
                            trg_r = col.GetResidue(0)
                            break
            if trg_r is not None:
                trg_patch_one.AddResidue(trg_r.handle,
                                         mol.ViewAddFlag.INCLUDE_ALL)
            else:
                full_trg_coverage = False
    
        trg_patch_two = self.target.CreateEmptyView()
        for r in mdl_patch_two.residues:
            trg_r = None
            mdl_cname = r.GetChain().GetName()
            if mdl_cname in flat_mapping:
                aln = self.mapping.alns[(flat_mapping[mdl_cname], mdl_cname)]
                for col in aln:
                    if col[0] != '-' and col[1] != '-':
                        if col.GetResidue(1).GetNumber() == r.GetNumber():
                            trg_r = col.GetResidue(0)
                            break
            if trg_r is not None:
                trg_patch_two.AddResidue(trg_r.handle,
                                         mol.ViewAddFlag.INCLUDE_ALL)
            else:
                full_trg_coverage = False
    
        return (full_trg_coverage, mdl_patch_one, mdl_patch_two,
                trg_patch_one, trg_patch_two)

    def _compute_patchqs_scores(self):
        self._patch_qs = dict()
        for cname, rnums in self.model_interface_residues.items():
            scores = list()
            for rnum in rnums:
                score = None
                r = self.model.FindResidue(cname, rnum)
                if r.IsValid() and r.GetChemType() == mol.ChemType.AMINOACIDS:
                    full_trg_coverage, mdl_patch_one, mdl_patch_two, \
                    trg_patch_one, trg_patch_two = \
                    self._get_interface_patches(cname, rnum)
                    if full_trg_coverage:
                        score = self._patchqs(mdl_patch_one, mdl_patch_two,
                                              trg_patch_one, trg_patch_two)
                scores.append(score)
            self._patch_qs[cname] = scores

    def _compute_patchdockq_scores(self):
        self._patch_dockq = dict()
        for cname, rnums in self.model_interface_residues.items():
            scores = list()
            for rnum in rnums:
                score = None
                r = self.model.FindResidue(cname, rnum)
                if r.IsValid() and r.GetChemType() == mol.ChemType.AMINOACIDS:
                    full_trg_coverage, mdl_patch_one, mdl_patch_two, \
                    trg_patch_one, trg_patch_two = \
                    self._get_interface_patches(cname, rnum)
                    if full_trg_coverage:
                        score = self._patchdockq(mdl_patch_one, mdl_patch_two,
                                                 trg_patch_one, trg_patch_two)
                scores.append(score)
            self._patch_dockq[cname] = scores

    def _patchqs(self, mdl_patch_one, mdl_patch_two, trg_patch_one, trg_patch_two):
        """ Score interface residue patches with QS-score
    
        In detail: Construct two entities with two chains each. First chain
        consists of residues from <x>_patch_one and second chain consists of
        <x>_patch_two. The returned score is the QS-score between the two
        entities
    
        :param mdl_patch_one: Interface patch representing scored residue
        :type mdl_patch_one: :class:`ost.mol.EntityView`
        :param mdl_patch_two: Interface patch representing scored residue
        :type mdl_patch_two: :class:`ost.mol.EntityView`
        :param trg_patch_one: Interface patch representing scored residue
        :type trg_patch_one: :class:`ost.mol.EntityView`
        :param trg_patch_two: Interface patch representing scored residue
        :type trg_patch_two: :class:`ost.mol.EntityView`
        :returns: PatchQS score
        """
        qs_ent_mdl = self._qs_ent_from_patches(mdl_patch_one, mdl_patch_two)
        qs_ent_trg = self._qs_ent_from_patches(trg_patch_one, trg_patch_two)
    
        alnA = seq.CreateAlignment()
        s = ''.join([r.one_letter_code for r in mdl_patch_one.residues])
        alnA.AddSequence(seq.CreateSequence("A", s))
        s = ''.join([r.one_letter_code for r in trg_patch_one.residues])
        alnA.AddSequence(seq.CreateSequence("A", s))
    
        alnB = seq.CreateAlignment()
        s = ''.join([r.one_letter_code for r in mdl_patch_two.residues])
        alnB.AddSequence(seq.CreateSequence("B", s))
        s = ''.join([r.one_letter_code for r in trg_patch_two.residues])
        alnB.AddSequence(seq.CreateSequence("B", s))
        alns = {("A", "A"): alnA, ("B", "B"): alnB}
    
        scorer = QSScorer(qs_ent_mdl, [["A"], ["B"]], qs_ent_trg, alns)
        score_result = scorer.Score([["A"], ["B"]])
    
        return score_result.QS_global
    
    def _patchdockq(self, mdl_patch_one, mdl_patch_two, trg_patch_one,
                    trg_patch_two):
        """ Score interface residue patches with DockQ
    
        In detail: Construct two entities with two chains each. First chain
        consists of residues from <x>_patch_one and second chain consists of
        <x>_patch_two. The returned score is the QS-score between the two
        entities
    
        :param mdl_patch_one: Interface patch representing scored residue
        :type mdl_patch_one: :class:`ost.mol.EntityView`
        :param mdl_patch_two: Interface patch representing scored residue
        :type mdl_patch_two: :class:`ost.mol.EntityView`
        :param trg_patch_one: Interface patch representing scored residue
        :type trg_patch_one: :class:`ost.mol.EntityView`
        :param trg_patch_two: Interface patch representing scored residue
        :type trg_patch_two: :class:`ost.mol.EntityView`
        :returns: DockQ score
        """
        m = self._qs_ent_from_patches(mdl_patch_one, mdl_patch_two)
        t = self._qs_ent_from_patches(trg_patch_one, trg_patch_two)
        dockq_result = dockq.DockQ(t, m, "A", "B", "A", "B")
        if dockq_result["nnat"] > 0:
            return dockq_result["DockQ"]
        return 0.0

    def _qs_ent_from_patches(self, patch_one, patch_two):
        """ Constructs Entity with two chains named "A" and "B""
    
        Blindly adds all residues from *patch_one* to chain A and residues from
        patch_two to chain B.
        """
        ent = mol.CreateEntity()
        ed = ent.EditXCS()
        added_ch = ed.InsertChain("A")
        for r in patch_one.residues:
            added_r = ed.AppendResidue(added_ch, r.GetName())
            added_r.SetChemClass(str(r.GetChemClass()))
            for a in r.atoms:
                ed.InsertAtom(added_r, a.handle)
        added_ch = ed.InsertChain("B")
        for r in patch_two.residues:
            added_r = ed.AppendResidue(added_ch, r.GetName())
            added_r.SetChemClass(str(r.GetChemClass()))
            for a in r.atoms:
                ed.InsertAtom(added_r, a.handle)
        return ent

    def _set_custom_mapping(self, mapping):
        """ sets self._mapping with a full blown MappingResult object

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
                    trg_view = chain_mapper.target.Select(f"cname='{ref_ch}'")
                    mdl_view = mdl.Select(f"cname='{mdl_ch}'")
                    aln.AttachView(0, trg_view)
                    aln.AttachView(1, mdl_view)
                    alns[(ref_ch, mdl_ch)] = aln

        self._mapping = chain_mapping.MappingResult(chain_mapper.target, mdl,
                                                    chain_mapper.chem_groups,
                                                    chem_mapping,
                                                    final_mapping, alns)

    def _compute_tmscore(self):
        res = None
        if self.usalign_exec is None:
            if self.oum:
                flat_mapping = self.mapping.GetFlatMapping()
                res = res = bindings.WrappedMMAlign(self.model, self.target,
                                                    mapping=flat_mapping)
            else:
                res = bindings.WrappedMMAlign(self.model, self.target)
        else:
            if self.oum:
                flat_mapping = self.mapping.GetFlatMapping()
                res = tmtools.USAlign(self.model, self.target,
                                      usalign = self.usalign_exec,
                                      custom_chain_mapping = flat_mapping)
            else:
                res = tmtools.USAlign(self.model, self.target,
                                      usalign = self.usalign_exec)

        self._tm_score = res.tm_score
        self._usalign_mapping = dict()
        for a,b in zip(res.ent1_mapped_chains, res.ent2_mapped_chains):
            self._usalign_mapping[b] = a
