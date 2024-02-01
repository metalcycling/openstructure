import numpy as np

from ost import mol
from ost import conop


class CustomCompound:
    """ Defines atoms for custom compounds

    lDDT requires the reference atoms of a compound which are typically
    extracted from a :class:`ost.conop.CompoundLib`. This lightweight
    container allows to handle arbitrary compounds which are not
    necessarily in the compound library.

    :param atom_names: Names of atoms of custom compound
    :type atom_names: :class:`list` of :class:`str`
    """
    def __init__(self, atom_names):
        self.atom_names = atom_names

    @staticmethod
    def FromResidue(res):
        """ Construct custom compound from residue

        :param res: Residue from which reference atom names are extracted,
                    hydrogen/deuterium atoms are filtered out
        :type res: :class:`ost.mol.ResidueView`/:class:`ost.mol.ResidueHandle`
        :returns: :class:`CustomCompound`
        """
        at_names = [a.name for a in res.atoms if a.element not in ["H", "D"]]
        if len(at_names) != len(set(at_names)):
            raise RuntimeError("Duplicate atoms detected in CustomCompound")
        compound = CustomCompound(at_names)
        return compound

class SymmetrySettings:
    """Container for symmetric compounds

    lDDT considers symmetries and selects the one resulting in the highest
    possible score.

    A symmetry is defined as a renaming operation on one or more atoms that
    leads to a chemically equivalent residue. Example would be OD1 and OD2 in
    ASP => renaming OD1 to OD2 and vice versa gives a chemically equivalent
    residue.

    Use :func:`AddSymmetricCompound` to define a symmetry which can then
    directly be accessed through the *symmetric_compounds* member.
    """
    def __init__(self):
        self.symmetric_compounds = dict()

    def AddSymmetricCompound(self, name, symmetric_atoms):
        """Adds symmetry for compound with *name*

        :param name: Name of compound with symmetry
        :type name: :class:`str`
        :param symmetric_atoms: Pairs of atom names that define renaming
                                operation, i.e. after applying all switches
                                defined in the tuples, the resulting residue
                                should be chemically equivalent. Atom names
                                must refer to the PDB component dictionary.
        :type symmetric_atoms: :class:`list` of :class:`tuple`
        """
        for pair in symmetric_atoms:
            if len(pair) != 2:
                raise RuntimeError("Expect pairs when defining symmetries")
        self.symmetric_compounds[name] = symmetric_atoms


def GetDefaultSymmetrySettings():
    """Constructs and returns :class:`SymmetrySettings` object for natural amino
    acids
    """
    symmetry_settings = SymmetrySettings()

    # ASP
    symmetry_settings.AddSymmetricCompound("ASP", [("OD1", "OD2")])

    # GLU
    symmetry_settings.AddSymmetricCompound("GLU", [("OE1", "OE2")])

    # LEU
    symmetry_settings.AddSymmetricCompound("LEU", [("CD1", "CD2")])

    # VAL
    symmetry_settings.AddSymmetricCompound("VAL", [("CG1", "CG2")])

    # ARG
    symmetry_settings.AddSymmetricCompound("ARG", [("NH1", "NH2")])

    # PHE
    symmetry_settings.AddSymmetricCompound(
        "PHE", [("CD1", "CD2"), ("CE1", "CE2")]
    )

    # TYR
    symmetry_settings.AddSymmetricCompound(
        "TYR", [("CD1", "CD2"), ("CE1", "CE2")]
    )

    return symmetry_settings


class lDDTScorer:
    """lDDT scorer object for a specific target

    Sets up everything to score models of that target. lDDT (local distance
    difference test) is defined as fraction of pairwise distances which exhibit
    a difference < threshold when considering target and model. In case of
    multiple thresholds, the average is returned. See

    V. Mariani, M. Biasini, A. Barbato, T. Schwede, lDDT : A local
    superposition-free score for comparing protein structures and models using
    distance difference tests, Bioinformatics, 2013

    :param target: The target
    :type target: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param compound_lib: Compound library from which a compound for each residue
                         is extracted based on its name. Uses
                         :func:`ost.conop.GetDefaultLib` if not given, raises
                         if this returns no valid compound library. Atoms
                         defined in the compound are searched in the residue and
                         build the reference for scoring. If the residue has
                         atoms with names ["A", "B", "C"] but the corresponding
                         compound only has ["A", "B"], "A" and "B" are
                         considered for scoring. If the residue has atoms
                         ["A", "B"] but the compound has ["A", "B", "C"], "C" is
                         considered missing and does not influence scoring, even
                         if present in the model.
    :param custom_compounds: Custom compounds defining reference atoms. If
                             given, *custom_compounds* take precedent over
                             *compound_lib*.
    :type custom_compounds: :class:`dict` with residue names (:class:`str`) as
                            key and :class:`CustomCompound` as value.
    :type compound_lib: :class:`ost.conop.CompoundLib`
    :param inclusion_radius: All pairwise distances < *inclusion_radius* are
                             considered for scoring
    :type inclusion_radius: :class:`float`
    :param sequence_separation: Only pairwise distances between atoms of
                                residues which are further apart than this
                                threshold are considered. Residue distance is
                                based on resnum. The default (0) considers all
                                pairwise distances except intra-residue
                                distances.
    :type sequence_separation: :class:`int`
    :param symmetry_settings: Define residues exhibiting internal symmetry, uses
                              :func:`GetDefaultSymmetrySettings` if not given.
    :type symmetry_settings: :class:`SymmetrySettings`
    :param seqres_mapping: Mapping of model residues at the scoring stage
                           happens with residue numbers defining their location
                           in a reference sequence (SEQRES) using one based
                           indexing. If the residue numbers in *target* don't
                           correspond to that SEQRES, you can specify the
                           mapping manually. You can provide a dictionary to
                           specify a reference sequence (SEQRES) for one or more
                           chain(s). Key: chain name, value: alignment
                           (seq1: SEQRES, seq2: sequence of residues in chain).
                           Example: The residues in a chain with name "A" have
                           sequence "YEAH" and residue numbers [42,43,44,45].
                           You can provide an alignment with seq1 "``HELLYEAH``"
                           and seq2 "``----YEAH``". "Y" gets assigned residue
                           number 5, "E" gets assigned 6 and so on no matter
                           what the original residue numbers were. 
    :type seqres_mapping: :class:`dict` (key: :class:`str`, value:
                          :class:`ost.seq.AlignmentHandle`)
    :param bb_only: Only consider atoms with name "CA" in case of amino acids and
                    "C3'" for Nucleotides. this invalidates *compound_lib*.
                    Raises if any residue in *target* is not
                    `r.chem_class.IsPeptideLinking()` or
                    `r.chem_class.IsNucleotideLinking()`
    :type bb_only: :class:`bool`
    :raises: :class:`RuntimeError` if *target* contains compound which is not in
             *compound_lib*, :class:`RuntimeError` if *symmetry_settings*
             specifies symmetric atoms that are not present in the according
             compound in *compound_lib*, :class:`RuntimeError` if
             *seqres_mapping* is not provided and *target* contains residue
             numbers with insertion codes or the residue numbers for each chain
             are not monotonically increasing, :class:`RuntimeError` if
             *seqres_mapping* is provided but an alignment is invalid
             (seq1 contains gaps, mismatch in seq1/seq2, seq2 does not match
             residues in corresponding chains).
    """
    def __init__(
        self,
        target,
        compound_lib=None,
        custom_compounds=None,
        inclusion_radius=15,
        sequence_separation=0,
        symmetry_settings=None,
        seqres_mapping=dict(),
        bb_only=False
    ):

        self.target = target
        self.inclusion_radius = inclusion_radius
        self.sequence_separation = sequence_separation
        if compound_lib is None:
            compound_lib = conop.GetDefaultLib()
        if compound_lib is None:
            raise RuntimeError("No compound_lib given and conop.GetDefaultLib "
                               "returns no valid compound library")
        self.compound_lib = compound_lib
        self.custom_compounds = custom_compounds
        if symmetry_settings is None:
            self.symmetry_settings = GetDefaultSymmetrySettings()
        else:
            self.symmetry_settings = symmetry_settings

        # whether to only consider atoms with name "CA" (amino acids) or C3'
        # (nucleotides), invalidates *compound_lib*
        self.bb_only=bb_only

        # names of heavy atoms of each unique compound present in *target* as
        # extracted from *compound_lib*, e.g.
        # self.compound_anames["GLY"] = ["N", "CA", "C", "O"]
        self.compound_anames = dict()

        # stores symmetry information for those compounds as defined in
        # *symmetry_settings*
        self.compound_symmetric_atoms = dict()

        # list of len(target.chains) containing all chain names in *target*
        self.chain_names = list()

        # list of len(target.residues) containing all compound names in *target*
        self.compound_names = list()

        # list of len(target.residues) defining start pos in internal reference
        # positions for each residue
        self.res_start_indices = list()

        # list of len(target.residues) defining residue numbers in target
        self.res_resnums = list()

        # list of len(target.chains) defining start pos in internal reference
        # positions for each chain     
        self.chain_start_indices = list()

        # list of len(target.chains) defining start pos in self.compound_names
        # for each chain     
        self.chain_res_start_indices = list()

        # maps residues in *target* to indices in
        # self.compound_names/self.res_start_indices. A residue gets identified
        # by a tuple (first element: chain name, second element: residue number,
        # residue number is either the actual residue number in *target* or
        # given by *seqres_mapping*)
        self.res_mapper = dict()

        # number of atoms as specified in compounds. not all are necessarily
        # covered by structure
        self.n_atoms = None

        # stores an index for each AtomHandle in *target*
        # (atom hashcode => index)
        self.atom_indices = dict()

        # store indices of all atoms that have symmetry properties
        self.symmetric_atoms = set()

        # setup members defined above
        self._SetupEnv(self.compound_lib, self.custom_compounds,
                       self.symmetry_settings, seqres_mapping, self.bb_only)

        # distance related members are lazily computed as they're affected
        # by different flavours of lDDT (e.g. lDDT including inter-chain
        # contacts or not etc.)

        # stores for each atom the other atoms within inclusion_radius
        self._ref_indices = None
        # the corresponding distances
        self._ref_distances = None

        # The following lists will be sparsely populated. We keep for each
        # symmetry related atom the distances towards all atoms which are NOT
        # affected by symmetry. So we can evaluate two symmetric versions
        # against the fixed stuff later on and select the better scoring one.
        self._sym_ref_indices = None
        self._sym_ref_distances = None

        # total number of distances
        self._n_distances = None

        # exactly the same as above but without interchain contacts
        # => single-chain (sc)
        self._ref_indices_sc = None
        self._ref_distances_sc = None
        self._sym_ref_indices_sc = None
        self._sym_ref_distances_sc = None
        self._n_distances_sc = None

        # exactly the same as above but without intrachain contacts
        # => inter-chain (ic)
        self._ref_indices_ic = None
        self._ref_distances_ic = None
        self._sym_ref_indices_ic = None
        self._sym_ref_distances_ic = None
        self._n_distances_ic = None

        # input parameter checking
        self._ProcessSequenceSeparation()

    @property
    def ref_indices(self):
        if self._ref_indices is None:
            self._SetupDistances()
        return self._ref_indices

    @property
    def ref_distances(self):
        if self._ref_distances is None:
            self._SetupDistances()
        return self._ref_distances
    
    @property
    def sym_ref_indices(self):
        if self._sym_ref_indices is None:
            self._SetupDistances()
        return self._sym_ref_indices

    @property
    def sym_ref_distances(self):
        if self._sym_ref_distances is None:
            self._SetupDistances()
        return self._sym_ref_distances

    @property
    def n_distances(self):
        if self._n_distances is None:
            self._n_distances = sum([len(x) for x in self.ref_indices])
        return self._n_distances

    @property
    def ref_indices_sc(self):
        if self._ref_indices_sc is None:
            self._SetupDistancesSC()
        return self._ref_indices_sc

    @property
    def ref_distances_sc(self):
        if self._ref_distances_sc is None:
            self._SetupDistancesSC()
        return self._ref_distances_sc
    
    @property
    def sym_ref_indices_sc(self):
        if self._sym_ref_indices_sc is None:
            self._SetupDistancesSC()
        return self._sym_ref_indices_sc

    @property
    def sym_ref_distances_sc(self):
        if self._sym_ref_distances_sc is None:
            self._SetupDistancesSC()
        return self._sym_ref_distances_sc

    @property
    def n_distances_sc(self):
        if self._n_distances_sc is None:
            self._n_distances_sc = sum([len(x) for x in self.ref_indices_sc])
        return self._n_distances_sc

    @property
    def ref_indices_ic(self):
        if self._ref_indices_ic is None:
            self._SetupDistancesIC()
        return self._ref_indices_ic

    @property
    def ref_distances_ic(self):
        if self._ref_distances_ic is None:
            self._SetupDistancesIC()
        return self._ref_distances_ic
    
    @property
    def sym_ref_indices_ic(self):
        if self._sym_ref_indices_ic is None:
            self._SetupDistancesIC()
        return self._sym_ref_indices_ic

    @property
    def sym_ref_distances_ic(self):
        if self._sym_ref_distances_ic is None:
            self._SetupDistancesIC()
        return self._sym_ref_distances_ic

    @property
    def n_distances_ic(self):
        if self._n_distances_ic is None:
            self._n_distances_ic = sum([len(x) for x in self.ref_indices_ic])
        return self._n_distances_ic

    def lDDT(self, model, thresholds = [0.5, 1.0, 2.0, 4.0],
             local_lddt_prop=None, local_contact_prop=None,
             chain_mapping=None, no_interchain=False,
             no_intrachain=False, penalize_extra_chains=False,
             residue_mapping=None, return_dist_test=False,
             check_resnames=True):
        """Computes lDDT of *model* - globally and per-residue

        :param model: Model to be scored - models are preferably scored upon
                      performing stereo-chemistry checks in order to punish for
                      non-sensical irregularities. This must be done separately
                      as a pre-processing step. Target contacts that are not
                      covered by *model* are considered not conserved, thus
                      decreasing lDDT score. This also includes missing model
                      chains or model chains for which no mapping is provided in
                      *chain_mapping*.
        :type model: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
        :param thresholds: Thresholds of distance differences to be considered
                           as correct - see docs in constructor for more info.
                           default: [0.5, 1.0, 2.0, 4.0]
        :type thresholds: :class:`list` of :class:`floats`
        :param local_lddt_prop: If set, per-residue scores will be assigned as
                                generic float property of that name
        :type local_lddt_prop: :class:`str`
        :param local_contact_prop: If set, number of expected contacts as well
                                   as number of conserved contacts will be
                                   assigned as generic int property.
                                   Excected contacts will be set as
                                   <local_contact_prop>_exp, conserved contacts
                                   as <local_contact_prop>_cons. Values
                                   are summed over all thresholds.
        :type local_contact_prop: :class:`str`
        :param chain_mapping: Mapping of model chains (key) onto target chains
                              (value). This is required if target or model have
                              more than one chain.
        :type chain_mapping: :class:`dict` with :class:`str` as keys/values
        :param no_interchain: Whether to exclude interchain contacts
        :type no_interchain: :class:`bool`
        :param no_intrachain: Whether to exclude intrachain contacts (i.e. only
                              consider interface related contacts)
        :type no_intrachain: :class:`bool`
        :param penalize_extra_chains: Whether to include a fixed penalty for
                                      additional chains in the model that are
                                      not mapped to the target. ONLY AFFECTS
                                      RETURNED GLOBAL SCORE. In detail: adds the
                                      number of intra-chain contacts of each
                                      extra chain to the expected contacts, thus
                                      adding a penalty.
        :param penalize_extra_chains: :class:`bool`
        :param residue_mapping: By default, residue mapping is based on residue
                                numbers. That means, a model chain and the
                                respective target chain map to the same
                                underlying reference sequence (SEQRES).
                                Alternatively, you can specify one or
                                several alignment(s) between model and target
                                chains by providing a dictionary. key: Name
                                of chain in model (respective target chain is
                                extracted from *chain_mapping*),
                                value: Alignment with first sequence
                                corresponding to target chain and second
                                sequence to model chain. There is NO reference
                                sequence involved, so the two sequences MUST
                                exactly match the actual residues observed in
                                the respective target/model chains (ATOMSEQ).
        :type residue_mapping: :class:`dict` with key: :class:`str`,
                               value: :class:`ost.seq.AlignmentHandle`
        :param return_dist_test: Whether to additionally return the underlying
                                 per-residue data for the distance difference
                                 test. Adds five objects to the return tuple.
                                 First: Number of total contacts summed over all
                                 thresholds
                                 Second: Number of conserved contacts summed
                                 over all thresholds
                                 Third: list with length of scored residues.
                                 Contains indices referring to model.residues.
                                 Fourth: numpy array of size
                                 len(scored_residues) containing the number of
                                 total contacts,
                                 Fifth: numpy matrix of shape 
                                 (len(scored_residues), len(thresholds))
                                 specifying how many for each threshold are
                                 conserved.
        :param check_resnames: On by default. Enforces residue name matches
                               between mapped model and target residues.
        :type check_resnames: :class:`bool`
        :returns: global and per-residue lDDT scores as a tuple -
                  first element is global lDDT score and second element
                  a list of per-residue scores with length len(*model*.residues)
                  None is assigned to residues that are not covered by target
        """
        if chain_mapping is None:
            if len(self.chain_names) > 1 or len(model.chains) > 1:
                raise NotImplementedError("Must provide chain mapping if "
                                          "target or model have > 1 chains.")
            chain_mapping = {model.chains[0].GetName(): self.chain_names[0]}
        else:
            # check whether chains specified in mapping exist
            for model_chain, target_chain in chain_mapping.items():
                if target_chain not in self.chain_names:
                    raise RuntimeError(f"Target chain specified in "
                                       f"chain_mapping ({target_chain}) does "
                                       f"not exist. Target has chains: "
                                       f"{self.chain_names}")
                ch = model.FindChain(model_chain)
                if not ch.IsValid():
                    raise RuntimeError(f"Model chain specified in "
                                       f"chain_mapping ({model_chain}) does "
                                       f"not exist. Model has chains: "
                                       f"{[c.GetName() for c in model.chains]}")

        # initialize positions with values far in nirvana. If a position is not
        # set, it should be far away from any position in model.
        max_pos = model.bounds.GetMax()
        max_coordinate = abs(max(max_pos[0], max_pos[1], max_pos[2]))
        max_coordinate += 42 * max(thresholds)
        pos = np.ones((self.n_atoms, 3), dtype=np.float32) * max_coordinate

        # for each scored residue in model a list of indices describing the
        # atoms from the reference that should be there
        res_ref_atom_indices = list()

        # for each scored residue in model a list of indices of atoms that are
        # actually there
        res_atom_indices = list()

        # indices of the scored residues
        res_indices = list()

        # Will contain one element per symmetry group
        symmetries = list()

        current_model_res_idx = -1
        for ch in model.chains:
            model_ch_name = ch.GetName()
            if model_ch_name not in chain_mapping:
                current_model_res_idx += len(ch.residues)
                continue # additional model chain which is not mapped
            target_ch_name = chain_mapping[model_ch_name]

            rnums = self._GetChainRNums(ch, residue_mapping, model_ch_name,
                                        target_ch_name)

            for r, rnum in zip(ch.residues, rnums):
                current_model_res_idx += 1
                res_mapper_key = (target_ch_name, rnum)
                if res_mapper_key not in self.res_mapper:
                    continue
                r_idx = self.res_mapper[res_mapper_key]
                if check_resnames and r.name != self.compound_names[r_idx]:
                    raise RuntimeError(
                        f"Residue name mismatch for {r}, "
                        f" expect {self.compound_names[r_idx]}"
                    )
                res_start_idx = self.res_start_indices[r_idx]
                rname = self.compound_names[r_idx]
                anames = self.compound_anames[rname]
                atoms = [r.FindAtom(aname) for aname in anames]
                res_ref_atom_indices.append(
                    list(range(res_start_idx, res_start_idx + len(anames)))
                )
                res_atom_indices.append(list())
                res_indices.append(current_model_res_idx)
                for a_idx, a in enumerate(atoms):
                    if a.IsValid():
                        p = a.GetPos()
                        pos[res_start_idx + a_idx][0] = p[0]
                        pos[res_start_idx + a_idx][1] = p[1]
                        pos[res_start_idx + a_idx][2] = p[2]
                        res_atom_indices[-1].append(res_start_idx + a_idx)
                if rname in self.compound_symmetric_atoms:
                    sym_indices = list()
                    for sym_tuple in self.compound_symmetric_atoms[rname]:
                        a_one = atoms[sym_tuple[0]]
                        a_two = atoms[sym_tuple[1]]
                        if a_one.IsValid() and a_two.IsValid():
                            sym_indices.append(
                                (
                                    res_start_idx + sym_tuple[0],
                                    res_start_idx + sym_tuple[1],
                                )
                            )
                    if len(sym_indices) > 0:
                        symmetries.append(sym_indices)

        if no_interchain and no_intrachain:
            raise RuntimeError("on_interchain and no_intrachain flags are "
                               "mutually exclusive")

        if no_interchain:
            sym_ref_indices = self.sym_ref_indices_sc
            sym_ref_distances = self.sym_ref_distances_sc
            ref_indices = self.ref_indices_sc
            ref_distances = self.ref_distances_sc
            n_distances = self.n_distances_sc
        elif no_intrachain:
            sym_ref_indices = self.sym_ref_indices_ic
            sym_ref_distances = self.sym_ref_distances_ic
            ref_indices = self.ref_indices_ic
            ref_distances = self.ref_distances_ic
            n_distances = self.n_distances_ic
        else:
            sym_ref_indices = self.sym_ref_indices
            sym_ref_distances = self.sym_ref_distances
            ref_indices = self.ref_indices
            ref_distances = self.ref_distances
            n_distances = self.n_distances

        self._ResolveSymmetries(pos, thresholds, symmetries, sym_ref_indices,
                                sym_ref_distances)

        per_res_exp = np.asarray([self._GetNExp(res_ref_atom_indices[idx],
            ref_indices) for idx in range(len(res_indices))], dtype=np.int32)
        per_res_conserved = self._EvalResidues(pos, thresholds,
                                               res_atom_indices,
                                               ref_indices, ref_distances)

        n_thresh = len(thresholds)

        # do per-residue scores
        per_res_lDDT = [None] * len(model.residues)
        for idx in range(len(res_indices)):
            n_exp = n_thresh * per_res_exp[idx]
            if n_exp > 0:
                score = np.sum(per_res_conserved[idx,:]) / n_exp
                per_res_lDDT[res_indices[idx]] = score
            else:
                per_res_lDDT[res_indices[idx]] = 0.0


        # do full model score
        if penalize_extra_chains:
            n_distances += self._GetExtraModelChainPenalty(model, chain_mapping)

        lDDT_tot = int(n_thresh * n_distances)
        lDDT_cons = int(np.sum(per_res_conserved))
        lDDT = None
        if lDDT_tot > 0:
            lDDT = float(lDDT_cons) / lDDT_tot

        # set properties if necessary
        if local_lddt_prop:
            residues = model.residues
            for idx in res_indices:
                residues[idx].SetFloatProp(local_lddt_prop, per_res_lDDT[idx])

        if local_contact_prop:
            residues = model.residues
            exp_prop = local_contact_prop + "_exp"
            conserved_prop = local_contact_prop + "_cons"

            for i, r_idx in enumerate(res_indices):
                residues[r_idx].SetIntProp(exp_prop,
                                           n_thresh * int(per_res_exp[i]))
                residues[r_idx].SetIntProp(conserved_prop,
                                           int(np.sum(per_res_conserved[i,:])))

        if return_dist_test:
            return lDDT, per_res_lDDT, lDDT_tot, lDDT_cons, res_indices, \
            per_res_exp, per_res_conserved
        else:
            return lDDT, per_res_lDDT

    def GetNChainContacts(self, target_chain, no_interchain=False):
        """Returns number of contacts expected for a certain chain in *target*

        :param target_chain: Chain in *target* for which you want the number
                             of expected contacts
        :type target_chain: :class:`str`
        :param no_interchain: Whether to exclude interchain contacts
        :type no_interchain: :class:`bool`
        :raises: :class:`RuntimeError` if specified chain doesnt exist
        """
        if target_chain not in self.chain_names:
            raise RuntimeError(f"Specified chain name ({target_chain}) not in "
                               f"target")
        ch_idx = self.chain_names.index(target_chain)
        s = self.chain_start_indices[ch_idx]
        e = self.n_atoms
        if ch_idx + 1 < len(self.chain_names):
            e = self.chain_start_indices[ch_idx+1]
        if no_interchain:
            return self._GetNExp(list(range(s, e)), self.ref_indices_sc)
        else:
            return self._GetNExp(list(range(s, e)), self.ref_indices)


    def _GetExtraModelChainPenalty(self, model, chain_mapping):
        """Counts n distances in extra model chains to be added as penalty
        """
        penalty = 0
        for chain in model.chains:
            ch_name = chain.GetName()
            if ch_name not in chain_mapping:
                sm = self.symmetry_settings
                mdl_sel = model.Select(f"cname={mol.QueryQuoteName(ch_name)}")
                dummy_scorer = lDDTScorer(mdl_sel, self.compound_lib,
                                          symmetry_settings = sm,
                                          inclusion_radius = self.inclusion_radius,
                                          bb_only = self.bb_only)
                penalty += dummy_scorer.n_distances
        return penalty

    def _GetChainRNums(self, ch, residue_mapping, model_ch_name,
                       target_ch_name):
        """Map residues in model chain to target residues

        There are two options: one is simply using residue numbers,
        the other is a custom mapping as given in *residue_mapping*
        """
        if residue_mapping and model_ch_name in residue_mapping:
            # extract residue numbers from target chain
            ch_idx = self.chain_names.index(target_ch_name)
            start_idx = self.chain_res_start_indices[ch_idx]
            if ch_idx < len(self.chain_names) - 1:
                end_idx = self.chain_res_start_indices[ch_idx+1]
            else:
                end_idx = len(self.compound_names)
            target_rnums = self.res_resnums[start_idx:end_idx]
            # get sequences from alignment and do consistency checks
            target_seq = residue_mapping[model_ch_name].GetSequence(0)
            model_seq = residue_mapping[model_ch_name].GetSequence(1)
            if len(target_seq.GetGaplessString()) != len(target_rnums):
                raise RuntimeError(f"Try to perform residue mapping for "
                                   f"model chain {model_ch_name} which "
                                   f"maps to {target_ch_name} in target. "
                                   f"Target sequence in alignment suggests "
                                   f"{len(target_seq.GetGaplessString())} "
                                   f"residues but {len(target_rnums)} are "
                                   f"expected.")
            if len(model_seq.GetGaplessString()) != len(ch.residues):
                raise RuntimeError(f"Try to perform residue mapping for "
                                   f"model chain {model_ch_name} which "
                                   f"maps to {target_ch_name} in target. "
                                   f"Model sequence in alignment suggests "
                                   f"{len(model_seq.GetGaplessString())} "
                                   f"residues but {len(ch.residues)} are "
                                   f"expected.")
            rnums = list()
            target_idx = -1
            for col in residue_mapping[model_ch_name]:
                if col[0] != '-':
                    target_idx += 1
                # handle match
                if col[0] != '-' and col[1] != '-':
                    rnums.append(target_rnums[target_idx])
                # insertion in model adds None to rnum
                if col[0] == '-' and col[1] != '-':
                    rnums.append(None)
        else:
            rnums = [r.GetNumber() for r in ch.residues]

        return rnums


    def _SetupEnv(self, compound_lib, custom_compounds, symmetry_settings,
                  seqres_mapping, bb_only):
        """Sets target related lDDTScorer members defined in constructor

        No distance related members - see _SetupDistances
        """
        residue_numbers = self._GetTargetResidueNumbers(self.target,
                                                        seqres_mapping)
        current_idx = 0
        for chain in self.target.chains:
            ch_name = chain.GetName()
            self.chain_names.append(ch_name)
            self.chain_start_indices.append(current_idx)
            self.chain_res_start_indices.append(len(self.compound_names))
            for r, rnum in zip(chain.residues, residue_numbers[ch_name]):
                if r.name not in self.compound_anames:
                    # sets compound info in self.compound_anames and
                    # self.compound_symmetric_atoms
                    self._SetupCompound(r, compound_lib, custom_compounds,
                                        symmetry_settings, bb_only)

                self.res_start_indices.append(current_idx)
                self.res_mapper[(ch_name, rnum)] = len(self.compound_names)
                self.compound_names.append(r.name)
                self.res_resnums.append(rnum)

                atoms = [r.FindAtom(an) for an in self.compound_anames[r.name]]
                for a in atoms:
                    if a.IsValid():
                        self.atom_indices[a.handle.GetHashCode()] = current_idx
                    current_idx += 1
                
                if r.name in self.compound_symmetric_atoms:
                    for sym_tuple in self.compound_symmetric_atoms[r.name]:
                        for a_idx in sym_tuple:
                            a = atoms[a_idx]
                            if a.IsValid():
                                hashcode = a.handle.GetHashCode()
                                self.symmetric_atoms.add(
                                    self.atom_indices[hashcode]
                                )
        self.n_atoms = current_idx


    def _GetTargetResidueNumbers(self, target, seqres_mapping):
        """Returns residue numbers for each chain in target as dict

        They're either directly extracted from the raw residue number
        from the structure or from user provided alignments
        """
        residue_numbers = dict()
        for ch in target.chains:
            ch_name = ch.GetName()
            rnums = list()
            if ch_name in seqres_mapping:
                seqres = seqres_mapping[ch_name].GetSequence(0).GetString()
                atomseq = seqres_mapping[ch_name].GetSequence(1).GetString()
                # SEQRES must not contain gaps
                if "-" in seqres:
                    raise RuntimeError(
                        "SEQRES in seqres_mapping must not " "contain gaps"
                    )
                atomseq_from_chain = [r.one_letter_code for r in ch.residues]
                if atomseq.replace("-", "") != atomseq_from_chain:
                    raise RuntimeError(
                        "ATOMSEQ in seqres_mapping must match "
                        "raw sequence extracted from chain "
                        "residues"
                    )
                rnum = 0
                for seqres_olc, atomseq_olc in zip(seqres, atomseq):
                    if seqres_olc != "-":
                        rnum += 1
                    if atomseq_olc != "-":
                        if seqres_olc != atomseq_olc:
                            raise RuntimeError(
                                f"Residue with number {rnum} in "
                                f"chain {ch_name} has SEQRES "
                                f"ATOMSEQ mismatch"
                            )
                        rnums.append(mol.ResNum(rnum))
            else:
                rnums = [r.GetNumber() for r in ch.residues]
            assert len(rnums) == len(ch.residues)
            residue_numbers[ch_name] = rnums
        return residue_numbers

    def _SetupCompound(self, r, compound_lib, custom_compounds,
                       symmetry_settings, bb_only):
        """fill self.compound_anames/self.compound_symmetric_atoms
        """
        if bb_only:
            # throw away compound_lib info
            if r.chem_class.IsPeptideLinking():
                self.compound_anames[r.name] = ["CA"]
            elif r.chem_class.IsNucleotideLinking():
                self.compound_anames[r.name] = ["C3'"]
            else:
                raise RuntimeError(f"Only support amino acids and nucleotides "
                                   f"if bb_only is True, failed with {str(r)}")
            self.compound_symmetric_atoms[r.name] = list()
        else:
            atom_names = list()
            symmetric_atoms = list()
            if custom_compounds is not None and r.GetName() in custom_compounds:
                atom_names = list(custom_compounds[r.GetName()].atom_names)
            else:
                compound = compound_lib.FindCompound(r.name)
                if compound is None:
                    raise RuntimeError(f"no entry for {r} in compound_lib")
                for atom_spec in compound.GetAtomSpecs():
                    if atom_spec.element not in ["H", "D"]:
                        atom_names.append(atom_spec.name)
            if r.name in symmetry_settings.symmetric_compounds:
                for pair in symmetry_settings.symmetric_compounds[r.name]:
                    try:
                        a = atom_names.index(pair[0])
                        b = atom_names.index(pair[1])
                    except:
                        msg = f"Could not find symmetric atoms "
                        msg += f"({pair[0]}, {pair[1]}) for {r.name} "
                        msg += f"as specified in SymmetrySettings in "
                        msg += f"compound from component dictionary. "
                        msg += f"Atoms in compound: {atom_names}"
                        raise RuntimeError(msg)
                    symmetric_atoms.append((a, b))
            self.compound_anames[r.name] = atom_names
            if len(symmetric_atoms) > 0:
                self.compound_symmetric_atoms[r.name] = symmetric_atoms

    def _SetupDistances(self):
        """Compute distance related members of lDDTScorer
        """
        # init
        self._ref_indices = [np.asarray([], dtype=np.int64) for idx in range(self.n_atoms)]
        self._ref_distances = [np.asarray([], dtype=np.float64) for idx in range(self.n_atoms)]
        self._sym_ref_indices = [np.asarray([], dtype=np.int64) for idx in range(self.n_atoms)]
        self._sym_ref_distances = [np.asarray([], dtype=np.float64) for idx in range(self.n_atoms)]

        # initialize positions with values far in nirvana. If a position is not
        # set, it should be far away from any position in target (or at least
        # more than inclusion_radius).
        max_pos = self.target.bounds.GetMax()
        max_coordinate = abs(max(max_pos[0], max_pos[1], max_pos[2]))
        max_coordinate += 2 * self.inclusion_radius

        pos = np.ones((self.n_atoms, 3), dtype=np.float32) * max_coordinate
        atom_indices = list()
        mask_start = list()
        mask_end = list()

        for r_idx, r in enumerate(self.target.residues):
            r_start_idx = self.res_start_indices[r_idx]
            r_n_atoms = len(self.compound_anames[r.name])
            r_end_idx = r_start_idx + r_n_atoms
            for a in r.atoms:
                if a.handle.GetHashCode() in self.atom_indices:
                    idx = self.atom_indices[a.handle.GetHashCode()]
                    p = a.GetPos()
                    pos[idx][0] = p[0]
                    pos[idx][1] = p[1]
                    pos[idx][2] = p[2]
                    atom_indices.append(idx)
                    mask_start.append(r_start_idx)
                    mask_end.append(r_end_idx)

        indices, distances = self._CloseStuff(pos, self.inclusion_radius,
                                              atom_indices, mask_start,
                                              mask_end)

        for i in range(len(atom_indices)):
            self._ref_indices[atom_indices[i]] = indices[i]
            self._ref_distances[atom_indices[i]] = distances[i]
        self._NonSymDistances(self._ref_indices, self._ref_distances,
                              self._sym_ref_indices,
                              self._sym_ref_distances)

    def _SetupDistancesSC(self):
        """Select subset of contacts only covering intra-chain contacts
        """
        # init
        self._ref_indices_sc = [np.asarray([], dtype=np.int64) for idx in range(self.n_atoms)]
        self._ref_distances_sc = [np.asarray([], dtype=np.float64) for idx in range(self.n_atoms)]
        self._sym_ref_indices_sc = [np.asarray([], dtype=np.int64) for idx in range(self.n_atoms)]
        self._sym_ref_distances_sc = [np.asarray([], dtype=np.float64) for idx in range(self.n_atoms)]

        # start from overall contacts
        ref_indices = self.ref_indices
        ref_distances = self.ref_distances
        sym_ref_indices = self.sym_ref_indices
        sym_ref_distances = self.sym_ref_distances

        n_chains = len(self.chain_start_indices)
        for ch_idx, ch in enumerate(self.target.chains):
            chain_s = self.chain_start_indices[ch_idx]
            chain_e = self.n_atoms
            if ch_idx + 1 < n_chains:
                chain_e = self.chain_start_indices[ch_idx+1]
            for i in range(chain_s, chain_e):
                if len(ref_indices[i]) > 0:
                    intra_idx = np.where(np.logical_and(ref_indices[i]>=chain_s,
                                                  ref_indices[i]<chain_e))[0]
                    self._ref_indices_sc[i] = ref_indices[i][intra_idx]
                    self._ref_distances_sc[i] = ref_distances[i][intra_idx]

        self._NonSymDistances(self._ref_indices_sc, self._ref_distances_sc,
                              self._sym_ref_indices_sc,
                              self._sym_ref_distances_sc)

    def _SetupDistancesIC(self):
        """Select subset of contacts only covering inter-chain contacts
        """
        # init
        self._ref_indices_ic = [np.asarray([], dtype=np.int64) for idx in range(self.n_atoms)]
        self._ref_distances_ic = [np.asarray([], dtype=np.float64) for idx in range(self.n_atoms)]
        self._sym_ref_indices_ic = [np.asarray([], dtype=np.int64) for idx in range(self.n_atoms)]
        self._sym_ref_distances_ic = [np.asarray([], dtype=np.float64) for idx in range(self.n_atoms)]

        # start from overall contacts
        ref_indices = self.ref_indices
        ref_distances = self.ref_distances
        sym_ref_indices = self.sym_ref_indices
        sym_ref_distances = self.sym_ref_distances

        n_chains = len(self.chain_start_indices)
        for ch_idx, ch in enumerate(self.target.chains):
            chain_s = self.chain_start_indices[ch_idx]
            chain_e = self.n_atoms
            if ch_idx + 1 < n_chains:
                chain_e = self.chain_start_indices[ch_idx+1]
            for i in range(chain_s, chain_e):
                if len(ref_indices[i]) > 0:
                    inter_idx = np.where(np.logical_or(ref_indices[i]<chain_s,
                                                  ref_indices[i]>=chain_e))[0]
                    self._ref_indices_ic[i] = ref_indices[i][inter_idx]
                    self._ref_distances_ic[i] = ref_distances[i][inter_idx]

        self._NonSymDistances(self._ref_indices_ic, self._ref_distances_ic,
                              self._sym_ref_indices_ic,
                              self._sym_ref_distances_ic)

    def _CloseStuff(self, pos, inclusion_radius, indices, mask_start, mask_end):
        """returns close stuff for positions specified by indices
        """
        # TODO: this function does brute force distance computation which has
        # quadratic complexity...
        close_indices = list()
        distances = list()
        # work with squared_inclusion_radius (sir) to save some square roots
        sir = inclusion_radius ** 2
        for idx, ms, me in zip(indices, mask_start, mask_end):
            p = pos[idx, :]
            tmp = pos - p[None, :]
            np.square(tmp, out=tmp)
            tmp = tmp.sum(axis=1)
            # mask out atoms of own residue => put them far away
            tmp[range(ms, me)] = 2 * sir
            close_indices.append(np.nonzero(tmp <= sir)[0])
            distances.append(np.sqrt(tmp[close_indices[-1]]))
        return (close_indices, distances)

    def _NonSymDistances(self, ref_indices, ref_distances,
                         sym_ref_indices, sym_ref_distances):
        """Transfer indices/distances of non-symmetric atoms in place
        """
        for idx in self.symmetric_atoms:
            indices = list()
            distances = list()
            for i, d in zip(ref_indices[idx], ref_distances[idx]):
                if i not in self.symmetric_atoms:
                    indices.append(i)
                    distances.append(d)
            sym_ref_indices[idx] = indices
            sym_ref_distances[idx] = np.asarray(distances)

    def _EvalAtom(self, pos, atom_idx, thresholds, ref_indices, ref_distances):
        """Computes number of distance differences within given thresholds

        returns np.array with len(thresholds) elements
        """
        a_p = pos[atom_idx, :]
        tmp = pos.take(ref_indices[atom_idx], axis=0)
        np.subtract(tmp, a_p[None, :], out=tmp)
        np.square(tmp, out=tmp)
        tmp = tmp.sum(axis=1)
        np.sqrt(tmp, out=tmp)  # distances against all relevant atoms
        np.subtract(ref_distances[atom_idx], tmp, out=tmp)
        np.absolute(tmp, out=tmp)  # absolute dist diffs
        return np.asarray([(tmp <= thresh).sum() for thresh in thresholds],
                          dtype=np.int32)

    def _EvalAtoms(
        self, pos, atom_indices, thresholds, ref_indices, ref_distances
    ):
        """Calls _EvalAtom for several atoms and sums up the computed number
        of distance differences within given thresholds

        returns numpy matrix of shape (n_atoms, len(threshold))
        """
        conserved = np.zeros((len(atom_indices), len(thresholds)),
                             dtype=np.int32)
        for a_idx, a in enumerate(atom_indices):
            conserved[a_idx, :] = self._EvalAtom(pos, a, thresholds,
                                                 ref_indices, ref_distances)
        return conserved

    def _EvalResidues(self, pos, thresholds, res_atom_indices, ref_indices,
                      ref_distances):
        """Calls _EvalAtoms for a bunch of residues

        residues are defined in *res_atom_indices* as lists of atom indices
        returns numpy matrix of shape (n_residues, len(thresholds)).
        """
        conserved = np.zeros((len(res_atom_indices), len(thresholds)),
                             dtype=np.int32)
        for rai_idx, rai in enumerate(res_atom_indices):
            conserved[rai_idx,:] = np.sum(self._EvalAtoms(pos, rai, thresholds,
                                          ref_indices, ref_distances), axis=0)
        return conserved

    def _ProcessSequenceSeparation(self):
        if self.sequence_separation != 0:
            raise NotImplementedError("Congratulations! You're the first one "
                                      "requesting a non-default "
                                      "sequence_separation in the new and "
                                      "awesome lDDT implementation. A crate of "
                                      "beer for Gabriel and he'll implement "
                                      "it.")

    def _GetNExp(self, atom_idx, ref_indices):
        """Returns number of close atoms around one or several atoms
        """
        if isinstance(atom_idx, int):
            return len(ref_indices[atom_idx])
        elif isinstance(atom_idx, list):
            return sum([len(ref_indices[idx]) for idx in atom_idx])
        else:
            raise RuntimeError("invalid input type")

    def _ResolveSymmetries(self, pos, thresholds, symmetries, sym_ref_indices,
                           sym_ref_distances):
        """Swaps symmetric positions in-place in order to maximize lDDT scores
        towards non-symmetric atoms.
        """
        for sym in symmetries:

            atom_indices = list()
            for sym_tuple in sym:
                atom_indices += [sym_tuple[0], sym_tuple[1]]
            tot = self._GetNExp(atom_indices, sym_ref_indices)

            if tot == 0:
                continue  # nothing to do

            # score as is
            sym_one_conserved = self._EvalAtoms(
                pos,
                atom_indices,
                thresholds,
                sym_ref_indices,
                sym_ref_distances,
            )

            # switch positions and score again
            for pair in sym:
                pos[[pair[0], pair[1]]] = pos[[pair[1], pair[0]]]

            sym_two_conserved = self._EvalAtoms(
                pos,
                atom_indices,
                thresholds,
                sym_ref_indices,
                sym_ref_distances,
            )

            sym_one_score = np.sum(sym_one_conserved) / (len(thresholds) * tot)
            sym_two_score = np.sum(sym_two_conserved) / (len(thresholds) * tot)

            if sym_one_score >= sym_two_score:
                # switch back, initial positions were better or equal
                # for the equal case: we still switch back to reproduce the old
                # lDDT behaviour
                for pair in sym:
                    pos[[pair[0], pair[1]]] = pos[[pair[1], pair[0]]]
