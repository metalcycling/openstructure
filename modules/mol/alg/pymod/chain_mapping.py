import itertools
import copy

import numpy as np

from scipy.special import factorial
from scipy.special import binom # as of Python 3.8, the math module implements
                                # comb, i.e. n choose k

from ost import seq
from ost import mol
from ost import geom

from ost.mol.alg import lddt

class ChainMapper:
    def __init__(self, target, resnum_alignments=False):
        """Object to compute chain mappings

        :param target: Target structure onto which models are mapped.
                       Computations happen on a selection only containing
                       polypeptides and polynucleotides.
        :type target: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param resnum_alignments: Use residue numbers instead of
                                  Needleman-Wunsch to compute pairwise
                                  alignments. Relevant for :attr:`~chem_groups` 
                                  and related attributes.
        :type resnum_alignments: :class:`bool`        
        """

        # target structure preprocessing
        self.target, self.polypep_seqs, self.polynuc_seqs = \
        _ProcessStructure(target)

        # helper class to generate pairwise alignments
        self.aligner = _Aligner(resnum_aln = resnum_alignments)

        # lazy computed attributes
        self._chem_groups = None
        self._chem_group_alignments = None
        self._chem_group_ref_seqs = None
        self._chem_group_types = None

    @property
    def chem_groups(self):
        """Groups of chemically equivalent poly-peptides/nucleotides in *target*

        Seq. id. > 95% is considered equivalent. First chain in group is the one
        with longest sequence.
      
        :getter: Computed on first use (cached)
        :type: :class:`list` of :class:`list` of :class:`str` (chain names)
        """
        if self._chem_groups is None:
            self._chem_groups = list()
            for a in self.chem_group_alignments:
                self._chem_groups.append([s.GetName() for s in a.sequences])
        return self._chem_groups
    
    @property
    def chem_group_alignments(self):
        """MSA for each group in :attr:`~chem_groups`

        Sequences in MSAs exhibit same order as in :attr:`~chem_groups`

        :getter: Computed on first use (cached)
        :type: :class:`ost.seq.AlignmentList`
        """
        if self._chem_group_alignments is None:
            self._chem_group_alignments, self._chem_group_types = \
            _GetChemGroupAlignments(self.polypep_seqs, self.polynuc_seqs,
                                    self.aligner)
        return self._chem_group_alignments

    @property
    def chem_group_ref_seqs(self):
        """Reference (longest) sequence for each group in :attr:`~chem_groups`

        :getter: Computed on first use (cached)
        :type: :class:`ost.seq.SequenceList`
        """
        if self._chem_group_ref_seqs is None:
            self._chem_group_ref_seqs = seq.CreateSequenceList()
            for a in self.chem_group_alignments:
                s = seq.CreateSequence(a.GetSequence(0).GetName(),
                                       a.GetSequence(0).GetGaplessString())
                if a.GetSequence(0).HasAttachedView():
                    s.AttachView(a.GetSequence(0).GetAttachedView())
                self._chem_group_ref_seqs.AddSequence(s)
        return self._chem_group_ref_seqs

    @property
    def chem_group_types(self):
        """ChemType of each group in :attr:`~chem_groups`

        Specifying if groups are poly-peptides/nucleotides, i.e. 
        :class:`ost.mol.ChemType.AMINOACIDS` or
        :class:`ost.mol.ChemType.NUCLEOTIDES`
        
        :getter: Computed on first use (cached)
        :type: :class:`list` of :class:`ost.mol.ChemType`
        """
        if self._chem_group_types is None:
            self._chem_group_alignments, self._chem_group_types = \
            _GetChemGroupAlignments(self.polypep_seqs, self.polynuc_seqs,
                                    self.aligner)            
        return self._chem_group_types
        
    def GetChemMapping(self, model):
        """Maps sequences in *model* to chem_groups of target

        :param model: Model from which to extract sequences, a
                      selection that only includes peptides and nucleotides
                      is performed and returned along other results.
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :returns: Tuple with two lists of length `len(self.chem_groups)` and
                  an :class:`ost.mol.EntityView` representing *model*:
                  1) Each element is a :class:`list` with mdl chain names that
                  map to the chem group at that position.
                  2) Each element is a :class:`ost.seq.AlignmentList` aligning
                  these mdl chain sequences to the chem group ref sequences.
                  3) A selection of *model* that only contains polypeptides and
                  polynucleotides whose ATOMSEQ exactly matches the sequence
                  info in the returned alignments.
        """
        mdl, mdl_pep_seqs, mdl_nuc_seqs = _ProcessStructure(model)

        mapping = [list() for x in self.chem_groups]
        alns = [seq.AlignmentList() for x in self.chem_groups]

        for s in mdl_pep_seqs:
            idx, aln = _MapSequence(self.chem_group_ref_seqs, 
                                    self.chem_group_types,
                                    s, mol.ChemType.AMINOACIDS,
                                    95., 0.1, self.aligner)
            if idx is not None:
                mapping[idx].append(s.GetName())
                alns[idx].append(aln)

        for s in mdl_nuc_seqs:
            idx, aln = _MapSequence(self.chem_group_ref_seqs, 
                                    self.chem_group_types,
                                    s, mol.ChemType.NUCLEOTIDES,
                                    95., 0.1, self.aligner)
            if idx is not None:
                mapping[idx].append(s.GetName())
                alns[idx].append(aln)

        return (mapping, alns, mdl)


    def GetNaivelDDTMapping(self, model, bb_only=False, inclusion_radius=15.0,
                            thresholds=[0.5, 1.0, 2.0, 4.0]):
        """Naively iterates all possible chain mappings and returns the best

        Maps *model* chain sequences to :attr:`~chem_groups` and performs all
        possible permutations. The best mapping is selected based on lDDT score.

        :param model: Model to map
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param bb_only: lDDT considers only atoms of name "CA" (peptides) or
                        "C3'" (nucleotides). Gives speed improvement but
                        sidechains are not considered anymore.
        :type bb_only: :class:`bool`
        :param inclusion_radius: Inclusion radius for lDDT
        :type inclusion_radius: :class:`float`
        :param thresholds: Thresholds for lDDT
        :type thresholds: :class:`list` of :class:`float`
        :returns: A :class:`list` of :class:`list` that reflects
                  :attr:`~chem_groups` but is filled with the respective model
                  chains. Target chains without mapped model chains are set to
                  None.
        """
        chem_mapping, chem_group_alns, mdl = self.GetChemMapping(model)

        # check for the simplest case
        one_to_one = _CheckOneToOneMapping(self.chem_groups, chem_mapping)
        if one_to_one is not None:
            return one_to_one

        # all possible ref/mdl alns given chem mapping
        ref_mdl_alns =  _GetRefMdlAlns(self.chem_groups,
                                       self.chem_group_alignments,
                                       chem_mapping,
                                       chem_group_alns)

        # Setup scoring
        lddt_scorer = lddt.lDDTScorer(self.target, bb_only = bb_only)
        best_mapping = None
        best_lddt = -1.0

        for mapping in _ChainMappings(self.chem_groups, chem_mapping):
            # chain_mapping and alns as input for lDDT computation
            lddt_chain_mapping = dict()
            lddt_alns = dict()

            for ref_chem_group, mdl_chem_group, ref_aln in \
            zip(self.chem_groups, mapping, self.chem_group_alignments):
                for ref_ch, mdl_ch in zip(ref_chem_group, mdl_chem_group):
                    # some mdl chains can be None
                    if mdl_ch is not None:
                        lddt_chain_mapping[mdl_ch] = ref_ch
                        lddt_alns[mdl_ch] = ref_mdl_alns[(ref_ch, mdl_ch)]

            lDDT, _ = lddt_scorer.lDDT(mdl, thresholds=thresholds,
                                       chain_mapping=lddt_chain_mapping,
                                       residue_mapping = lddt_alns,
                                       check_resnames = False)

            if lDDT > best_lddt:
                best_mapping = mapping
                best_lddt = lDDT

        return best_mapping

    def GetDecomposerlDDTMapping(self, model, inclusion_radius=15.0,
                                 thresholds=[0.5, 1.0, 2.0, 4.0]):
        """Naively iterates all possible chain mappings and returns the best

        Maps *model* chain sequences to :attr:`~chem_groups` and performs all
        possible permutations. The best mapping is selected based on lDDT score
        which is computed only on backbone atoms (CA for amino acids C3' for
        Nucleotides).

        :param model: Model to map
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param inclusion_radius: Inclusion radius for lDDT
        :type inclusion_radius: :class:`float`
        :param thresholds: Thresholds for lDDT
        :type thresholds: :class:`list` of :class:`float`
        :returns: A :class:`list` of :class:`list` that reflects
                  :attr:`~chem_groups` but is filled with the respective model
                  chains. Target chains without mapped model chains are set to
                  None.
        """
        chem_mapping, chem_group_alns, mdl = self.GetChemMapping(model)

        # check for the simplest case
        one_to_one = _CheckOneToOneMapping(self.chem_groups, chem_mapping)
        if one_to_one is not None:
            return one_to_one

        ref_mdl_alns =  _GetRefMdlAlns(self.chem_groups,
                                       self.chem_group_alignments,
                                       chem_mapping,
                                       chem_group_alns)

        # Setup scoring
        lddt_scorer = _lDDTDecomposer(self.target, mdl, ref_mdl_alns,
                                      inclusion_radius=inclusion_radius,
                                      thresholds = thresholds)
        best_mapping = None
        best_lddt = -1.0

        for mapping in _ChainMappings(self.chem_groups, chem_mapping):

            lDDT = lddt_scorer.lDDT(self.chem_groups, mapping)
            if lDDT > best_lddt:
                best_mapping = mapping
                best_lddt = lDDT

        return best_mapping


    def GetGreedylDDTMapping(self, model, inclusion_radius=15.0,
                             thresholds=[0.5, 1.0, 2.0, 4.0],
                             seed_strategy="fast", steep_opt_rate = None,
                             full_n_mdl_chains = None, block_seed_size = None,
                             block_n_mdl_chains = None):
        """Heuristic to lower the complexity of naive iteration

        Maps *model* chain sequences to :attr:`~chem_groups` and extends these
        start mappings (seeds) in a greedy way. In each iteration, the
        one-to-one mapping that leads to highest increase in number of conserved
        contacts is added with the additional requirement that this added
        mapping must have non-zero interface counts towards the already mapped
        chains. So basically we're "growing" the mapped structure by only adding
        connected stuff.

        Several strategies exist to identify the start seed(s):

        * fast: perform all vs. all single chain lDDTs within the respective
          ref/mdl chem groups. The mapping with highest number of conserved
          contacts is selected as seed for extension

        * full: try multiple seeds, i.e. try all ref/mdl chain combinations
          within the respective chem groups and retain the mapping leading to
          the best lDDT. Optionally, you can reduce the number of mdl chains
          per ref chain to the *full_n_mdl_chains* best scoring ones.

        * block: try multiple seeds, i.e. try all ref/mdl chain combinations
          within the respective chem groups but only extend these seeds by
          *block_seed_size* chains. The highest scoring block for every ref
          chain is extended exhaustively to identify the best scoring initial
          block.

        :param model: Model to map
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param inclusion_radius: Inclusion radius for lDDT
        :type inclusion_radius: :class:`float`
        :param thresholds: Thresholds for lDDT
        :type thresholds: :class:`list` of :class:`float`
        :param seed_strategy: Strategy to pick starting seeds for expansion
        :type seed_strategy: :class:`str`
        :param steep_opt_rate: If set, every *steep_opt_rate* mappings, a simple
                               optimization is executed with the goal of
                               avoiding local minima. The optimization
                               iteratively checks all possible swaps of mappings
                               within their respective chem groups and accepts
                               swaps that improve lDDT score. Iteration stops as
                               soon as no improvement can be achieved anymore.
        :type stepp_opt_rate: :class:`int`
        :param full_n_mdl_chains: Param for *full* seed strategy - Max number of
                                  mdl chains that are tried per ref chain. The
                                  default (None) tries all of them.
        :type full_n_mdl_chains: :class:`int`
        :param block_seed_size: Param for *block* seed strategy - Initial seeds
                                are extended by that number of chains. The
                                default (None) performs full extensions and you
                                get equivalent behaviour as in *full* strategy.
        :type block_seed_size: :class:`int`
        :param block_n_mdl_chains: Equivalent of *full_n_mdl_chains* but for
                                   *block* seed strategy.
        :type block_n_mdl_chains: :class:`int`
        :returns: A :class:`list` of :class:`list` that reflects
                  :attr:`~chem_groups` but is filled with the respective model
                  chains. Target chains without mapped model chains are set to
                  None.
        """

        seed_strategies = ["fast", "full", "block"]
        if seed_strategy not in seed_strategies:
            raise RuntimeError(f"Seed strategy must be in {seed_strategies}")

        chem_mapping, chem_group_alns, mdl = self.GetChemMapping(model)

        # check for the simplest case
        only_one_to_one = True
        for ref_chains, mdl_chains in zip(self.chem_groups, chem_mapping):
            if len(ref_chains) != 1 or len(mdl_chains) not in [0, 1]:
                only_one_to_one = False
                break
        if only_one_to_one:
            # skip ref chem groups with no mapped mdl chain
            return [(a,b) for a,b in zip(self.chem_groups, chem_mapping) if len(b) == 1]

        ref_mdl_alns =  _GetRefMdlAlns(self.chem_groups,
                                       self.chem_group_alignments,
                                       chem_mapping,
                                       chem_group_alns)

        # setup greedy searcher
        the_greed = _GreedySearcher(self.target, mdl, self.chem_groups,
                                    chem_mapping, ref_mdl_alns,
                                    inclusion_radius=inclusion_radius,
                                    thresholds=thresholds,
                                    steep_opt_rate=steep_opt_rate)

        if seed_strategy == "fast":
            return _FastGreedy(the_greed)
        elif seed_strategy == "full":
            return _FullGreedy(the_greed, full_n_mdl_chains)
        elif seed_strategy == "block":
            return _BlockGreedy(the_greed, block_seed_size, block_n_mdl_chains)


    def GetRigidMapping(self, model, single_chain_gdtts_thresh=0.5,
                        subsampling=None):
        """Identify chain mapping based on rigid superposition

        Superposition and scoring is based on CA/C3' positions which are present
        in all chains of a :attr:`ChainMapper.chem_group` as well as the *model*
        chains which are mapped to that respective chem group.

        Transformations to superpose *model* onto :attr:`ChainMapper.target`
        are estimated using all possible combinations of target and model chains
        within the same chem groups. For each transformation, the mapping is
        extended by iteratively adding the model/target chain pair that adds
        the most conserved contacts based on the GDT-TS metric (Number of
        CA/C3' atoms within [8, 4, 2, 1] Angstrom). The mapping with highest
        GDT-TS score is returned. However, that mapping is not guaranteed to be
        complete (see *single_chain_gdtts_thresh*).

        :param model: Model to map
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param single_chain_gdtts_thresh: Minimal GDT-TS score for model/target
                                          chain pair to be added to mapping.
                                          Mapping extension for a given
                                          transform stops when no pair fulfills
                                          this threshold, potentially leading to
                                          an incomplete mapping.
        :type single_chain_gdtts_thresh: :class:`float`
        :param subsampling: If given, only use an equally distributed subset
                            of all CA/C3' positions for superposition/scoring.
        :type subsampling: :class:`int`
        :returns: A :class:`list` of :class:`list` that reflects
                  :attr:`~chem_groups` but is filled with the respective model
                  chains. Target chains without mapped model chains are set to
                  None.
        """

        chem_mapping, chem_group_alns, mdl = self.GetChemMapping(model)

        trg_group_pos, mdl_group_pos = _GetRefPos(self.target, mdl,
                                                  self.chem_group_alignments,
                                                  chem_group_alns,
                                                  max_pos = subsampling)

        # get transforms of any mdl chain onto any trg chain in same chem group
        transforms = list()
        for trg_pos, mdl_pos in zip(trg_group_pos, mdl_group_pos):
            for t in trg_pos:
                for m in mdl_pos:
                    if len(t) >= 3 and len(m) >= 3:
                        assert(len(t) == len(m))
                        res = mol.alg.SuperposeSVD(m, t)
                        transforms.append(res.transformation)

        best_mapping = None
        best_gdt = 0
        for transform in transforms:
            mapping = dict()
            gdt_cache = dict() # to cache non-normalized gdt scores
            for trg_pos, trg_chains, mdl_pos, mdl_chains in zip(trg_group_pos,
                                                                self.chem_groups,
                                                                mdl_group_pos,
                                                                chem_mapping):
                if len(trg_pos) == 0 or len(trg_pos[0]) == 0:
                    continue # cannot compute valid gdt

                n_contacts = 4 * len(trg_pos[0]) # for gdt normalization
                mapped_mdl_chains = set()
                while True:
                    best_sc_mapping = None
                    best_sc_gdt = 0.0
                    for m_pos, m in zip(mdl_pos, mdl_chains):
                        if m not in mapped_mdl_chains:
                            for t_pos, t in zip(trg_pos, trg_chains):
                                if t not in mapping:
                                    if (t,m) in gdt_cache:
                                        gdt = gdt_cache[(t,m)]
                                    else:
                                        t_m_pos = geom.Vec3List(m_pos)
                                        t_m_pos.ApplyTransform(transform)
                                        gdt = t_pos.GetGDTTS(t_m_pos, norm=False)
                                        gdt_cache[(t,m)] = gdt
                                    if gdt > best_sc_gdt:
                                        best_sc_gdt = gdt
                                        best_sc_mapping = (t, m)
                    if best_sc_gdt/n_contacts > single_chain_gdtts_thresh:
                        mapping[best_sc_mapping[0]] = best_sc_mapping[1]
                        mapped_mdl_chains.add(best_sc_mapping[1])
                    else:
                        break

            # compute overall gdt for this transform
            gdt = 0
            for t,m in mapping.items():
                gdt += gdt_cache[(t,m)]

            if gdt > best_gdt:
                best_gdt = gdt
                best_mapping = mapping

        # translate mapping format and return
        final_mapping = list()
        for ref_chains in self.chem_groups:
            mapped_mdl_chains = list()
            for ref_ch in ref_chains:
                if ref_ch in mapping:
                    mapped_mdl_chains.append(best_mapping[ref_ch])
                else:
                    mapped_mdl_chains.append(None)
            final_mapping.append(mapped_mdl_chains)

        return final_mapping


    def GetNMappings(self, model):
        """ Returns number of possible mappings

        :param model: Model with chains that are mapped onto
                      :attr:`ChainMapper.chem_groups`
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        """
        chem_mapping, chem_group_alns, mdl = self.GetChemMapping(model)
        return _NMappings(self.chem_groups, chem_mapping)


def _FastGreedy(the_greed):

    something_happened = True
    mapping = dict()

    while something_happened:
        something_happened = False
        # search for best scoring starting point
        n_best = 0
        best_seed = None
        mapped_ref_chains = set(mapping.keys())
        mapped_mdl_chains = set(mapping.values())
        for ref_chains, mdl_chains in zip(the_greed.ref_chem_groups,
                                          the_greed.mdl_chem_groups):
            for ref_ch in ref_chains:
                if ref_ch not in mapped_ref_chains:
                    for mdl_ch in mdl_chains:
                        if mdl_ch not in mapped_mdl_chains:
                            n = the_greed.SCCounts(ref_ch, mdl_ch)
                            if n > n_best:
                                n_best = n
                                best_seed = (ref_ch, mdl_ch)
        if n_best == 0:
            break # no proper seed found anymore...
        # add seed to mapping and start the greed
        mapping[best_seed[0]] = best_seed[1]
        mapping = the_greed.ExtendMapping(mapping)
        something_happened = True


    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping


def _FullGreedy(the_greed, n_mdl_chains):
    """ Uses each reference chain as starting point for expansion

    However, not all mdl chain are mapped onto these reference chains,
    that's controlled by *n_mdl_chains*
    """

    if n_mdl_chains is not None and n_mdl_chains < 1:
        raise RuntimeError("n_mdl_chains must be None or >= 1")

    something_happened = True
    mapping = dict()

    while something_happened:
        something_happened = False
        # Try all possible starting points and keep the one giving the best lDDT
        best_lddt = 0.0
        best_mapping = None
        mapped_ref_chains = set(mapping.keys())
        mapped_mdl_chains = set(mapping.values())
        for ref_chains, mdl_chains in zip(the_greed.ref_chem_groups,
                                          the_greed.mdl_chem_groups):
            for ref_ch in ref_chains:
                if ref_ch not in mapped_ref_chains:
                    seeds = list()
                    for mdl_ch in mdl_chains:
                        if mdl_ch not in mapped_mdl_chains:
                            seeds.append((ref_ch, mdl_ch))
                    if n_mdl_chains is not None and n_mdl_chains < len(seeds):
                        counts = [the_greed.SCCounts(s[0], s[1]) for s in seeds]
                        tmp = [(a,b) for a,b in zip(counts, seeds)]
                        tmp.sort(reverse=True)
                        seeds = [item[1] for item in tmp[:n_mdl_chains]]
                    for seed in seeds:
                        tmp_mapping = dict(mapping)
                        tmp_mapping[seed[0]] = seed[1]
                        tmp_mapping = the_greed.ExtendMapping(tmp_mapping)
                        tmp_lddt = the_greed.lDDTFromFlatMap(tmp_mapping)
                        if tmp_lddt > best_lddt:
                            best_lddt = tmp_lddt
                            best_mapping = tmp_mapping

        if best_lddt == 0.0:
            break # no proper mapping found anymore...

        something_happened = True
        mapping = best_mapping

    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping


def _BlockGreedy(the_greed, seed_size, n_mdl_chains):
    """ Uses each reference chain as starting point for expansion

    Tries to map all mdl chains (optionally up to *n_mdl_chains* best ones)
    to these references but initially does not perform full expansion but only
    up to *seed_size*. The best scoring block for each reference is then used
    for full expansion.
    """

    if seed_size is not None and seed_size < 1:
        raise RuntimeError("seed_size must be None or >= 1")

    if n_mdl_chains is not None and n_mdl_chains < 1:
        raise RuntimeError("n_mdl_chains must be None or >= 1")

    ref_chem_groups = copy.deepcopy(the_greed.ref_chem_groups)
    mdl_chem_groups = copy.deepcopy(the_greed.mdl_chem_groups)

    mapping = dict()

    something_happened = True
    while something_happened:
        something_happened = False

        # one block per ref chain, i.e. a mapping that is extended by seed_size
        starting_blocks = dict()
        for ref_chains, mdl_chains in zip(ref_chem_groups, mdl_chem_groups):
            if len(mdl_chains) == 0:
                continue # nothing to map
            for ref_ch in ref_chains:
                best_lddt = 0.0
                best_mapping = None
                seeds = [(ref_ch, mdl_ch) for mdl_ch in mdl_chains]
                if n_mdl_chains is not None and n_mdl_chains < len(seeds):
                    counts = [the_greed.SCCounts(s[0], s[1]) for s in seeds]
                    tmp = [(a,b) for a,b in zip(counts, seeds)]
                    tmp.sort(reverse=True)
                    seeds = [item[1] for item in tmp[:n_mdl_chains]]
                for s in seeds:
                    seed = dict(mapping)
                    seed.update({s[0]: s[1]})
                    # max_ext = seed_size-1 => start seed already has size 1
                    seed = the_greed.ExtendMapping(seed, max_ext = seed_size-1)
                    seed_lddt = the_greed.lDDTFromFlatMap(seed)
                    if seed_lddt > best_lddt:
                        best_lddt = seed_lddt
                        best_mapping = seed
                starting_blocks[ref_ch] = best_mapping

        best_lddt = 0.0
        best_mapping = None
        for ref_ch, seed in starting_blocks.items():
            seed = the_greed.ExtendMapping(seed)
            seed_lddt = the_greed.lDDTFromFlatMap(seed)
            if seed_lddt > best_lddt:
                best_lddt = seed_lddt
                best_mapping = seed

        if best_lddt == 0.0:
            break # no proper mapping found anymore

        something_happened = True
        mapping.update(best_mapping)
        for ref_ch, mdl_ch in best_mapping.items():
            for group_idx in range(len(ref_chem_groups)):
                if ref_ch in ref_chem_groups[group_idx]:
                    ref_chem_groups[group_idx].remove(ref_ch)
                if mdl_ch in mdl_chem_groups[group_idx]:
                    mdl_chem_groups[group_idx].remove(mdl_ch)

    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping


def _GetRefPos(trg, mdl, trg_msas, mdl_alns, max_pos = None):
    """ Extracts reference positions which are present in trg and mdl
    """

    # select only backbone atoms, makes processing simpler later on
    # (just select res.atoms[0].GetPos() as ref pos)
    bb_trg = trg.Select("aname=\"CA\",\"C3'\"")
    bb_mdl = mdl.Select("aname=\"CA\",\"C3'\"")

    # mdl_alns are pairwise, let's construct MSAs
    mdl_msas = list()
    for aln_list in mdl_alns:
        if len(aln_list) > 0:
            tmp = aln_list[0].GetSequence(0)
            ref_seq = seq.CreateSequence(tmp.GetName(), tmp.GetGaplessString())
            mdl_msas.append(seq.alg.MergePairwiseAlignments(aln_list, ref_seq))
        else:
            mdl_msas.append(seq.CreateAlignment())

    trg_pos = list()
    mdl_pos = list()

    for trg_msa, mdl_msa in zip(trg_msas, mdl_msas):

        # make sure they have the same ref sequence (should be a given...)
        assert(trg_msa.GetSequence(0).GetGaplessString() == \
               mdl_msa.GetSequence(0).GetGaplessString())

        # check which columns in MSAs are fully covered (indices relative to
        # first sequence)
        trg_indices = _GetFullyCoveredIndices(trg_msa)
        mdl_indices = _GetFullyCoveredIndices(mdl_msa)

        # get indices where both, mdl and trg, are fully covered
        indices = sorted(list(trg_indices.intersection(mdl_indices)))

        # subsample if necessary
        if max_pos is not None and len(indices) > max_pos:
            step = int(len(indices)/max_pos)
            indices = [indices[i] for i in range(0, len(indices), step)]

        # translate to column indices in the respective MSAs
        trg_indices = _RefIndicesToColumnIndices(trg_msa, indices)
        mdl_indices = _RefIndicesToColumnIndices(mdl_msa, indices)

        # extract positions
        trg_pos.append(list())
        mdl_pos.append(list())
        for s_idx in range(trg_msa.GetCount()):
            trg_pos[-1].append(_ExtractMSAPos(trg_msa, s_idx, trg_indices,
                                              bb_trg))
        # first seq in mdl_msa is ref sequence in trg and does not belong to mdl
        for s_idx in range(1, mdl_msa.GetCount()):
            mdl_pos[-1].append(_ExtractMSAPos(mdl_msa, s_idx, mdl_indices,
                                              bb_mdl))

    return (trg_pos, mdl_pos)

def _GetFullyCoveredIndices(msa):
    """ Helper for _GetRefPos

    Returns a set containing the indices relative to first sequence in msa which
    are fully covered in all other sequences

    --AA-A-A
    -BBBB-BB
    CCCC-C-C

    => (0,1,3)
    """
    indices = set()
    ref_idx = 0
    for col in msa:
        if sum([1 for olc in col if olc != '-']) == col.GetRowCount():
            indices.add(ref_idx)
        if col[0] != '-':
            ref_idx += 1
    return indices

def _RefIndicesToColumnIndices(msa, indices):
    """ Helper for _GetRefPos

    Returns a list of mapped indices. indices refer to non-gap one letter
    codes in the first msa sequence. The returnes mapped indices are translated
    to the according msa column indices
    """
    ref_idx = 0
    mapping = dict()
    for col_idx, col in enumerate(msa):
        if col[0] != '-':
            mapping[ref_idx] = col_idx
            ref_idx += 1
    return [mapping[i] for i in indices]

def _ExtractMSAPos(msa, s_idx, indices, view):
    """ Helper for _GetRefPos

    Returns a geom.Vec3List containing positions refering to given msa sequence.
    => Chain with corresponding name is mapped onto sequence and the position of
    the first atom of each residue specified in indices is extracted.
    Indices refers to column indices in msa!
    """
    s = msa.GetSequence(s_idx)
    s_v = view.Select(f"cname={s.GetName()}")

    # sanity check
    assert(len(s.GetGaplessString()) == len(s_v.residues))

    residue_idx = [s.GetResidueIndex(i) for i in indices]
    return geom.Vec3List([s_v.residues[i].atoms[0].pos for i in residue_idx])

def _CheckOneToOneMapping(ref_chains, mdl_chains):
    """ Checks whether we already have a perfect one to one mapping

    That means each list in *ref_chains* has exactly one element and each
    list in *mdl_chains* has either one element (it's mapped) or is empty
    (ref chain has no mapped mdl chain). Returns None if no such mapping
    can be found.

    :param ref_chains: corresponds to :attr:`ChainMapper.chem_groups`
    :type ref_chains: :class:`list` of :class:`list` of :class:`str`
    :param mdl_chains: mdl chains mapped to chem groups in *ref_chains*, i.e.
                       the return value of :func:`ChainMapper.GetChemMapping`
    :type mdl_chains: class:`list` of :class:`list` of :class:`str`
    :returns: A :class:`list` of :class:`list` if a one to one mapping is found,
              None otherwise
    """
    only_one_to_one = True
    one_to_one = list()
    for ref, mdl in zip(ref_chains, mdl_chains):
        if len(ref) == 1 and len(mdl) == 1:
            one_to_one.append(mdl[0])
        elif len(ref) == 1 and len(mdl) == 0:
            one_to_one.append(None)
        else:
            only_one_to_one = False
            break
    if only_one_to_one:
        return one_to_one
    else:
        return None

def _GetChemGroupAlignments(pep_seqs, nuc_seqs, aligner, seqid_thr=95.,
                            gap_thr=0.1):
    """Returns alignments with groups of chemically equivalent chains

    :param pep_seqs: List of polypeptide sequences
    :type pep_seqs: :class:`seq.SequenceList`
    :param nuc_seqs: List of polynucleotide sequences
    :type nuc_seqs: :class:`seq.SequenceList` 
    :param seqid_thr: Threshold used to decide when two chains are identical.
                      95 percent tolerates the few mutations crystallographers
                      like to do.
    :type seqid_thr:  :class:`float`
    :param gap_thr: Additional threshold to avoid gappy alignments with high
                    seqid. The reason for not just normalizing seqid by the
                    longer sequence is that one sequence might be a perfect
                    subsequence of the other but only cover half of it. This
                    threshold checks for a maximum allowed fraction of gaps
                    in any of the two sequences after stripping terminal gaps.
    :type gap_thr: :class:`float`
    :param aligner: Helper class to generate pairwise alignments
    :type aligner: :class:`_Aligner`
    :returns: Tuple with first element being an AlignmentList. Each alignment
              represents a group of chemically equivalent chains and the first
              sequence is the longest. Second element is a list of equivalent
              length specifying the types of the groups. List elements are in
              [:class:`ost.ChemType.AMINOACIDS`,
              :class:`ost.ChemType.NUCLEOTIDES`] 
    """
    pep_groups = _GroupSequences(pep_seqs, seqid_thr, gap_thr, aligner,
                                 mol.ChemType.AMINOACIDS)
    nuc_groups = _GroupSequences(nuc_seqs, seqid_thr, gap_thr, aligner,
                                 mol.ChemType.NUCLEOTIDES)
    group_types = [mol.ChemType.AMINOACIDS] * len(pep_groups)
    group_types += [mol.ChemType.NUCLEOTIDES] * len(nuc_groups)
    groups = pep_groups
    groups.extend(nuc_groups)
    return (groups, group_types)

def _ProcessStructure(ent):
    """ Entity processing for chain mapping

    * Selects view containing peptide and nucleotide residues
    * Extracts atom sequences for each chain in that view
    * Attaches corresponding :class:`ost.mol.EntityView` to each sequence
    * Ensures strictly increasing residue numbers without insertion codes
      in each chain

    :param ent: Entity to process
    :type ent: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :returns: Tuple with 3 elements: 1) :class:`EntityView` containing
              peptide and nucleotide residues 2) :class:`seq.SequenceList`
              containing ATOMSEQ sequences for each polypeptide chain in
              returned view, sequences have :class:`ost.mol.EntityView` of
              according chains attached 3) same for polynucleotide chains
    """
    view = _StructureSelection(ent)
    polypep_seqs, polynuc_seqs = _GetAtomSeqs(view)

    for s in polypep_seqs:
        _AttachSeq(view, s)

    for s in polynuc_seqs:
        _AttachSeq(view, s)

    return (view, polypep_seqs, polynuc_seqs)

def _AttachSeq(view, s):
    """ Helper for _ProcessStructure

    Attaches view of chain with name s.GetName() to s and does sanity checks on
    residue numbers in that view
    """
    chain_view = view.Select(f"cname={s.GetName()}")

    # check if no insertion codes are present in residue numbers
    ins_codes = [r.GetNumber().GetInsCode() for r in chain_view.residues]
    if len(set(ins_codes)) != 1 or ins_codes[0] != '\0':
        raise RuntimeError("Residue numbers in input structures must not "
                           "contain insertion codes")

    # check if residue numbers are strictly increasing
    nums = [r.GetNumber().GetNum() for r in chain_view.residues]
    if not all(i < j for i, j in zip(nums, nums[1:])):
        raise RuntimeError("Residue numbers in input structures must be "
                           "strictly increasing for each chain")

    # attach
    s.AttachView(chain_view)

def _StructureSelection(ent):
    """ Helper for _ProcessStructure

    Selects structure to only contain peptide and nucleotide residues

    Additionally selects for residues that also contain the backbone
    representatives (CA for peptide and C3' for nucleotides) to avoid ATOMSEQ
    alignment inconsistencies when switching between all atom and backbone only
    representations.
    """
    query = "peptide=true or nucleotide=true and aname=\"CA\",\"C3'\""
    sel = ent.Select(query)
    view = ent.CreateEmptyView()
    for r in sel.residues:
        view.AddResidue(r.handle, mol.INCLUDE_ALL)
    return view

def _GetAtomSeqs(ent):
    """Extracts and returns atomseqs for polypeptide/nucleotide chains in ent

    :param ent: Entity for which you want to extract atomseqs
    :type ent: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :returns: Two lists, first containing atomseqs for all polypeptide chains,
              the second is the same for polynucleotides
    :raises: :class:`RuntimeError` if ent contains any residue which doesn't
             evaluate True for `r.IsPeptideLinking() or r.IsNucleotideLinking()`
             or when these two types are not strictly separated in different
             chains.
    """
    polypep_seqs = seq.CreateSequenceList()
    polynuc_seqs = seq.CreateSequenceList()

    for ch in ent.chains:
        n_res = len(ch.residues)
        n_pep = sum([r.IsPeptideLinking() for r in ch.residues])
        n_nuc = sum([r.IsNucleotideLinking() for r in ch.residues])
        
        # guarantee that we have either pep or nuc (no mix of the two)
        if n_pep > 0 and n_nuc > 0:
            raise RuntimeError(f"Must not mix peptide and nucleotide linking "
                               f"residues in same chain ({ch.GetName()})")

        if (n_pep + n_nuc) != n_res:
            raise RuntimeError("All residues must either be peptide_linking "
                               "or nucleotide_linking")

        s = ''.join([r.one_letter_code for r in ch.residues])
        if n_pep == n_res:
            polypep_seqs.AddSequence(seq.CreateSequence(ch.GetName(), s))
        elif n_nuc == n_res:
            polynuc_seqs.AddSequence(seq.CreateSequence(ch.GetName(), s))
        else:
            raise RuntimeError("This shouldnt happen")

    return (polypep_seqs, polynuc_seqs)

def _GroupSequences(seqs, seqid_thr, gap_thr, aligner, chem_type):
    """Get list of alignments representing groups of equivalent sequences

    :param seqid_thr: Threshold used to decide when two chains are identical.
                      95 percent tolerates the few mutations crystallographers
                      like to do.
    :type seqid_thr:  :class:`float`
    :param gap_thr: Additional threshold to avoid gappy alignments with high
                    seqid. The reason for not just normalizing seqid by the
                    longer sequence is that one sequence might be a perfect
                    subsequence of the other but only cover half of it. This
                    threshold checks for a maximum allowed fraction of gaps
                    in any of the two sequences after stripping terminal gaps.
    :type gap_thr: :class:`float`
    :param aligner: Helper class to generate pairwise alignments
    :type aligner: :class:`_Aligner`
    :param chem_type: ChemType of seqs which is passed to *aligner*, must be in
                      [:class:`ost.mol.ChemType.AMINOACIDS`,
                      :class:`ost.mol.ChemType.NUCLEOTIDES`]
    :type chem_type: :class:`ost.mol.ChemType` 
    :returns: A list of alignments, one alignment for each group
              with longest sequence (reference) as first sequence.
    :rtype: :class:`ost.seq.AlignmentList`
    """
    sim_seqs = list()
    aln_cache = dict() # cache alignments to generate an MSA for each group
                       # at the very end
    for i, j in itertools.combinations(range(len(seqs)), 2):
        aln  = aligner.Align(seqs[i], seqs[j], chem_type)
        seqid, gap_frac_i, gap_frac_j = _GetAlnProps(aln)
        if seqid >= seqid_thr and gap_frac_i < gap_thr and gap_frac_j < gap_thr:
            sim_seqs.append((i, j))
        aln_cache[(i,j)] = aln

    # trivial case: no matching pairs
    if len(sim_seqs) == 0:
        aln_list = seq.AlignmentList()
        for s in seqs:
            aln_list.append(seq.CreateAlignment(s))
        return aln_list

    # merge transitive pairs
    groups = list()
    for i, j in sim_seqs:
        found = False
        for g in groups:
            if i in g or j in g:
                found = True
                g.add(i)
                g.add(j)
        if not found:
            groups.append(set([i, j]))

    # sort based on sequence length
    sorted_groups = list()
    for g in groups:
        tmp = sorted([[len(seqs[i]), i] for i in g], reverse=True)
        sorted_groups.append([x[1] for x in tmp])

    # add all single chains
    for i in range(len(seqs)):
        if not any(i in g for g in sorted_groups):
            sorted_groups.append([i])

    # translate from indices back to sequences and directly generate alignments
    # of the groups with the longest (first) sequence as reference
    aln_list = seq.AlignmentList()
    for g in sorted_groups:
        if len(g) == 1:
            # aln with one single sequence
            aln_list.append(seq.CreateAlignment(seqs[g[0]]))
        else:
            # obtain pairwise aln of first sequence (reference) to all others
            alns = seq.AlignmentList()
            i = g[0]
            for j in g[1:]:
                if (i,j) in aln_cache:
                    alns.append(aln_cache[(i,j)])
                else:
                    # need to revert
                    aln = aln_cache[(j,i)]
                    alns.append(seq.CreateAlignment(aln.GetSequence(1),
                                                    aln.GetSequence(0)))
            # and merge
            aln_list.append(seq.alg.MergePairwiseAlignments(alns, seqs[i]))

    # transfer attached views
    seq_dict = {s.GetName(): s for s in seqs}
    for aln_idx in range(len(aln_list)):
        for aln_s_idx in range(aln_list[aln_idx].GetCount()):
            s_name = aln_list[aln_idx].GetSequence(aln_s_idx).GetName()
            s = seq_dict[s_name]
            if s.HasAttachedView():
                aln_list[aln_idx].AttachView(aln_s_idx, s.GetAttachedView())

    return aln_list

def _MapSequence(ref_seqs, ref_types, s, s_type, seqid_thr, gap_thr, aligner):
    """Tries top map *s* onto any of the sequences in *ref_seqs* of equal type

    :param ref_seqs: Reference sequences 
    :type ref_seqs: :class:`ost.seq.SequenceList`
    :param ref_types: Types of reference sequences, e.g.
                      ost.mol.ChemType.AminoAcids
    :type ref_types: :class:`list` of :class:`ost.mol.ChemType`
    :param s: Sequence to map
    :type s: :class:`ost.seq.SequenceHandle`
    :param s_type: Type of *s*, only try mapping to sequences in *ref_seqs*
                   with equal type as defined in *ref_types*
    :param seqid_thr: Threshold used to decide whether sequence is mapped.
    :type seqid_thr:  :class:`float`
    :param gap_thr: Additional threshold to avoid gappy alignments with high
                    seqid. The reason for not just normalizing seqid by the
                    longer sequence is that one sequence might be a perfect
                    subsequence of the other but only cover half of it. This
                    threshold checks for a maximum allowed fraction of gaps
                    in any of the two sequences after stripping terminal gaps.
    :type gap_thr: :class:`float`
    :param aligner: Helper class to generate pairwise alignments
    :type aligner: :class:`_Aligner`
    :returns: Tuple with two elements. 1) index of sequence in *ref_seqs* to
              which *s* can be mapped 2) Pairwise sequence alignment with 
              sequence from *ref_seqs* as first sequence. Both elements are
              None if no mapping can be found.
    :raises: :class:`RuntimeError` if mapping is ambiguous, i.e. *s*
             successfully maps to more than one sequence in *ref_seqs* 
    """
    map_idx = None
    map_aln = None
    for ref_idx, ref_seq in enumerate(ref_seqs):
        aln = aligner.Align(ref_seq, s, s_type)
        seqid, gap_frac_i, gap_frac_j = _GetAlnProps(aln)
        if seqid >= seqid_thr and gap_frac_i < gap_thr and gap_frac_j < gap_thr:
            if map_idx is not None:
                # match is ambiguous!
                raise RuntimeError(f"Ambiguous mapping for mdl chain "
                                   f"{s.GetName()}. Maps to chemically "
                                   f"equivalent groups with ref sequences "
                                   f"{ref_seqs[map_idx].GetName()} and "
                                   f"{ref_seqs[ref_idx].GetName()}.")
            map_idx = ref_idx
            map_aln = aln
    return (map_idx, map_aln)

def _GetRefMdlAlns(ref_chem_groups, ref_chem_group_msas, mdl_chem_groups,
                   mdl_chem_group_alns):
    """ Get all possible ref/mdl chain alignments given chem group mapping

    :param ref_chem_groups: :attr:`ChainMapper.chem_groups`
    :type ref_chem_groups: :class:`list` of :class:`list` of :class:`str`
    :param ref_chem_group_msas: :attr:`ChainMapper.chem_group_alignments`
    :type ref_chem_group_msas: :class:`ost.seq.AlignmentList`
    :param mdl_chem_groups: Groups of model chains that are mapped to
                            *ref_chem_groups*. Return value of
                            :func:`ChainMapper.GetChemMapping`.
    :type mdl_chem_groups: :class:`list` of :class:`list` of :class:`str`
    :param mdl_chem_group_alns: A pairwise sequence alignment for every chain
                                in *mdl_chem_groups* that aligns these sequences
                                to the respective reference sequence.
                                Return values of
                                :func:`ChainMapper.GetChemMapping`.
    :type mdl_chem_group_alns: :class:`list` of :class:`ost.seq.AlignmentList`
    :returns: A dictionary holding all possible ref/mdl chain alignments. Keys
              in that dictionary are tuples of the form (ref_ch, mdl_ch) and
              values are the respective pairwise alignments with first sequence
              being from ref, the second from mdl.
    """
    # alignment of each model chain to chem_group reference sequence
    mdl_alns = dict()
    for alns in mdl_chem_group_alns:
        for aln in alns:
            mdl_chain_name = aln.GetSequence(1).GetName()
            mdl_alns[mdl_chain_name] = aln

    # generate all alignments between ref/mdl chain atomseqs that we will
    # ever observe
    ref_mdl_alns = dict()
    for ref_chains, mdl_chains, ref_aln in zip(ref_chem_groups, mdl_chem_groups,
                                               ref_chem_group_msas):
        for ref_ch in ref_chains:
            for mdl_ch in mdl_chains:
                # obtain alignments of mdl and ref chains towards chem
                # group ref sequence and merge them
                aln_list = seq.AlignmentList()
                # do ref aln
                s1 = ref_aln.GetSequence(0)
                s2 = ref_aln.GetSequence(ref_chains.index(ref_ch))
                aln_list.append(seq.CreateAlignment(s1, s2))
                # do mdl aln
                aln_list.append(mdl_alns[mdl_ch])
                # merge
                ref_seq = seq.CreateSequence(s1.GetName(),
                                             s1.GetGaplessString())
                merged_aln = seq.alg.MergePairwiseAlignments(aln_list,
                                                             ref_seq)
                # merged_aln:
                # seq1: ref seq of chem group
                # seq2: seq of ref chain
                # seq3: seq of mdl chain
                # => we need the alignment between seq2 and seq3
                aln = seq.CreateAlignment(merged_aln.GetSequence(1),
                                          merged_aln.GetSequence(2))
                ref_mdl_alns[(ref_ch, mdl_ch)] = aln

    return ref_mdl_alns

class _Aligner:
    def __init__(self, pep_subst_mat = seq.alg.BLOSUM100, pep_gap_open = -5,
                 pep_gap_ext = -2, nuc_subst_mat = seq.alg.NUC44,
                 nuc_gap_open = -4, nuc_gap_ext = -4, resnum_aln = False):
        """ Helper class to compute alignments

        Sets default values for substitution matrix, gap open and gap extension
        penalties. They are only used in default mode (Needleman-Wunsch aln).
        If *resnum_aln* is True, only residue numbers of views that are attached
        to input sequences are considered. 
        """
        self.pep_subst_mat = pep_subst_mat
        self.pep_gap_open = pep_gap_open
        self.pep_gap_ext = pep_gap_ext
        self.nuc_subst_mat = nuc_subst_mat
        self.nuc_gap_open = nuc_gap_open
        self.nuc_gap_ext = nuc_gap_ext
        self.resnum_aln = resnum_aln

    def Align(self, s1, s2, chem_type=None):
        if self.resnum_aln:
            return self.ResNumAlign(s1, s2)
        else:
            if chem_type is None:
                raise RuntimeError("Must specify chem_type for NW alignment")
            return self.NWAlign(s1, s2, chem_type) 


    def NWAlign(self, s1, s2, chem_type):
        """ Returns pairwise alignment using Needleman-Wunsch algorithm
    
        :param s1: First sequence to align
        :type s1: :class:`ost.seq.SequenceHandle`
        :param s2: Second sequence to align
        :type s2: :class:`ost.seq.SequenceHandle`
        :param chem_type: Must be in [:class:`ost.mol.ChemType.AMINOACIDS`,
                          :class:`ost.mol.ChemType.NUCLEOTIDES`], determines
                          substitution matrix and gap open/extension penalties
        :type chem_type: :class:`ost.mol.ChemType`
        :returns: Alignment with s1 as first and s2 as second sequence 
        """
        if chem_type == mol.ChemType.AMINOACIDS:
            return seq.alg.GlobalAlign(s1, s2, self.pep_subst_mat,
                                       gap_open=self.pep_gap_open,
                                       gap_ext=self.pep_gap_ext)[0]
        elif chem_type == mol.ChemType.NUCLEOTIDES:
            return seq.alg.GlobalAlign(s1, s2, self.nuc_subst_mat,
                                       gap_open=self.nuc_gap_open,
                                       gap_ext=self.nuc_gap_ext)[0]
        else:
            raise RuntimeError("Invalid ChemType")
        return aln

    def ResNumAlign(self, s1, s2):
        """ Returns pairwise alignment using residue numbers of attached views
    
        :param s1: First sequence to align, must have :class:`ost.mol.EntityView`
                   attached
        :type s1: :class:`ost.seq.SequenceHandle`
        :param s2: Second sequence to align, must have :class:`ost.mol.EntityView`
                   attached
        :type s2: :class:`ost.seq.SequenceHandle`
        """
        assert(s1.HasAttachedView())
        assert(s2.HasAttachedView())
        v1 = s1.GetAttachedView()
        rnums1 = [r.GetNumber().GetNum() for r in v1.residues]
        v2 = s2.GetAttachedView()
        rnums2 = [r.GetNumber().GetNum() for r in v2.residues]

        min_num = min(rnums1[0], rnums2[0])
        max_num = max(rnums1[-1], rnums2[-1])
        aln_length = max_num - min_num + 1

        aln_s1 = ['-'] * aln_length
        for r, rnum in zip(v1.residues, rnums1):
            aln_s1[rnum-min_num] = r.one_letter_code

        aln_s2 = ['-'] * aln_length
        for r, rnum in zip(v2.residues, rnums2):
            aln_s2[rnum-min_num] = r.one_letter_code

        aln = seq.CreateAlignment()
        aln.AddSequence(seq.CreateSequence(s1.GetName(), ''.join(aln_s1)))
        aln.AddSequence(seq.CreateSequence(s2.GetName(), ''.join(aln_s2)))
        return aln

def _GetAlnProps(aln):
    """Returns basic properties of *aln*

    :param aln: Alignment to compute properties
    :type aln: :class:`seq.AlignmentHandle`
    :returns: Tuple with 3 elements. 1) sequence identify in range [0, 100] 
              considering aligned columns 2) Fraction of gaps between
              first and last aligned column in s1 4) same for s2. 
    """
    assert(aln.GetCount() == 2)
    seqid = seq.alg.SequenceIdentity(aln)
    n_gaps_1 = str(aln.GetSequence(0)).strip('-').count('-')
    n_gaps_2 = str(aln.GetSequence(1)).strip('-').count('-')
    gap_frac_1 = float(n_gaps_1)/len(aln.GetSequence(0).GetGaplessString())
    gap_frac_2 = float(n_gaps_2)/len(aln.GetSequence(1).GetGaplessString())
    return (seqid, gap_frac_1, gap_frac_2)    

def _NChemGroupMappings(ref_chains, mdl_chains):
    """ Number of mappings within one chem group

    :param ref_chains: Reference chains
    :type ref_chains: :class:`list` of :class:`str`
    :param mdl_chains: Model chains that are mapped onto *ref_chains*
    :type mdl_chains: :class:`list` of :class:`str`
    :returns: Number of possible mappings of *mdl_chains* onto *ref_chains*
    """
    n_ref = len(ref_chains)
    n_mdl = len(mdl_chains)
    if n_ref == n_mdl:
        return factorial(n_ref)
    elif n_ref > n_mdl:
        n_choose_k = binom(n_ref, n_mdl)
        return n_choose_k * factorial(n_mdl)
    else:
        n_choose_k = binom(n_mdl, n_ref)
        return n_choose_k * factorial(n_ref)

def _NMappings(ref_chains, mdl_chains):
    """ Number of mappings for a full chem mapping

    :param ref_chains: Chem groups of reference
    :type ref_chains: :class:`list` of :class:`list` of :class:`str`
    :param mdl_chains: Model chains that map onto those chem groups
    :type mdl_chains: :class:`list` of :class:`list` of :class:`str`
    :returns: Number of possible mappings of *mdl_chains* onto *ref_chains*
    """
    assert(len(ref_chains) == len(mdl_chains))
    n = 1
    for a,b in zip(ref_chains, mdl_chains):
        n *= _NChemGroupMappings(a,b)
    return n

def _NMappingsWithin(ref_chains, mdl_chains, max_mappings):
    """ Check whether total number of mappings is smaller than given maximum

    In principle the same as :func:`_NMappings` but it stops as soon as the
    maximum is hit.

    :param ref_chains: Chem groups of reference
    :type ref_chains: :class:`list` of :class:`list` of :class:`str`
    :param mdl_chains: Model chains that map onto those chem groups
    :type mdl_chains: :class:`list` of :class:`list` of :class:`str`
    :param max_mappings: Number of max allowed mappings
    :returns: Whether number of possible mappings of *mdl_chains* onto
              *ref_chains* is below or equal *max_mappings*.
    """
    assert(len(ref_chains) == len(mdl_chains))
    n = 1
    for a,b in zip(ref_chains, mdl_chains):
        n *= _NChemGroupMappings(a,b)
        if n > max_mappings:
            return False
    return True

def _RefSmallerGenerator(ref_chains, mdl_chains):
    """ Returns all possible ways to map mdl_chains onto ref_chains

    Specific for the case where len(ref_chains) < len(mdl_chains)
    """
    for c in itertools.combinations(mdl_chains, len(ref_chains)):
        for p in itertools.permutations(c):
            yield p

def _RefLargerGenerator(ref_chains, mdl_chains):
    """ Returns all possible ways to map mdl_chains onto ref_chains

    Specific for the case where len(ref_chains) > len(mdl_chains)
    Ref chains without mapped mdl chain are assigned None
    """
    n_ref = len(ref_chains)
    n_mdl = len(mdl_chains)
    for c in itertools.combinations(range(n_ref), n_mdl):
        for p in itertools.permutations(mdl_chains):
            ret_list = [None] * n_ref
            for idx, ch in zip(c, p):
                ret_list[idx] = ch
            yield ret_list

def _ChainMappings(ref_chains, mdl_chains):
    """Returns all possible ways to map *mdl_chains* onto fixed *ref_chains*

    :param ref_chains: List of list of chemically equivalent chains in reference
    :type ref_chains: :class:`list` of :class:`list`
    :param mdl_chains: Equally long list of list of chemically equivalent chains
                       in model that map on those ref chains.
    :type mdl_chains: :class:`list` of :class:`list`
    :returns: Iterator over all possible mappings of *mdl_chains* onto fixed
              *ref_chains*. Potentially contains None as padding when number of
              model chains for a certain mapping is smaller than the according
              reference chains.
              Example: _ChainMappings([['A', 'B', 'C'], ['D', 'E']],
                                      [['x', 'y'], ['i', 'j']])
              gives an iterator over: [(['x', 'y', None], ('i', 'j')),
                                       (['x', 'y', None], ('j', 'i')),
                                       (['y', 'x', None], ('i', 'j')),
                                       (['y', 'x', None], ('j', 'i')),
                                       (['x', None, 'y'], ('i', 'j')),
                                       (['x', None, 'y'], ('j', 'i')),
                                       (['y', None, 'x'], ('i', 'j')),
                                       (['y', None, 'x'], ('j', 'i')),
                                       ([None, 'x', 'y'], ('i', 'j')),
                                       ([None, 'x', 'y'], ('j', 'i')),
                                       ([None, 'y', 'x'], ('i', 'j')),
                                       ([None, 'y', 'x'], ('j', 'i'))]
    """
    assert(len(ref_chains) == len(mdl_chains))
    # one iterator per mapping representing all mdl combinations relative to
    # reference
    iterators = list()
    for ref, mdl in zip(ref_chains, mdl_chains):
        if len(ref) == len(mdl):
            iterators.append(itertools.permutations(mdl))
        elif len(ref) < len(mdl):
            iterators.append(_RefSmallerGenerator(ref, mdl))
        else:
            iterators.append(_RefLargerGenerator(ref, mdl))

    return itertools.product(*iterators)


class _lDDTDecomposer:

    def __init__(self, ref, mdl, ref_mdl_alns, inclusion_radius = 15.0,
                 thresholds = [0.5, 1.0, 2.0, 4.0]):

        self.ref = ref
        self.mdl = mdl
        self.ref_mdl_alns = ref_mdl_alns
        self.inclusion_radius = inclusion_radius
        self.thresholds = thresholds

        # keep track of single chains and interfaces in ref
        self.ref_chains = list() # e.g. ['A', 'B', 'C']
        self.ref_interfaces = list() # e.g. [('A', 'B'), ('A', 'C')]

        # holds lDDT scorer for each chain in ref
        # key: chain name, value: scorer
        self.single_chain_scorer = dict()

        # cache for single chain conserved contacts
        # key: tuple (ref_ch, mdl_ch) value: number of conserved contacts
        self.single_chain_cache = dict()

        # holds lDDT scorer for each pairwise interface in target
        # key: tuple (ref_ch1, ref_ch2), value: scorer
        self.interface_scorer = dict()

        # cache for interface conserved contacts
        # key: tuple of tuple ((ref_ch1, ref_ch2),((mdl_ch1, mdl_ch2))
        # value: number of conserved contacts
        self.interface_cache = dict()

        self.n = 0

        self._SetupScorer()

    def _SetupScorer(self):
        for ch in self.ref.chains:
            # Select everything close to that chain
            query = f"{self.inclusion_radius} <> [cname={ch.GetName()}] "
            query += f"and cname!={ch.GetName()}"
            for close_ch in self.ref.Select(query).chains:
                k1 = (ch.GetName(), close_ch.GetName())
                k2 = (close_ch.GetName(), ch.GetName())
                if k1 not in self.interface_scorer and \
                k2 not in self.interface_scorer:
                    dimer_ref = self.ref.Select(f"cname={k1[0]},{k1[1]}")
                    s = lddt.lDDTScorer(dimer_ref, bb_only=True)
                    self.interface_scorer[k1] = s
                    self.interface_scorer[k2] = s
                    self.n += self.interface_scorer[k1].n_distances_ic
                    self.ref_interfaces.append(k1)
                    # single chain scorer are actually interface scorers to save
                    # some distance calculations
                    if ch.GetName() not in self.single_chain_scorer:
                        self.single_chain_scorer[ch.GetName()] = s
                        self.n += s.GetNChainContacts(ch.GetName(),
                                                      no_interchain=True)
                        self.ref_chains.append(ch.GetName())
                    if close_ch.GetName() not in self.single_chain_scorer:
                        self.single_chain_scorer[close_ch.GetName()] = s
                        self.n += s.GetNChainContacts(close_ch.GetName(),
                                                      no_interchain=True)
                        self.ref_chains.append(close_ch.GetName())

        # add any missing single chain scorer
        for ch in self.ref.chains:
            if ch.GetName() not in self.single_chain_scorer:
                single_chain_ref = self.ref.Select(f"cname={ch.GetName()}")
                self.single_chain_scorer[ch.GetName()] = \
                lddt.lDDTScorer(single_chain_ref, bb_only = True)
                self.n += self.single_chain_scorer[ch.GetName()].n_distances
                self.ref_chains.append(ch.GetName())

    def lDDT(self, ref_chain_groups, mdl_chain_groups):

        flat_map = dict()
        for ref_chains, mdl_chains in zip(ref_chain_groups, mdl_chain_groups):
            for ref_ch, mdl_ch in zip(ref_chains, mdl_chains):
                flat_map[ref_ch] = mdl_ch

        return self.lDDTFromFlatMap(flat_map)


    def lDDTFromFlatMap(self, flat_map):
        conserved = 0

        # do single chain scores
        for ref_ch in self.ref_chains:
            if ref_ch in flat_map and flat_map[ref_ch] is not None:
                conserved += self.SCCounts(ref_ch, flat_map[ref_ch])

        # do interfaces
        for ref_ch1, ref_ch2 in self.ref_interfaces:
            if ref_ch1 in flat_map and ref_ch2 in flat_map:
                mdl_ch1 = flat_map[ref_ch1]
                mdl_ch2 = flat_map[ref_ch2]
                if mdl_ch1 is not None and mdl_ch2 is not None:
                    conserved += self.IntCounts(ref_ch1, ref_ch2, mdl_ch1,
                                                mdl_ch2)

        return conserved / (len(self.thresholds) * self.n)

    def SCCounts(self, ref_ch, mdl_ch):
        if not (ref_ch, mdl_ch) in self.single_chain_cache:
            alns = dict()
            alns[mdl_ch] = self.ref_mdl_alns[(ref_ch, mdl_ch)]
            mdl_sel = self.mdl.Select(f"cname={mdl_ch}")
            s = self.single_chain_scorer[ref_ch]
            _,_,_,conserved,_,_,_ = s.lDDT(mdl_sel,
                                           residue_mapping=alns,
                                           return_dist_test=True,
                                           no_interchain=True,
                                           chain_mapping={mdl_ch: ref_ch},
                                           check_resnames=False)
            self.single_chain_cache[(ref_ch, mdl_ch)] = conserved
        return self.single_chain_cache[(ref_ch, mdl_ch)]

    def IntCounts(self, ref_ch1, ref_ch2, mdl_ch1, mdl_ch2):
        k1 = ((ref_ch1, ref_ch2),(mdl_ch1, mdl_ch2))
        k2 = ((ref_ch2, ref_ch1),(mdl_ch2, mdl_ch1))
        if k1 not in self.interface_cache and k2 not in self.interface_cache:
            alns = dict()
            alns[mdl_ch1] = self.ref_mdl_alns[(ref_ch1, mdl_ch1)]
            alns[mdl_ch2] = self.ref_mdl_alns[(ref_ch2, mdl_ch2)]
            mdl_sel = self.mdl.Select(f"cname={mdl_ch1},{mdl_ch2}")
            s = self.interface_scorer[(ref_ch1, ref_ch2)]
            _,_,_,conserved,_,_,_ = s.lDDT(mdl_sel,
                                           residue_mapping=alns,
                                           return_dist_test=True,
                                           no_intrachain=True,
                                           chain_mapping={mdl_ch1: ref_ch1,
                                                          mdl_ch2: ref_ch2},
                                           check_resnames=False)
            self.interface_cache[k1] = conserved
            self.interface_cache[k2] = conserved
        return self.interface_cache[k1]


class _GreedySearcher(_lDDTDecomposer):
    def __init__(self, ref, mdl, ref_chem_groups, mdl_chem_groups,
                 ref_mdl_alns, inclusion_radius = 15.0,
                 thresholds = [0.5, 1.0, 2.0, 4.0],
                 steep_opt_rate = None):
        super().__init__(ref, mdl, ref_mdl_alns,
                         inclusion_radius = inclusion_radius,
                         thresholds = thresholds)
        self.steep_opt_rate = steep_opt_rate
        self.neighbors = {k: set() for k in self.ref_chains}
        for k in self.interface_scorer.keys():
            self.neighbors[k[0]].add(k[1])
            self.neighbors[k[1]].add(k[0])

        assert(len(ref_chem_groups) == len(mdl_chem_groups))
        self.ref_chem_groups = ref_chem_groups
        self.mdl_chem_groups = mdl_chem_groups
        self.ref_ch_group_mapper = dict()
        self.mdl_ch_group_mapper = dict()
        for g_idx, (ref_g, mdl_g) in enumerate(zip(ref_chem_groups,
                                                   mdl_chem_groups)):
            for ch in ref_g:
                self.ref_ch_group_mapper[ch] = g_idx
            for ch in mdl_g:
                self.mdl_ch_group_mapper[ch] = g_idx

        # keep track of mdl chains that potentially give lDDT contributions,
        # i.e. they have locations within inclusion_radius + max(thresholds)
        self.mdl_neighbors = dict()
        d = self.inclusion_radius + max(self.thresholds)
        for ch in self.mdl.chains:
            ch_name = ch.GetName()
            self.mdl_neighbors[ch_name] = set()
            query = f"{d} <> [cname={ch_name}] and cname !={ch_name}"
            for close_ch in self.mdl.Select(query).chains:
                self.mdl_neighbors[ch_name].add(close_ch.GetName())


    def ExtendMapping(self, mapping, max_ext = None):

        if len(mapping) == 0:
            raise RuntimError("Mapping must contain a starting point")

        for ref_ch, mdl_ch in mapping.items():
            assert(ref_ch in self.ref_ch_group_mapper)
            assert(mdl_ch in self.mdl_ch_group_mapper)
            assert(self.ref_ch_group_mapper[ref_ch] == \
                   self.mdl_ch_group_mapper[mdl_ch])

        # Ref chains onto which we can map. The algorithm starts with a mapping
        # on ref_ch. From there we can start to expand to connected neighbors.
        # All neighbors that we can reach from the already mapped chains are
        # stored in this set which will be updated during runtime
        map_targets = set()
        for ref_ch in mapping.keys():
            map_targets.update(self.neighbors[ref_ch])

        # remove the already mapped chains
        for ref_ch in mapping.keys():
            map_targets.discard(ref_ch)

        if len(map_targets) == 0:
            return mapping # nothing to extend

        # keep track of what model chains are not yet mapped for each chem group
        free_mdl_chains = list()
        for chem_group in self.mdl_chem_groups:
            tmp = [x for x in chem_group if x not in mapping.values()]
            free_mdl_chains.append(set(tmp))

        # keep track of what ref chains got a mapping
        newly_mapped_ref_chains = list()

        something_happened = True
        while something_happened:
            something_happened=False

            if self.steep_opt_rate is not None:
                n_chains = len(newly_mapped_ref_chains)
                if n_chains > 0 and n_chains % self.steep_opt_rate == 0:
                    mapping = self._SteepOpt(mapping, newly_mapped_ref_chains)

            if max_ext is not None and len(newly_mapped_ref_chains) >= max_ext:
                break

            max_n = 0
            max_mapping = None
            for ref_ch in map_targets:
                chem_group_idx = self.ref_ch_group_mapper[ref_ch]
                for mdl_ch in free_mdl_chains[chem_group_idx]:
                    # single chain score
                    n_single = self.SCCounts(ref_ch, mdl_ch)
                    # scores towards neighbors that are already mapped
                    n_inter = 0
                    for neighbor in self.neighbors[ref_ch]:
                        if neighbor in mapping and mapping[neighbor] in \
                        self.mdl_neighbors[mdl_ch]:
                            n_inter += self.IntCounts(ref_ch, neighbor, mdl_ch,
                                                      mapping[neighbor])
                    n = n_single + n_inter

                    if n_inter > 0 and n > max_n:
                        # Only accept a new solution if its actually connected
                        # i.e. n_inter > 0. Otherwise we could just map a big
                        # fat mdl chain sitting somewhere in Nirvana
                        max_n = n
                        max_mapping = (ref_ch, mdl_ch)
     
            if max_n > 0:
                something_happened = True
                # assign new found mapping
                mapping[max_mapping[0]] = max_mapping[1]

                # add all neighboring chains to map targets as they are now
                # reachable
                for neighbor in self.neighbors[max_mapping[0]]:
                    if neighbor not in mapping:
                        map_targets.add(neighbor)

                # remove the ref chain from map targets
                map_targets.remove(max_mapping[0])

                # remove the mdl chain from free_mdl_chains - its taken...
                chem_group_idx = self.ref_ch_group_mapper[max_mapping[0]]
                free_mdl_chains[chem_group_idx].remove(max_mapping[1])

                # keep track of what ref chains got a mapping
                newly_mapped_ref_chains.append(max_mapping[0])

        return mapping


    def _SteepOpt(self, mapping, chains_to_optimize=None):

        # just optimize ALL ref chains if nothing specified
        if chains_to_optimize is None:
            chains_to_optimize = mapping.keys()

        # make sure that we only have ref chains which are actually mapped
        ref_chains = [x for x in chains_to_optimize if mapping[x] is not None]

        # group ref chains to be optimized into chem groups
        tmp = dict()
        for ch in ref_chains:
            chem_group_idx = self.ref_ch_group_mapper[ch] 
            if chem_group_idx in tmp:
                tmp[chem_group_idx].append(ch)
            else:
                tmp[chem_group_idx] = [ch]
        chem_groups = list(tmp.values())

        # try all possible mapping swaps. Swaps that improve the score are
        # immediately accepted and we start all over again
        current_lddt = self.lDDTFromFlatMap(mapping)
        something_happened = True
        while something_happened:
            something_happened = False
            for chem_group in chem_groups:
                if something_happened:
                    break
                for ch1, ch2 in itertools.combinations(chem_group, 2):
                    swapped_mapping = dict(mapping)
                    swapped_mapping[ch1] = mapping[ch2]
                    swapped_mapping[ch2] = mapping[ch1]
                    score = self.lDDTFromFlatMap(swapped_mapping)
                    if score > current_lddt:
                        something_happened = True
                        mapping = swapped_mapping
                        current_lddt = score
                        break        

        return mapping
