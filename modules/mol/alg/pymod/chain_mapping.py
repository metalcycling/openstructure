import itertools

from ost import seq
from ost import mol

from ost.mol.alg import lddt

class ChainMapper:
    def __init__(self, target):
        """Object to compute chain mappings

        :param target: Target structure onto which models are mapped.
                       Computations happen on a selection only containing
                       polypeptides and polynucleotides.
        :type target: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        """
        self.target = target.Select("peptide=true or nucleotide=true")

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
            _GetChemGroupAlignments(self.target)
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
                self._chem_group_ref_seqs.AddSequence(s)
        return self._chem_group_ref_seqs

    @property
    def chem_group_types(self):
        """ChemType of each group in :attr:`~chem_groups`

        Specifying if groups are poly-peptides/nucleotides, i.e. 
        :class:`ost.mol.ChemType.AMINOACIDS` or
        :class:`ost.mol.ChemType.AMINOACIDS` 
        
        :getter: Computed on first use (cached)
        :type: :class:`list` of :class:`ost.mol.ChemType`
        """
        if self._chem_group_types is None:
            self._chem_group_alignments, self._chem_group_types = \
            _GetChemGroupAlignments(self.target)            
        return self._chem_group_types
        
    def GetChemMapping(self, model):
        """Maps sequences in *model* to chem_groups of target

        :param model: Model from which to extract sequences, a
                      selection that only includes peptides and 
                      nucleotides is performed.
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :returns: Tuple with two lists of length `len(self.chem_groups)`:
                  1) Each element is a :class:`list` with mdl chain names that
                  map to the chem group at that position.
                  2) Each element is a :class:`ost.seq.AlignmentList` aligning
                  these mdl chain sequences to the chem group ref sequences.
        """
        mdl = model.Select("peptide=true or nucleotide=true")
        mdl_pep_seqs, mdl_nuc_seqs = _GetAtomSeqs(mdl)

        mapping = [list() for x in self.chem_groups]
        alns = [seq.AlignmentList() for x in self.chem_groups]

        for s in mdl_pep_seqs:
            idx, aln = _MapSequence(self.chem_group_ref_seqs, 
                                    self.chem_group_types,
                                    s, mol.ChemType.AMINOACIDS,
                                    95., 0.1, seq.alg.BLOSUM100,
                                    -5, -2)
            if idx is not None:
                mapping[idx].append(s.GetName())
                alns[idx].append(aln)

        for s in mdl_nuc_seqs:
            idx, aln = _MapSequence(self.chem_group_ref_seqs, 
                                    self.chem_group_types,
                                    s, mol.ChemType.NUCLEOTIDES,
                                    95., 0.1, seq.alg.NUC44,
                                    -4, -4)
            if idx is not None:
                mapping[idx].append(s.GetName())
                alns[idx].append(aln)

        return (mapping, alns)


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

        mdl = model.Select("peptide=true or nucleotide=true")
        chem_mapping, chem_group_alns = self.GetChemMapping(mdl)

        # check for the simplest case
        only_one_to_one = True
        for ref_chains, mdl_chains in zip(self.chem_groups, chem_mapping):
            if len(ref_chains) != 1 or len(mdl_chains) not in [0, 1]:
                only_one_to_one = False
                break
        if only_one_to_one:
            # skip ref chem groups with no mapped mdl chain
            return [(a,b) for a,b in zip(self.chem_groups, chem_mapping) if len(b) == 1]

        # alignment of each model chain to chem_group reference sequence
        mdl_alns = dict()
        for alns in chem_group_alns:
            for aln in alns:
                mdl_chain_name = aln.GetSequence(1).GetName()
                mdl_alns[mdl_chain_name] = aln

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

                        # obtain alignments of mdl and ref chains towards chem
                        # group ref sequence and merge them
                        aln_list = seq.AlignmentList()
                        # do ref aln
                        s1 = ref_aln.GetSequence(0)
                        s2 = ref_aln.GetSequence(ref_chem_group.index(ref_ch))
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
                        lddt_alns[mdl_ch] = aln

            lDDT, _ = lddt_scorer.lDDT(mdl, thresholds=thresholds,
                                       chain_mapping=lddt_chain_mapping,
                                       residue_mapping = lddt_alns,
                                       check_resnames = False)
            if lDDT > best_lddt:
                best_mapping = mapping
                best_lddt = lDDT

        return best_mapping

def _GetChemGroupAlignments(ent, seqid_thr=95., gap_thr=0.1):
    """Returns alignments with groups of chemically equivalent chains

    :param ent: Entity to process
    :type ent:  :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
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
    :returns: Tuple with first element being an AlignmentList. Each alignment
              represents a group of chemically equivalent chains and the first
              sequence is the longest. Second element is a list of equivalent
              length specifying the types of the groups. List elements are in
              [:class:`ost.ChemType.AMINOACIDS`,
              :class:`ost.ChemType.NUCLEOTIDES`] 
    :raises: :class:`RuntimeError` if ent contains any residue which doesn't
             evaluate True for `r.IsPeptideLinking() or r.IsNucleotideLinking()`
             or when these two types are not strictly separated to different
             chains.
    """
    pep_seqs, nuc_seqs = _GetAtomSeqs(ent)
    pep_groups = _GroupSequences(pep_seqs, seqid_thr, gap_thr,
                                 seq.alg.BLOSUM100, -5, -2)
    nuc_groups = _GroupSequences(nuc_seqs, seqid_thr, gap_thr,
                                 seq.alg.NUC44, -4, -4)
    group_types = [mol.ChemType.AMINOACIDS] * len(pep_groups)
    group_types += [mol.ChemType.NUCLEOTIDES] * len(nuc_groups)
    groups = pep_groups
    groups.extend(nuc_groups)
    return (groups, group_types)

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

def _GroupSequences(seqs, seqid_thr, gap_thr, subst_mat, gap_open, gap_ext):
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
    :param subst_mat: Substitution matrix for standard Needleman-Wunsch
    :type subst_mat: :class:`ost.seq.alg.SubstWeightMatrix`
    :param gap_open: Gap open penalty for standard Needleman-Wunsch
    :type gap_open: :class:`int`
    :param gap_ext: Gap extension penalty for standard Needleman-Wunsch
    :type gap_ext: :class:`int`

    :returns: A list of alignments, one alignment for each group
              with longest sequence (reference) as first sequence.
    :rtype: :class:`ost.seq.AlignmentList`
    """
    sim_seqs = list()
    aln_cache = dict() # cache alignments to generate an MSA for each group
                       # at the very end
    for i, j in itertools.combinations(range(len(seqs)), 2):
        aln, seqid, gap_frac_i, gap_frac_j = _Align(seqs[i], seqs[j], subst_mat,
                                                    gap_open, gap_ext)
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

    return aln_list

def _MapSequence(ref_seqs, ref_types, s, s_type, seqid_thr, gap_thr, subst_mat,
                 gap_open, gap_ext):
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
    :param subst_mat: Substitution matrix for standard Needleman-Wunsch
    :type subst_mat: :class:`ost.seq.alg.SubstWeightMatrix`
    :param gap_open: Gap open penalty for standard Needleman-Wunsch
    :type gap_open: :class:`int`
    :param gap_ext: Gap extension penalty for standard Needleman-Wunsch
    :type gap_ext: :class:`int`
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
        aln, seqid, gap_frac_i, gap_frac_j = _Align(ref_seq, s, subst_mat,
                                                    gap_open, gap_ext)
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

def _Align(s1, s2, subst_mat, gap_open, gap_ext):
    """Returns alignment and some distance metrics of that alignment

    :param s1: First sequence to align
    :type s1: :class:`ost.seq.SequenceHandle`
    :param s2: Second sequence to align
    :type s2: :class:`ost.seq.SequenceHandle`
    :param subst_mat: Substitution matrix for standard Needleman-Wunsch
    :type subst_mat: :class:`ost.seq.alg.SubstWeightMatrix`
    :param gap_open: Gap open penalty for standard Needleman-Wunsch
    :type gap_open: :class:`int`
    :param gap_open: Gap extension penalty for standard Needleman-Wunsch
    :type gap_open: :class:`int`
    :returns: Tuple with 4 elements. 1) aln with s1 as first and s2 as
              second sequence 2) sequence identify in range [0, 100] 
              considering aligned columns 3) Fraction of gaps between
              first and last aligned column in s1 4) same for s2.
    """
    aln = seq.alg.GlobalAlign(s1, s2, subst_mat, gap_open=gap_open,
                              gap_ext=gap_ext)[0]
    seqid = seq.alg.SequenceIdentity(aln)
    n_gaps_1 = str(aln.GetSequence(0)).strip('-').count('-')
    n_gaps_2 = str(aln.GetSequence(1)).strip('-').count('-')
    gap_frac_1 = float(n_gaps_1)/len(aln.GetSequence(0).GetGaplessString())
    gap_frac_2 = float(n_gaps_2)/len(aln.GetSequence(1).GetGaplessString())
    return (aln, seqid, gap_frac_1, gap_frac_2)

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
