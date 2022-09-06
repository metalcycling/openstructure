import itertools
import numpy as np
from scipy.spatial import distance

import time
from ost import mol

class QSEntity:
    """ Helper object for QS-score computation

    Holds structural information and getters for interacting chains, i.e.
    interfaces. Peptide residues are represented by their CB position
    (CA for GLY) and nucleotides by C3'.

    :param ent: Structure for QS score computation
    :type ent: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param contact_d: Pairwise distance of residues to be considered as contacts
    :type contact_d: :class:`float`
    """
    def __init__(self, ent, contact_d = 12.0):
        pep_query = "(peptide=true and (aname=\"CB\" or (rname=\"GLY\" and aname=\"CA\")))"
        nuc_query = "(nucleotide=True and aname=\"C3'\")"
        self._view = ent.Select(" or ".join([pep_query, nuc_query]))
        self._contact_d = contact_d

        # the following attributes will be lazily evaluated
        self._chain_names = None
        self._interacting_chains = None
        self._sequence = dict()
        self._pos = dict()
        self._pair_dist = dict()

    @property
    def view(self):
        """ Processed structure

        View that only contains representative atoms. That's CB for peptide
        residues (CA for GLY) and C3' for nucleotides.

        :type: :class:`ost.mol.EntityView`
        """
        return self._view

    @property
    def contact_d(self):
        """ Pairwise distance of residues to be considered as contacts

        Given at :class:`QSEntity` construction

        :type: :class:`float`
        """
        return self._contact_d

    @property
    def chain_names(self):
        """ Chain names in :attr:`~view`
 
        Names are sorted

        :type: :class:`list` of :class:`str`
        """
        if self._chain_names is None:
            self._chain_names = sorted([ch.name for ch in self.view.chains])
        return self._chain_names

    @property
    def interacting_chains(self):
        """ Pairs of chains in :attr:`~view` with at least one contact

        :type: :class:`list` of :class:`tuples`
        """
        if self._interacting_chains is None:
            self._interacting_chains = list()
            for x in itertools.combinations(self.chain_names, 2):
                if np.count_nonzero(self.PairDist(x[0], x[1]) < self.contact_d):
                    self._interacting_chains.append(x)
        return self._interacting_chains
    
    def GetChain(self, chain_name):
        """ Get chain by name

        :param chain_name: Chain in :attr:`~view`
        :type chain_name: :class:`str`
        """ 
        chain = self.view.FindChain(chain_name)
        if not chain.IsValid():
            raise RuntimeError(f"view has no chain named \"{chain_name}\"")
        return chain

    def GetSequence(self, chain_name):
        """ Get sequence of chain

        Returns sequence of specified chain as raw :class:`str`

        :param chain_name: Chain in :attr:`~view`
        :type chain_name: :class:`str`
        """
        if chain_name not in self._sequence:
            ch = self.GetChain(chain_name)
            s = ''.join([r.one_letter_code for r in ch.residues])
            self._sequence[chain_name] = s
        return self._sequence[chain_name]

    def GetPos(self, chain_name):
        """ Get representative positions of chain

        That's CB positions for peptide residues (CA for GLY) and C3' for
        nucleotides. Returns positions as a numpy array of shape
        (n_residues, 3).

        :param chain_name: Chain in :attr:`~view`
        :type chain_name: :class:`str`
        """
        if chain_name not in self._pos:
            ch = self.GetChain(chain_name)
            pos = np.zeros((len(ch.residues), 3))
            for i, r in enumerate(ch.residues):
                pos[i,:] = r.atoms[0].GetPos().data
            self._pos[chain_name] = pos
        return self._pos[chain_name]

    def PairDist(self, chain_name_one, chain_name_two):
        """ Get pairwise distances between two chains

        Returns distances as numpy array of shape (a, b).
        Where a is the number of residues of the chain that comes BEFORE the
        other in :attr:`~chain_names` 
        """
        key = (min(chain_name_one, chain_name_two),
               max(chain_name_one, chain_name_two))
        if key not in self._pair_dist:
            self._pair_dist[key] = distance.cdist(self.GetPos(key[0]),
                                                  self.GetPos(key[1]),
                                                  'euclidean')
        return self._pair_dist[key]

class QSScorer:
    """ Helper object to compute QS-score

    It's conceptually tightly integrated into the mechanisms from the 
    chain_mapping module. The prefered way to derive an object of type
    :class:`QSScorer` is through the static constructor:
    :func:`~FromMappingResult`. Example score computation including mapping:

    ::

        from ost.mol.alg.qsscore import QSScorer
        from ost.mol.alg.chain_mapping import ChainMapper

        ent_1 = io.LoadPDB("path_to_assembly_1.pdb")
        ent_2 = io.LoadPDB("path_to_assembly_2.pdb")

        chain_mapper = ChainMapper(ent_1)
        mapping_result = chain_mapper.GetlDDTMapping(ent_2)
        qs_scorer = QSScorer.FromMappingResult(mapping_result)
        print("score:", qs_scorer.GetQSScore(mapping_result.mapping))


    Formulas for QS scores:

    ::
  
      - QS_global = weighted_scores / (weight_sum + weight_extra_all)
      -> weighted_scores = sum(w(min(d1,d2)) * (1 - abs(d1-d2)/12)) for shared
      -> weight_sum = sum(w(min(d1,d2))) for shared
      -> weight_extra_all = sum(w(d)) for all non-shared
      -> w(d) = 1 if d <= 5, exp(-2 * ((d-5.0)/4.28)^2) else
  
    In the formulas above:

    * "d": CA/CB-CA/CB distance of an "inter-chain contact" ("d1", "d2" for
      "shared" contacts).
    * "mapped": we could map chains of two structures and align residues in
      :attr:`alignments`.
    * "shared": pairs of residues which are "mapped" and have
      "inter-chain contact" in both structures.
    * "inter-chain contact": CB-CB pairs (CA for GLY) with distance <= 12 A
      (fallback to CA-CA if :attr:`calpha_only` is True).
    * "w(d)": weighting function (prob. of 2 res. to interact given CB distance)
      from `Xu et al. 2009 <https://dx.doi.org/10.1016%2Fj.jmb.2008.06.002>`_.

    The actual QS-score computation in :func:`QSScorer.GetQSScore` implements
    caching. Repeated computations with alternative chain mappings thus become
    faster.

    :param target: Structure designated as "target". Can be fetched from
                   :class:`ost.mol.alg.chain_mapping.MappingResult`
    :type target: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param chem_groups: Groups of chemically equivalent chains in *target*.
                        Can be fetched from
                        :class:`ost.mol.alg.chain_mapping.MappingResult`
    :type chem_groups: :class:`list` of :class:`list` of :class:`str`
    :param model: Structure designated as "model". Can be fetched from
                  :class:`ost.mol.alg.chain_mapping.MappingResult`
    :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param alns: Each alignment is accessible with ``alns[(t_chain,m_chain)]``.
                 First sequence is the sequence of the respective chain in
                 :attr:`~qsent1`, second sequence the one from :attr:`~qsent2`.
                 Can be fetched from
                 :class:`ost.mol.alg.chain_mapping.MappingResult`
    :type alns: :class:`dict` with key: :class:`tuple` of :class:`str`, value:
                :class:`ost.seq.AlignmentHandle`
    """
    def __init__(self, target, chem_groups, model, alns, contact_d = 12.0):

        self._qsent1 = QSEntity(target, contact_d = contact_d)

        # ensure that target chain names match the ones in chem_groups
        chem_group_ch_names = list(itertools.chain.from_iterable(chem_groups))
        if self._qsent1.chain_names != sorted(chem_group_ch_names):
            raise RuntimeError(f"Expect exact same chain names in chem_groups "
                               f"and in target (which is processed to only "
                               f"contain peptides/nucleotides). target: "
                               f"{self._qsent1.chain_names}, chem_groups: "
                               f"{chem_group_ch_names}")

        self._chem_groups = chem_groups
        self._qsent2 = QSEntity(model, contact_d = contact_d)
        self._alns = alns

        # cache for mapped interface scores
        # key: tuple of tuple ((qsent1_ch1, qsent1_ch2),
        #                     ((qsent2_ch1, qsent2_ch2))
        # value: tuple with two numbers referring to QS-score formalism
        #        (equation 6 in Bertoni et al., 2017):
        #        1: contribution of that interface to nominator
        #        2: contribution of that interface to denominator
        self._mapped_cache = dict()

        # cache for non-mapped interfaces in qsent1
        # key: tuple (qsent1_ch1, qsent1_ch2)
        # value: contribution of that interface to second term in denominator in
        #        equation 6 in Bertoni et al., 2017.
        #        (sum of non-shared penalty of qsent1)
        self._qsent_1_penalties = dict()

        # same for qsent2
        self._qsent_2_penalties = dict()

    @staticmethod
    def FromMappingResult(mapping_result):
        """ The preferred way to get a :class:`QSScorer`

        Static constructor that derives an object of type :class:`QSScorer`
        using a :class:`ost.mol.alg.chain_mapping.MappingResult`

        :param mapping_result: Data source
        :type mapping_result: :class:`ost.mol.alg.chain_mapping.MappingResult`
        """
        qs_scorer = QSScorer(mapping_result.target, mapping_result.chem_groups,
                             mapping_result.model, alns = mapping_result.alns)
        return qs_scorer

    @property
    def qsent1(self):
        """ Represents *target*

        :type: :class:`QSEntity`
        """
        return self._qsent1

    @property
    def chem_groups(self):
        """ Groups of chemically equivalent chains in *target*

        Provided at object construction

        :type: :class:`list` of :class:`list` of :class:`str`
        """
        return self._chem_groups

    @property
    def qsent2(self):
        """ Represents *model*

        :type: :class:`QSEntity`
        """
        return self._qsent2

    @property
    def alns(self):
        """ Alignments between chains in :attr:`~qsent1` and :attr:`~qsent2`

        Provided at object construction. Each alignment is accessible with
        ``alns[(t_chain,m_chain)]``. First sequence is the sequence of the
        respective chain in :attr:`~qsent1`, second sequence the one from
        :attr:`~qsent2`.

        :type: :class:`dict` with key: :class:`tuple` of :class:`str`, value:
               :class:`ost.seq.AlignmentHandle`
        """
        return self._alns
    
    def GetQSScore(self, mapping, check=True):
        """ Computes QS-score given chain mapping

        Again, the preferred way is to get this chain mapping is from an object
        of type :class:`ost.mol.alg.chain_mapping.MappingResult`.

        :param mapping: see 
                        :attr:`ost.mol.alg.chain_mapping.MappingResult.mapping`
        :type mapping: :class:`list` of :class:`list` of :class:`str`
        :param check: Perform input checks, can be disabled for speed purposes
                      if you know what you're doing.
        :type check: :class:`bool`
        :returns: The QS-Score
        """

        if check:
            # ensure that dimensionality of mapping matches self.chem_groups
            if len(self.chem_groups) != len(mapping):
                raise RuntimeError("Dimensions of self.chem_groups and mapping "
                                   "must match")
            for a,b in zip(self.chem_groups, mapping):
                if len(a) != len(b):
                    raise RuntimeError("Dimensions of self.chem_groups and "
                                       "mapping must match")
            # ensure that chain names in mapping are all present in qsent2
            for name in itertools.chain.from_iterable(mapping):
                if name is not None and name not in self.qsent2.chain_names:
                    raise RuntimeError(f"Each chain in mapping must be present "
                                       f"in self.qsent2. No match for "
                                       f"\"{name}\"")

        flat_mapping = dict()
        for a, b in zip(self.chem_groups, mapping):
            flat_mapping.update({x: y for x, y in zip(a, b) if y is not None})

        # refers to equation 6 in Bertoni et al., 2017
        nominator, denominator = self._FromFlatMapping(flat_mapping)

        if denominator > 0.0:
            return nominator / denominator
        else:
            return 0.0

    def _FromFlatMapping(self, flat_mapping):

        # refers to equation 6 in Bertoni et al., 2017
        nominator = 0.0
        denominator= 0.0

        # keep track of processed interfaces in qsent2
        processed_qsent2_interfaces = set()

        for int1 in self.qsent1.interacting_chains:
            if int1[0] in flat_mapping and int1[1] in flat_mapping:
                int2 = (flat_mapping[int1[0]], flat_mapping[int1[1]])
                a,b = self._MappedInterfaceScores(int1, int2)
                nominator += a
                denominator += b
                processed_qsent2_interfaces.add((min(int2[0], int2[1]),
                                                 max(int2[0], int2[1])))
            else:
                denominator += self._InterfacePenalty1(int1)

        # add penalties for non-mapped interfaces in qsent2
        for int2 in self.qsent2.interacting_chains:
            if int2 not in processed_qsent2_interfaces:
                denominator += self._InterfacePenalty2(int2)

        return (nominator, denominator)

    def _MappedInterfaceScores(self, int1, int2):
        key_one = (int1, int2)
        if key_one in self._mapped_cache:
            return self._mapped_cache[key_one]
        key_two = ((int1[1], int1[0]), (int2[1], int2[0]))
        if key_two in self._mapped_cache:
            return self._mapped_cache[key_two]
        nominator, denominator = self._InterfaceScores(int1, int2)
        self._mapped_cache[key_one] = (nominator, denominator)
        return (nominator, denominator) 

    def _InterfaceScores(self, int1, int2):

        d1 = self.qsent1.PairDist(int1[0], int1[1])
        d2 = self.qsent2.PairDist(int2[0], int2[1])

        # given two chain names a and b: if a < b, shape of pairwise distances is
        # (len(a), len(b)). However, if b > a, its (len(b), len(a)) => transpose
        if int1[0] > int1[1]:
            d1 = d1.transpose()
        if int2[0] > int2[1]:
            d2 = d2.transpose()

        # indices of the first chain in the two interfaces
        mapped_indices_1_1, mapped_indices_1_2 = \
        self._IndexMapping(int1[0], int2[0])
        # indices of the second chain in the two interfaces
        mapped_indices_2_1, mapped_indices_2_2 = \
        self._IndexMapping(int1[1], int2[1])

        # get shared_masks - for this we first need to select the actual
        # mapped positions to get a one-to-one relationshop and map it back
        # to the original mask size
        mapped_idx_grid_1 = np.ix_(mapped_indices_1_1, mapped_indices_2_1)
        mapped_idx_grid_2 = np.ix_(mapped_indices_1_2, mapped_indices_2_2)
        mapped_d1 = d1[mapped_idx_grid_1]
        mapped_d2 = d2[mapped_idx_grid_2]
        assert(mapped_d1.shape == mapped_d2.shape)
        assert(self.qsent1.contact_d == self.qsent2.contact_d)
        contact_d = self.qsent1.contact_d
        shared_mask = np.logical_and(mapped_d1 < contact_d,
                                     mapped_d2 < contact_d)
        shared_mask_d1 = np.full(d1.shape, False, dtype=bool)
        shared_mask_d1[mapped_idx_grid_1] = shared_mask
        shared_mask_d2 = np.full(d2.shape, False, dtype=bool)
        shared_mask_d2[mapped_idx_grid_2] = shared_mask

        # nominator and denominator contributions from shared contacts
        shared_d1 = d1[shared_mask_d1]
        shared_d2 = d2[shared_mask_d2]
        shared_min = np.minimum(shared_d1, shared_d2)
        shared_abs_diff_div_12 = np.abs(np.subtract(shared_d1, shared_d2))/12.0
        weight_term = np.ones(shared_min.shape[0])
        bigger_5_mask = shared_min > 5.0
        weights = np.exp(-2.0*np.square((shared_min[bigger_5_mask]-5.0)/4.28))
        weight_term[bigger_5_mask] = weights
        diff_term = np.subtract(np.ones(weight_term.shape[0]),
                                shared_abs_diff_div_12)
        nominator = np.sum(np.multiply(weight_term, diff_term))
        denominator = np.sum(weight_term)

        # denominator contribution from non shared contacts in interface one
        contact_distances = d1[np.logical_and(np.logical_not(shared_mask_d1),
                                              d1 < contact_d)]
        bigger_5 = contact_distances[contact_distances > 5]
        denominator += np.sum(np.exp(-2.0*np.square((bigger_5-5.0)/4.28)))
        # add 1.0 for all contact distances <= 5.0
        denominator += contact_distances.shape[0] - bigger_5.shape[0]

        # same for interface two
        contact_distances = d2[np.logical_and(np.logical_not(shared_mask_d2),
                                              d2 < contact_d)]
        bigger_5 = contact_distances[contact_distances > 5]
        denominator += np.sum(np.exp(-2.0*np.square((bigger_5-5.0)/4.28)))
        # add 1.0 for all contact distances <= 5.0
        denominator += contact_distances.shape[0] - bigger_5.shape[0]

        return (nominator, denominator)

    def _IndexMapping(self, ch1, ch2):
        """ Fetches aln and returns indices of (non-)aligned residues

        returns 2 numpy arrays containing the indices of residues in
        ch1 and ch2 which are aligned
        """
        mapped_indices_1 = list()
        mapped_indices_2 = list()
        idx_1 = 0
        idx_2 = 0
        for col in self.alns[(ch1, ch2)]:
            if col[0] != '-' and col[1] != '-':
                mapped_indices_1.append(idx_1)
                mapped_indices_2.append(idx_2)
            if col[0] != '-':
                idx_1 +=1
            if col[1] != '-':
                idx_2 +=1
        return (np.array(mapped_indices_1), np.array(mapped_indices_2))

    def _InterfacePenalty1(self, interface):
        if interface not in self._qsent_1_penalties:
            self._qsent_1_penalties[interface] = \
            self._InterfacePenalty(self.qsent1, interface)
        return self._qsent_1_penalties[interface]

    def _InterfacePenalty2(self, interface):
        if interface not in self._qsent_2_penalties:
            self._qsent_2_penalties[interface] = \
            self._InterfacePenalty(self.qsent2, interface)
        return self._qsent_2_penalties[interface]

    def _InterfacePenalty(self, qsent, interface):
        d = qsent.PairDist(interface[0], interface[1])
        contact_distances = d[d < qsent.contact_d]
        bigger_5 = contact_distances[contact_distances > 5]
        penalty = np.sum(np.exp(-2.0*np.square((bigger_5-5.0)/4.28)))
        # add 1.0 for all contact distances <= 5.0
        penalty += contact_distances.shape[0] - bigger_5.shape[0]
        return penalty

# specify public interface
__all__ = ('QSEntity', 'QSScorer')
