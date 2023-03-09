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

class QSScorerResult:
    """
    Holds data relevant for QS-score computation. Formulas for QS scores:

    ::

      - QS_best = weighted_scores / (weight_sum + weight_extra_mapped)
      - QS_global = weighted_scores / (weight_sum + weight_extra_all)
      -> weighted_scores = sum(w(min(d1,d2)) * (1 - abs(d1-d2)/12)) for shared
      -> weight_sum = sum(w(min(d1,d2))) for shared
      -> weight_extra_mapped = sum(w(d)) for all mapped but non-shared
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
    """
    def __init__(self, weighted_scores, weight_sum, weight_extra_mapped,
                 weight_extra_all):
        self._weighted_scores = weighted_scores
        self._weight_sum = weight_sum
        self._weight_extra_mapped = weight_extra_mapped
        self._weight_extra_all = weight_extra_all

    @property
    def weighted_scores(self):
        """ weighted_scores attribute as described in formula section above

        :type: :class:`float`
        """
        return self._weighted_scores

    @property
    def weight_sum(self):
        """ weight_sum attribute as described in formula section above

        :type: :class:`float`
        """
        return self._weight_sum

    @property
    def weight_extra_mapped(self):
        """ weight_extra_mapped attribute as described in formula section above

        :type: :class:`float`
        """
        return self._weight_extra_mapped

    @property
    def weight_extra_all(self):
        """ weight_extra_all attribute as described in formula section above

        :type: :class:`float`
        """
        return self._weight_extra_all

    @property
    def QS_best(self):
        """ QS_best - the actual score as described in formula section above

        :type: :class:`float`
        """
        nominator = self.weighted_scores
        denominator = self.weight_sum + self.weight_extra_mapped
        if denominator != 0.0:
            return nominator/denominator
        else:
            return 0.0

    @property
    def QS_global(self):
        """ QS_global - the actual score as described in formula section above

        :type: :class:`float`
        """
        nominator = self.weighted_scores
        denominator = self.weight_sum + self.weight_extra_all
        if denominator != 0.0:
            return nominator/denominator
        else:
            return 0.0


class QSScorer:
    """ Helper object to compute QS-score

    Tightly integrated into the mechanisms from the chain_mapping module.
    The prefered way to derive an object of type :class:`QSScorer` is through
    the static constructor: :func:`~FromMappingResult`. Example score
    computation including mapping:

    ::

        from ost.mol.alg.qsscore import QSScorer
        from ost.mol.alg.chain_mapping import ChainMapper

        ent_1 = io.LoadPDB("path_to_assembly_1.pdb")
        ent_2 = io.LoadPDB("path_to_assembly_2.pdb")

        chain_mapper = ChainMapper(ent_1)
        mapping_result = chain_mapper.GetlDDTMapping(ent_2)
        qs_scorer = QSScorer.FromMappingResult(mapping_result)
        score_result = qs_scorer.Score(mapping_result.mapping)
        print("score:", score_result.QS_global)

    QS-score computation in :func:`QSScorer.Score` implements caching.
    Repeated computations with alternative chain mappings thus become faster.

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
        # value: tuple with four numbers referring to QS-score formalism
        #        1: weighted_scores
        #        2: weight_sum
        #        3: weight_extra_mapped
        #        4: weight_extra_all
        self._mapped_cache = dict()

        # cache for non-mapped interfaces in qsent1
        # key: tuple (qsent1_ch1, qsent1_ch2)
        # value: contribution of that interface to weight_extra_all
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
    
    def Score(self, mapping, check=True):
        """ Computes QS-score given chain mapping

        Again, the preferred way is to get *mapping* is from an object
        of type :class:`ost.mol.alg.chain_mapping.MappingResult`.

        :param mapping: see 
                        :attr:`ost.mol.alg.chain_mapping.MappingResult.mapping`
        :type mapping: :class:`list` of :class:`list` of :class:`str`
        :param check: Perform input checks, can be disabled for speed purposes
                      if you know what you're doing.
        :type check: :class:`bool`
        :returns: Result object of type :class:`QSScorerResult`
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

        return self.FromFlatMapping(flat_mapping)

    def ScoreInterface(self, trg_ch1, trg_ch2, mdl_ch1, mdl_ch2):
        """ Computes QS-score only considering one interface

        This only works for interfaces that are computed in :func:`Score`, i.e.
        interfaces for which the alignments are set up correctly.

        :param trg_ch1: Name of first interface chain in target
        :type trg_ch1: :class:`str`
        :param trg_ch2: Name of second interface chain in target
        :type trg_ch2: :class:`str`
        :param mdl_ch1: Name of first interface chain in model
        :type mdl_ch1: :class:`str`
        :param mdl_ch2: Name of second interface chain in model
        :type mdl_ch2: :class:`str`
        :returns: Result object of type :class:`QSScorerResult`
        :raises: :class:`RuntimeError` if no aln for trg_ch1/mdl_ch1 or
                 trg_ch2/mdl_ch2 is available.
        """
        if (trg_ch1, mdl_ch1) not in self.alns:
            raise RuntimeError(f"No aln between trg_ch1 ({trg_ch1}) and "
                               f"mdl_ch1 ({mdl_ch1}) available. Did you "
                               f"construct the QSScorer object from a "
                               f"MappingResult and are trg_ch1 and mdl_ch1 "
                               f"mapped to each other?")
        if (trg_ch2, mdl_ch2) not in self.alns:
            raise RuntimeError(f"No aln between trg_ch1 ({trg_ch1}) and "
                               f"mdl_ch1 ({mdl_ch1}) available. Did you "
                               f"construct the QSScorer object from a "
                               f"MappingResult and are trg_ch1 and mdl_ch1 "
                               f"mapped to each other?")
        trg_int = (trg_ch1, trg_ch2)
        mdl_int = (mdl_ch1, mdl_ch2)
        a, b, c, d = self._MappedInterfaceScores(trg_int, mdl_int)
        return QSScorerResult(a, b, c, d)

    def FromFlatMapping(self, flat_mapping):
        """ Same as :func:`Score` but with flat mapping

        :param flat_mapping: Dictionary with target chain names as keys and
                             the mapped model chain names as value
        :type flat_mapping: :class:`dict` with :class:`str` as key and value
        :returns: Result object of type :class:`QSScorerResult`
        """

        weighted_scores = 0.0
        weight_sum = 0.0
        weight_extra_mapped = 0.0
        weight_extra_all = 0.0

        # keep track of processed interfaces in qsent2
        processed_qsent2_interfaces = set()

        for int1 in self.qsent1.interacting_chains:
            if int1[0] in flat_mapping and int1[1] in flat_mapping:
                int2 = (flat_mapping[int1[0]], flat_mapping[int1[1]])
                a, b, c, d = self._MappedInterfaceScores(int1, int2)
                weighted_scores += a
                weight_sum += b
                weight_extra_mapped += c
                weight_extra_all += d
                processed_qsent2_interfaces.add((min(int2[0], int2[1]),
                                                 max(int2[0], int2[1])))
            else:
                weight_extra_all += self._InterfacePenalty1(int1)

        # process interfaces that only exist in qsent2
        r_flat_mapping = {v:k for k,v in flat_mapping.items()} # reverse mapping...
        for int2 in self.qsent2.interacting_chains:
            if int2 not in processed_qsent2_interfaces:
                if int2[0] in r_flat_mapping and int2[1] in r_flat_mapping:
                    int1 = (r_flat_mapping[int2[0]], r_flat_mapping[int2[1]])
                    a, b, c, d = self._MappedInterfaceScores(int1, int2)
                    weighted_scores += a
                    weight_sum += b
                    weight_extra_mapped += c
                    weight_extra_all += d
                else:
                    weight_extra_all += self._InterfacePenalty2(int2)

        return QSScorerResult(weighted_scores, weight_sum, weight_extra_mapped,
                              weight_extra_all)

    def _MappedInterfaceScores(self, int1, int2):
        key_one = (int1, int2)
        if key_one in self._mapped_cache:
            return self._mapped_cache[key_one]
        key_two = ((int1[1], int1[0]), (int2[1], int2[0]))
        if key_two in self._mapped_cache:
            return self._mapped_cache[key_two]

        weighted_scores, weight_sum, weight_extra_mapped, weight_extra_all = \
        self._InterfaceScores(int1, int2)
        self._mapped_cache[key_one] = (weighted_scores, weight_sum, weight_extra_mapped,
                                       weight_extra_all)
        return (weighted_scores, weight_sum, weight_extra_mapped, weight_extra_all)

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
        # mapped positions to get a one-to-one relationship and map it back
        # to the original mask size
        assert(self.qsent1.contact_d == self.qsent2.contact_d)
        contact_d = self.qsent1.contact_d
        mapped_idx_grid_1 = np.ix_(mapped_indices_1_1, mapped_indices_2_1)
        mapped_idx_grid_2 = np.ix_(mapped_indices_1_2, mapped_indices_2_2)
        mapped_d1_contacts = d1[mapped_idx_grid_1] < contact_d
        mapped_d2_contacts = d2[mapped_idx_grid_2] < contact_d
        assert(mapped_d1_contacts.shape == mapped_d2_contacts.shape)
        shared_mask = np.logical_and(mapped_d1_contacts, mapped_d2_contacts)
        shared_mask_d1 = np.full(d1.shape, False, dtype=bool)
        shared_mask_d1[mapped_idx_grid_1] = shared_mask
        shared_mask_d2 = np.full(d2.shape, False, dtype=bool)
        shared_mask_d2[mapped_idx_grid_2] = shared_mask

        # get mapped but nonshared masks
        mapped_nonshared_mask_d1 = np.full(d1.shape, False, dtype=bool)
        mapped_nonshared_mask_d1[mapped_idx_grid_1] = \
        np.logical_and(np.logical_not(shared_mask), mapped_d1_contacts)
        mapped_nonshared_mask_d2 = np.full(d2.shape, False, dtype=bool)
        mapped_nonshared_mask_d2[mapped_idx_grid_2] = \
        np.logical_and(np.logical_not(shared_mask), mapped_d2_contacts)

        # contributions from shared contacts
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
        weighted_scores = np.sum(np.multiply(weight_term, diff_term))
        weight_sum = np.sum(weight_term)

        # do weight_extra_all for interface one
        nonshared_contact_mask_d1 = np.logical_and(np.logical_not(shared_mask_d1),
                                                   d1 < contact_d)
        contact_distances = d1[nonshared_contact_mask_d1]
        bigger_5 = contact_distances[contact_distances > 5]
        weight_extra_all = np.sum(np.exp(-2.0*np.square((bigger_5-5.0)/4.28)))
        # add 1.0 for all contact distances <= 5.0
        weight_extra_all += contact_distances.shape[0] - bigger_5.shape[0]
        # same for interface two
        nonshared_contact_mask_d2 = np.logical_and(np.logical_not(shared_mask_d2),
                                                   d2 < contact_d)
        contact_distances = d2[nonshared_contact_mask_d2]
        bigger_5 = contact_distances[contact_distances > 5]
        weight_extra_all += np.sum(np.exp(-2.0*np.square((bigger_5-5.0)/4.28)))
        # add 1.0 for all contact distances <= 5.0
        weight_extra_all += contact_distances.shape[0] - bigger_5.shape[0]

        # do weight_extra_mapped for interface one
        contact_distances = d1[mapped_nonshared_mask_d1]
        bigger_5 = contact_distances[contact_distances > 5]
        weight_extra_mapped = np.sum(np.exp(-2.0*np.square((bigger_5-5.0)/4.28)))
        # add 1.0 for all contact distances <= 5.0
        weight_extra_mapped += contact_distances.shape[0] - bigger_5.shape[0]
        # same for interface two
        contact_distances = d2[mapped_nonshared_mask_d2]
        bigger_5 = contact_distances[contact_distances > 5]
        weight_extra_mapped += np.sum(np.exp(-2.0*np.square((bigger_5-5.0)/4.28)))
        # add 1.0 for all contact distances <= 5.0
        weight_extra_mapped += contact_distances.shape[0] - bigger_5.shape[0]

        return (weighted_scores, weight_sum, weight_extra_mapped, weight_extra_all)

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
