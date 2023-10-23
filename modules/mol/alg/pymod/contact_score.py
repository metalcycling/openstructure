import itertools
import numpy as np

import time
from ost import mol
from ost import geom
from ost import io

class ContactEntity:
    """ Helper object for Contact-score computation
    """
    def __init__(self, ent, contact_d = 5.0, contact_mode="aa"):

        if contact_mode not in ["aa", "repr"]:
            raise RuntimeError("contact_mode must be in [\"aa\", \"repr\"]")

        if contact_mode == "repr":
            for r in ent.residues:
                repr_at = None
                if r.IsPeptideLinking():
                    cb = r.FindAtom("CB")
                    if cb.IsValid():
                        repr_at = cb
                    elif r.GetName() == "GLY":
                        ca = r.FindAtom("CA")
                        if ca.IsValid():
                            repr_at = ca
                elif r.IsNucleotideLinking():
                    c3 = r.FindAtom("C3'")
                    if c3.IsValid():
                        repr_at = c3
                else:
                    raise RuntimeError(f"Only support peptide and nucleotide "
                                       f"residues in \"repr\" contact mode. "
                                       f"Problematic residue: {r}")
                if repr_at is None:
                    raise RuntimeError(f"Residue {r} has no required "
                                       f"representative atom (CB for peptide "
                                       f"residues (CA for GLY) C3' for "
                                       f"nucleotide residues.")

        self._contact_mode = contact_mode

        if self.contact_mode == "aa":
            self._view = ent.CreateFullView()
        elif self.contact_mode == "repr":
            pep_query = "(peptide=true and (aname=\"CB\" or (rname=\"GLY\" and aname=\"CA\")))"
            nuc_query = "(nucleotide=True and aname=\"C3'\")"
            self._view = ent.Select(" or ".join([pep_query, nuc_query]))
        self._contact_d = contact_d

        # the following attributes will be lazily evaluated
        self._chain_names = None
        self._interacting_chains = None
        self._sequence = dict()
        self._contacts = None
        self._hr_contacts = None
        self._interface_residues = None
        self._hr_interface_residues = None

    @property
    def view(self):
        """ The structure depending on *contact_mode*

        Full view in case of "aa", view that only contains representative
        atoms in case of "repr".

        :type: :class:`ost.mol.EntityView`
        """
        return self._view

    @property
    def contact_mode(self):
        """ The contact mode

        Can either be "aa", meaning that all atoms are considered to identify
        contacts, or "repr" which only considers distances between
        representative atoms. For peptides thats CB (CA for GLY), for
        nucleotides thats C3'.

        :type: :class:`str`
        """
        return self._contact_mode
    
    @property
    def contact_d(self):
        """ Pairwise distance of residues to be considered as contacts

        Given at :class:`ContactScorer` construction

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
            self._interacting_chains = list(self.contacts.keys())
        return self._interacting_chains
    
    @property
    def contacts(self):
        """ Interchain contacts

        Organized as :class:`dict` with key (cname1, cname2) and values being
        a set of tuples with the respective residue indices. 
        cname1 < cname2 evaluates to True.
        """
        if self._contacts is None:
            self._SetupContacts()
        return self._contacts

    @property
    def hr_contacts(self):
        """ Human readable interchain contacts

        Human readable version of :attr:`~contacts`. Simple list with tuples
        containing two strings specifying the residues in contact. Format:
        <cname>.<rnum>.<ins_code>
        """
        if self._hr_contacts is None:
            self._SetupContacts()
        return self._hr_contacts

    @property
    def interface_residues(self):
        """ Interface residues

        Residues in each chain that are in contact with any other chain.
        Organized as :class:`dict` with key cname and values the respective
        residue indices in a :class:`set`.
        """
        if self._interface_residues is None:
            self._SetupInterfaceResidues()
        return self._interface_residues

    @property
    def hr_interface_residues(self):
        """ Human readable interface residues

        Human readable version of :attr:`interface_residues`. :class:`list` of
        strings specifying the interface residues in format:
        <cname>.<rnum>.<ins_code>
        """
        if self._interface_residues is None:
            self._SetupHRInterfaceResidues()
        return self._hr_interface_residues
    
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

    def _SetupContacts(self):
        # this function is incredibly inefficient... if performance is an issue,
        # go ahead and optimize
        self._contacts = dict()
        self._hr_contacts = list()

        # set indices relative to full view 
        for ch in self.view.chains:
            for r_idx, r in enumerate(ch.residues):
                r.SetIntProp("contact_idx", r_idx)

        for cname in self.chain_names:
            # q1 selects stuff in current chain that is close to any other chain
            q1 = f"cname={cname} and {self.contact_d} <> [cname!={cname}]"
            # q2 selects stuff in other chains that is close to current chain
            q2 = f"cname!={cname} and {self.contact_d} <> [cname={cname}]"
            v1 = self.view.Select(q1)
            v2 = self.view.Select(q2)
            v1_p = [geom.Vec3List([a.pos for a in r.atoms]) for r in v1.residues]
            for r1, p1 in zip(v1.residues, v1_p):
                for ch2 in v2.chains:
                    cname2 = ch2.GetName()
                    if cname2 > cname:
                        v2_p = [geom.Vec3List([a.pos for a in r.atoms]) for r in ch2.residues]
                        for r2, p2 in zip(ch2.residues, v2_p):
                            if p1.IsWithin(p2, self.contact_d):
                                cname_key = (cname, cname2)
                                if cname_key not in self._contacts:
                                    self._contacts[cname_key] = set()
                                self._contacts[cname_key].add((r1.GetIntProp("contact_idx"),
                                                               r2.GetIntProp("contact_idx")))
                                rnum1 = r1.GetNumber()
                                hr1 = f"{cname}.{rnum1.num}.{rnum1.ins_code}"
                                rnum2 = r2.GetNumber()
                                hr2 = f"{cname2}.{rnum2.num}.{rnum2.ins_code}"
                                self._hr_contacts.append((hr1.strip("\u0000"),
                                                          hr2.strip("\u0000")))

    def _SetupInterfaceResidues(self):
        self._interface_residues = {cname: set() for cname in self.chain_names}
        for k,v in self.contacts.items():
            for item in v:
                self._interface_residues[k[0]].add(item[0])
                self._interface_residues[k[1]].add(item[1])

    def _SetupHRInterfaceResidues(self):
        interface_residues = set()
        for item in self.hr_contacts:
            interface_residues.add(item[0])
            interface_residues.add(item[1])
        self._hr_interface_residues = list(interface_residues)


class ContactScorerResultICS:
    """
    Holds data relevant to compute ics
    """
    def __init__(self, n_trg_contacts, n_mdl_contacts, n_union, n_intersection):
        self._n_trg_contacts = n_trg_contacts
        self._n_mdl_contacts = n_mdl_contacts
        self._n_union = n_union
        self._n_intersection = n_intersection

    @property
    def n_trg_contacts(self):
        """ Number of contacts in target

        :type: :class:`int`
        """
        return self._n_trg_contacts

    @property
    def n_mdl_contacts(self):
        """ Number of contacts in model

        :type: :class:`int`
        """
        return self._n_mdl_contacts

    @property
    def precision(self):
        """ Precision of model contacts

        The fraction of model contacts that are also present in target

        :type: :class:`int`
        """
        if self._n_mdl_contacts != 0:
            return self._n_intersection / self._n_mdl_contacts
        else:
            return 0.0

    @property
    def recall(self):
        """ Recall of model contacts

        The fraction of target contacts that are also present in model

        :type: :class:`int`
        """
        if self._n_trg_contacts != 0:
            return self._n_intersection / self._n_trg_contacts
        else:
            return 0.0

    @property
    def ics(self):
        """ The Interface Contact Similarity score (ICS)

        Combination of :attr:`precision` and :attr:`recall` using the F1-measure

        :type: :class:`float`
        """
        p = self.precision
        r = self.recall
        nominator = p*r
        denominator = p + r
        if denominator != 0.0:
            return 2*nominator/denominator
        else:
            return 0.0

class ContactScorerResultIPS:
    """
    Holds data relevant to compute ips
    """
    def __init__(self, n_trg_int_res, n_mdl_int_res, n_union, n_intersection):
        self._n_trg_int_res = n_trg_int_res
        self._n_mdl_int_res = n_mdl_int_res
        self._n_union = n_union
        self._n_intersection = n_intersection

    @property
    def n_trg_int_res(self):
        """ Number of interface residues in target

        :type: :class:`int`
        """
        return self._n_trg_contacts

    @property
    def n_mdl_int_res(self):
        """ Number of interface residues in model

        :type: :class:`int`
        """
        return self._n_mdl_int_res

    @property
    def precision(self):
        """ Precision of model interface residues

        The fraction of model interface residues that are also interface
        residues in target

        :type: :class:`int`
        """
        if self._n_mdl_int_res != 0:
            return self._n_intersection / self._n_mdl_int_res
        else:
            return 0.0

    @property
    def recall(self):
        """ Recall of model interface residues

        The fraction of target interface residues that are also interface
        residues in model

        :type: :class:`int`
        """
        if self._n_trg_int_res != 0:
            return self._n_intersection / self._n_trg_int_res
        else:
            return 0.0

    @property
    def ips(self):
        """ The Interface Patch Similarity score (IPS)

        Jaccard coefficient of interface residues in model/target.
        Technically thats :attr:`intersection`/:attr:`union` 

        :type: :class:`float`
        """
        if(self._n_union > 0):
            return self._n_intersection/self._n_union
        return 0.0

class ContactScorer:
    """ Helper object to compute Contact scores

    Tightly integrated into the mechanisms from the chain_mapping module.
    The prefered way to derive an object of type :class:`ContactScorer` is
    through the static constructor: :func:`~FromMappingResult`.

    Usage is the same as for :class:`ost.mol.alg.QSScorer`
    """

    def __init__(self, target, chem_groups, model, alns,
                 contact_mode="aa", contact_d=5.0):
        self._cent1 = ContactEntity(target, contact_mode = contact_mode,
                                    contact_d = contact_d)
        # ensure that target chain names match the ones in chem_groups
        chem_group_ch_names = list(itertools.chain.from_iterable(chem_groups))
        if self._cent1.chain_names != sorted(chem_group_ch_names):
            raise RuntimeError(f"Expect exact same chain names in chem_groups "
                               f"and in target (which is processed to only "
                               f"contain peptides/nucleotides). target: "
                               f"{self._cent1.chain_names}, chem_groups: "
                               f"{chem_group_ch_names}")

        self._chem_groups = chem_groups
        self._cent2 = ContactEntity(model, contact_mode = contact_mode,
                                    contact_d = contact_d)
        self._alns = alns

        # cache for mapped interface scores
        # relevant to compute ICS
        # key: tuple of tuple ((qsent1_ch1, qsent1_ch2),
        #                     ((qsent2_ch1, qsent2_ch2))
        # value: tuple with two numbers required for computation of ICS
        #        1: n_union
        #        2: n_intersection
        self._mapped_cache_ics = dict()

        # cache for mapped scores
        # relevant to compute IPS
        # key: tuple: (qsent1_ch, qsent2_ch)
        # value: tuple with two numbers required for computation of IPS
        #        1: n_union
        #        2: n_intersection
        self._mapped_cache_ips = dict()

    @staticmethod
    def FromMappingResult(mapping_result, contact_mode="aa", contact_d = 5.0):
        """ The preferred way to get a :class:`ContactScorer`

        Static constructor that derives an object of type :class:`ContactScorer`
        using a :class:`ost.mol.alg.chain_mapping.MappingResult`

        :param mapping_result: Data source
        :type mapping_result: :class:`ost.mol.alg.chain_mapping.MappingResult`
        """
        contact_scorer = ContactScorer(mapping_result.target,
                                       mapping_result.chem_groups,
                                       mapping_result.model,
                                       mapping_result.alns,
                                       contact_mode = contact_mode,
                                       contact_d = contact_d)
        return contact_scorer

    @property
    def cent1(self):
        """ Represents *target*

        :type: :class:`ContactEntity`
        """
        return self._cent1

    @property
    def chem_groups(self):
        """ Groups of chemically equivalent chains in *target*

        Provided at object construction

        :type: :class:`list` of :class:`list` of :class:`str`
        """
        return self._chem_groups

    @property
    def cent2(self):
        """ Represents *model*

        :type: :class:`ContactEntity`
        """
        return self._cent2

    @property
    def alns(self):
        """ Alignments between chains in :attr:`~cent1` and :attr:`~cent2`

        Provided at object construction. Each alignment is accessible with
        ``alns[(t_chain,m_chain)]``. First sequence is the sequence of the
        respective chain in :attr:`~cent1`, second sequence the one from
        :attr:`~cent2`.

        :type: :class:`dict` with key: :class:`tuple` of :class:`str`, value:
               :class:`ost.seq.AlignmentHandle`
        """
        return self._alns

    def ScoreICS(self, mapping, check=True):
        """ Computes ICS given chain mapping

        Again, the preferred way is to get *mapping* is from an object
        of type :class:`ost.mol.alg.chain_mapping.MappingResult`.

        :param mapping: see 
                        :attr:`ost.mol.alg.chain_mapping.MappingResult.mapping`
        :type mapping: :class:`list` of :class:`list` of :class:`str`
        :param check: Perform input checks, can be disabled for speed purposes
                      if you know what you're doing.
        :type check: :class:`bool`
        :returns: Result object of type :class:`ContactScorerResultICS`
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
            # ensure that chain names in mapping are all present in cent2
            for name in itertools.chain.from_iterable(mapping):
                if name is not None and name not in self.cent2.chain_names:
                    raise RuntimeError(f"Each chain in mapping must be present "
                                       f"in self.cent2. No match for "
                                       f"\"{name}\"")

        flat_mapping = dict()
        for a, b in zip(self.chem_groups, mapping):
            flat_mapping.update({x: y for x, y in zip(a, b) if y is not None})

        return self.ICSFromFlatMapping(flat_mapping)

    def ScoreICSInterface(self, trg_ch1, trg_ch2, mdl_ch1, mdl_ch2):
        """ Computes ICS scores only considering one interface

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
        :returns: Result object of type :class:`ContactScorerResultICS`
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
        trg_int_r = (trg_ch2, trg_ch1)
        mdl_int_r = (mdl_ch2, mdl_ch1)

        if trg_int in self.cent1.contacts:
            n_trg = len(self.cent1.contacts[trg_int])
        elif trg_int_r in self.cent1.contacts:
            n_trg = len(self.cent1.contacts[trg_int_r])
        else:
            n_trg = 0

        if mdl_int in self.cent2.contacts:
            n_mdl = len(self.cent2.contacts[mdl_int])
        elif mdl_int_r in self.cent2.contacts:
            n_mdl = len(self.cent2.contacts[mdl_int_r])
        else:
            n_mdl = 0

        n_union, n_intersection = self._MappedInterfaceScores(trg_int, mdl_int)
        return ContactScorerResultICS(n_trg, n_mdl, n_union, n_intersection)

    def ICSFromFlatMapping(self, flat_mapping):
        """ Same as :func:`ScoreICS` but with flat mapping

        :param flat_mapping: Dictionary with target chain names as keys and
                             the mapped model chain names as value
        :type flat_mapping: :class:`dict` with :class:`str` as key and value
        :returns: Result object of type :class:`ContactScorerResultICS`
        """
        n_trg = sum([len(x) for x in self.cent1.contacts.values()])
        n_mdl = sum([len(x) for x in self.cent2.contacts.values()])
        n_union = 0
        n_intersection = 0

        processed_cent2_interfaces = set()
        for int1 in self.cent1.interacting_chains:
            if int1[0] in flat_mapping and int1[1] in flat_mapping:
                int2 = (flat_mapping[int1[0]], flat_mapping[int1[1]])
                a, b = self._MappedInterfaceScores(int1, int2)
                n_union += a
                n_intersection += b
                processed_cent2_interfaces.add((min(int2), max(int2)))

        # process interfaces that only exist in qsent2
        r_flat_mapping = {v:k for k,v in flat_mapping.items()} # reverse mapping
        for int2 in self.cent2.interacting_chains:
            if int2 not in processed_cent2_interfaces:
                if int2[0] in r_flat_mapping and int2[1] in r_flat_mapping:
                    int1 = (r_flat_mapping[int2[0]], r_flat_mapping[int2[1]])
                    a, b = self._MappedInterfaceScores(int1, int2)
                    n_union += a
                    n_intersection += b

        return ContactScorerResultICS(n_trg, n_mdl,
                                      n_union, n_intersection)

    def ScoreIPS(self, mapping, check=True):
        """ Computes IPS given chain mapping

        Again, the preferred way is to get *mapping* is from an object
        of type :class:`ost.mol.alg.chain_mapping.MappingResult`.

        :param mapping: see 
                        :attr:`ost.mol.alg.chain_mapping.MappingResult.mapping`
        :type mapping: :class:`list` of :class:`list` of :class:`str`
        :param check: Perform input checks, can be disabled for speed purposes
                      if you know what you're doing.
        :type check: :class:`bool`
        :returns: Result object of type :class:`ContactScorerResultIPS`
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
            # ensure that chain names in mapping are all present in cent2
            for name in itertools.chain.from_iterable(mapping):
                if name is not None and name not in self.cent2.chain_names:
                    raise RuntimeError(f"Each chain in mapping must be present "
                                       f"in self.cent2. No match for "
                                       f"\"{name}\"")

        flat_mapping = dict()
        for a, b in zip(self.chem_groups, mapping):
            flat_mapping.update({x: y for x, y in zip(a, b) if y is not None})

        return self.IPSFromFlatMapping(flat_mapping)

    def IPSFromFlatMapping(self, flat_mapping):
        """ Same as :func:`ScoreIPS` but with flat mapping

        :param flat_mapping: Dictionary with target chain names as keys and
                             the mapped model chain names as value
        :type flat_mapping: :class:`dict` with :class:`str` as key and value
        :returns: Result object of type :class:`ContactScorerResultIPS`
        """
        n_trg = sum([len(x) for x in self.cent1.interface_residues.values()])
        n_mdl = sum([len(x) for x in self.cent2.interface_residues.values()])
        n_union = 0
        n_intersection = 0

        processed_cent2_chains = set()
        for trg_ch in self.cent1.chain_names:
            if trg_ch in flat_mapping:
                a, b = self._MappedSCScores(trg_ch, flat_mapping[trg_ch])
                n_union += a
                n_intersection += b
                processed_cent2_chains.add(flat_mapping[trg_ch])
            else:
                n_union += len(self.cent1.interface_residues[trg_ch])

        for mdl_ch in self._cent2.chain_names:
            if mdl_ch not in processed_cent2_chains:
                n_union += len(self.cent2.interface_residues[mdl_ch])

        return ContactScorerResultIPS(n_trg, n_mdl,
                                      n_union, n_intersection)


    def _MappedInterfaceScores(self, int1, int2):
        key_one = (int1, int2)
        if key_one in self._mapped_cache_ics:
            return self._mapped_cache_ics[key_one]
        key_two = ((int1[1], int1[0]), (int2[1], int2[0]))
        if key_two in self._mapped_cache_ics:
            return self._mapped_cache_ics[key_two]

        n_union, n_intersection = self._InterfaceScores(int1, int2)
        self._mapped_cache_ics[key_one] = (n_union, n_intersection)
        return (n_union, n_intersection)

    def _InterfaceScores(self, int1, int2):
        if int1 in self.cent1.contacts:
            ref_contacts = self.cent1.contacts[int1]
        elif (int1[1], int1[0]) in self.cent1.contacts:
            ref_contacts = self.cent1.contacts[(int1[1], int1[0])]
            # need to reverse contacts
            ref_contacts = set([(x[1], x[0]) for x in ref_contacts])
        else:
            ref_contacts = set() # no contacts at all

        if int2 in self.cent2.contacts:
            mdl_contacts = self.cent2.contacts[int2]
        elif (int2[1], int2[0]) in self.cent2.contacts:
            mdl_contacts = self.cent2.contacts[(int2[1], int2[0])]
            # need to reverse contacts
            mdl_contacts = set([(x[1], x[0]) for x in mdl_contacts])
        else:
            mdl_contacts = set() # no contacts at all

        # indices in contacts lists are specific to the respective
        # structures, need manual mapping from alignments
        ch1_aln = self.alns[(int1[0], int2[0])]
        ch2_aln = self.alns[(int1[1], int2[1])]
        mapped_ref_contacts = set()
        mapped_mdl_contacts = set()
        for c in ref_contacts:
            mapped_c = (ch1_aln.GetPos(0, c[0]), ch2_aln.GetPos(0, c[1]))
            mapped_ref_contacts.add(mapped_c)
        for c in mdl_contacts:
            mapped_c = (ch1_aln.GetPos(1, c[0]), ch2_aln.GetPos(1, c[1]))
            mapped_mdl_contacts.add(mapped_c)

        return (len(mapped_ref_contacts.union(mapped_mdl_contacts)),
                len(mapped_ref_contacts.intersection(mapped_mdl_contacts)))

    def _MappedSCScores(self, ref_ch, mdl_ch):
        if (ref_ch, mdl_ch) in self._mapped_cache_ips:
            return self._mapped_cache_ips[(ref_ch, mdl_ch)]
        n_union, n_intersection = self._SCScores(ref_ch, mdl_ch)
        self._mapped_cache_ips[(ref_ch, mdl_ch)] = (n_union, n_intersection)
        return (n_union, n_intersection)

    def _SCScores(self, ch1, ch2):
        ref_int_res = self.cent1.interface_residues[ch1]
        mdl_int_res = self.cent2.interface_residues[ch2]
        aln = self.alns[(ch1, ch2)]
        mapped_ref_int_res = set()
        mapped_mdl_int_res = set()
        for r_idx in ref_int_res:
            mapped_ref_int_res.add(aln.GetPos(0, r_idx))
        for r_idx in mdl_int_res:
            mapped_mdl_int_res.add(aln.GetPos(1, r_idx))
        return(len(mapped_ref_int_res.union(mapped_mdl_int_res)),
               len(mapped_ref_int_res.intersection(mapped_mdl_int_res)))

# specify public interface
__all__ = ('ContactEntity', 'ContactScorerResultICS', 'ContactScorerResultIPS', 'ContactScorer')
