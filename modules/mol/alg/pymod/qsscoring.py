"""
Scoring of quaternary structures (QS). The QS scoring is according to the paper
by `Bertoni et al. <https://dx.doi.org/10.1038/s41598-017-09654-8>`_.

.. note ::

  Requirements for use:

  - A default :class:`compound library <ost.conop.CompoundLib>` must be defined
    and accessible via :func:`~ost.conop.GetDefaultLib`. This is set by default
    when executing scripts with ``ost``. Otherwise, you must set this with
    :func:`~ost.conop.SetDefaultLib`.
  - ClustalW must be installed (unless you provide chain mappings)
  - Python modules `numpy` and `scipy` must be installed and available
    (e.g. use ``pip install scipy numpy``)
"""
# Original authors: Gerardo Tauriello, Martino Bertoni

from ost import mol, geom, conop, seq, settings, PushVerbosityLevel
from ost import LogError, LogWarning, LogScript, LogInfo, LogVerbose, LogDebug
from ost.bindings.clustalw import ClustalW
from ost.mol.alg import lDDTScorer
from ost.seq.alg.renumber import Renumber
import numpy as np
from scipy.special import factorial
from scipy.special import binom
from scipy.cluster.hierarchy import fclusterdata
import itertools

###############################################################################
# QS scoring
###############################################################################

class QSscoreError(Exception):
  """Exception to be raised for "acceptable" exceptions in QS scoring.

  Those are cases we might want to capture for default behavior.
  """
  pass

class QSscorer:
  """Object to compute QS scores.

  Simple usage without any precomputed contacts, symmetries and mappings:

  .. code-block:: python

    import ost
    from ost.mol.alg import qsscoring

    # load two biounits to compare
    ent_full = ost.io.LoadPDB('3ia3', remote=True)
    ent_1 = ent_full.Select('cname=A,D')
    ent_2 = ent_full.Select('cname=B,C')
    # get score
    ost.PushVerbosityLevel(3)
    try:
      qs_scorer = qsscoring.QSscorer(ent_1, ent_2)
      ost.LogScript('QSscore:', str(qs_scorer.global_score))
      ost.LogScript('Chain mapping used:', str(qs_scorer.chain_mapping))
      # commonly you want the QS global score as output
      qs_score = qs_scorer.global_score
    except qsscoring.QSscoreError as ex:
      # default handling: report failure and set score to 0
      ost.LogError('QSscore failed:', str(ex))
      qs_score = 0

  For maximal performance when computing QS scores of the same entity with many
  others, it is advisable to construct and reuse :class:`QSscoreEntity` objects.

  Any known / precomputed information can be filled into the appropriate
  attribute here (no checks done!). Otherwise most quantities are computed on
  first access and cached (lazy evaluation). Setters are provided to set values
  with extra checks (e.g. :func:`SetSymmetries`).

  All necessary seq. alignments are done by global BLOSUM62-based alignment. A
  multiple sequence alignment is performed with ClustalW unless
  :attr:`chain_mapping` is provided manually. You will need to have an
  executable ``clustalw`` or ``clustalw2`` in your ``PATH`` or you must set
  :attr:`clustalw_bin` accordingly. Otherwise an exception
  (:class:`ost.settings.FileNotFound`) is thrown.

  Formulas for QS scores:

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
  
  :param ent_1: First structure to be scored.
  :type ent_1:  :class:`QSscoreEntity`, :class:`~ost.mol.EntityHandle` or
                :class:`~ost.mol.EntityView`
  :param ent_2: Second structure to be scored.
  :type ent_2:  :class:`QSscoreEntity`, :class:`~ost.mol.EntityHandle` or
                :class:`~ost.mol.EntityView`
  :param res_num_alignment: Sets :attr:`res_num_alignment`

  :raises: :class:`QSscoreError` if input structures are invalid or are monomers
           or have issues that make it impossible for a QS score to be computed.

  .. attribute:: qs_ent_1

    :class:`QSscoreEntity` object for *ent_1* given at construction.
    If entity names (:attr:`~QSscoreEntity.original_name`) are not unique, we
    set it to 'pdb_1' using :func:`~QSscoreEntity.SetName`.

  .. attribute:: qs_ent_2

    :class:`QSscoreEntity` object for *ent_2* given at construction.
    If entity names (:attr:`~QSscoreEntity.original_name`) are not unique, we
    set it to 'pdb_2' using :func:`~QSscoreEntity.SetName`.

  .. attribute:: calpha_only

    True if any of the two structures is CA-only (after cleanup).

    :type: :class:`bool`

  .. attribute:: max_ca_per_chain_for_cm

    Maximal number of CA atoms to use in each chain to determine chain mappings.
    Setting this to -1 disables the limit. Limiting it speeds up determination
    of symmetries and chain mappings. By default it is set to 100.

    :type: :class:`int`

  .. attribute:: max_mappings_extensive

    Maximal number of chain mappings to test for 'extensive'
    :attr:`chain_mapping_scheme`. The extensive chain mapping search must in the
    worst case check O(N^2) * O(N!) possible mappings for complexes with N
    chains. Two octamers without symmetry would require 322560 mappings to be
    checked. To limit computations, a :class:`QSscoreError` is thrown if we try
    more than the maximal number of chain mappings.
    The value must be set before the first use of :attr:`chain_mapping`.
    By default it is set to 100000.

    :type: :class:`int`

  .. attribute:: res_num_alignment

    Forces each alignment in :attr:`alignments` to be based on residue numbers
    instead of using a global BLOSUM62-based alignment.

    :type: :class:`bool`
  """
  def __init__(self, ent_1, ent_2, res_num_alignment=False):
    # generate QSscoreEntity objects?
    if isinstance(ent_1, QSscoreEntity):
      self.qs_ent_1 = ent_1
    else:
      self.qs_ent_1 = QSscoreEntity(ent_1)
    if isinstance(ent_2, QSscoreEntity):
      self.qs_ent_2 = ent_2
    else:
      self.qs_ent_2 = QSscoreEntity(ent_2)
    # check validity of inputs
    if not self.qs_ent_1.is_valid or not self.qs_ent_2.is_valid:
      raise QSscoreError("Invalid input in QSscorer!")
    # set names for structures
    if self.qs_ent_1.original_name == self.qs_ent_2.original_name:
      self.qs_ent_1.SetName('pdb_1')
      self.qs_ent_2.SetName('pdb_2')
    else:
      self.qs_ent_1.SetName(self.qs_ent_1.original_name)
      self.qs_ent_2.SetName(self.qs_ent_2.original_name)
    # set other public attributes
    self.res_num_alignment = res_num_alignment
    self.calpha_only = self.qs_ent_1.calpha_only or self.qs_ent_2.calpha_only
    self.max_ca_per_chain_for_cm = 100
    self.max_mappings_extensive = 100000
    # init cached stuff
    self._chem_mapping = None
    self._ent_to_cm_1 = None
    self._ent_to_cm_2 = None
    self._symm_1 = None
    self._symm_2 = None
    self._chain_mapping = None
    self._chain_mapping_scheme = None
    self._alignments = None
    self._mapped_residues = None
    self._global_score = None
    self._best_score = None
    self._superposition = None
    self._clustalw_bin = None

  @property
  def chem_mapping(self):
    """Inter-complex mapping of chemical groups.

    Each group (see :attr:`QSscoreEntity.chem_groups`) is mapped according to
    highest sequence identity. Alignment is between longest sequences in groups.

    Limitations:

    - If different numbers of groups, we map only the groups for the complex
      with less groups (rest considered unmapped and shown as warning)
    - The mapping is forced: the "best" mapping will be chosen independently of
      how low the seq. identity may be

    :getter: Computed on first use (cached)
    :type: :class:`dict` with key = :class:`tuple` of chain names in
           :attr:`qs_ent_1` and value = :class:`tuple` of chain names in
           :attr:`qs_ent_2`.

    :raises: :class:`QSscoreError` if we end up having no chains for either
             entity in the mapping (can happen if chains do not have CA atoms).
    """
    if self._chem_mapping is None:
      self._chem_mapping = _GetChemGroupsMapping(self.qs_ent_1, self.qs_ent_2)
    return self._chem_mapping

  @chem_mapping.setter
  def chem_mapping(self, chem_mapping):
    self._chem_mapping = chem_mapping

  @property
  def ent_to_cm_1(self):
    """Subset of :attr:`qs_ent_1` used to compute chain mapping and symmetries.

    Properties:

    - Includes only residues aligned according to :attr:`chem_mapping`
    - Includes only 1 CA atom per residue
    - Has at least 5 and at most :attr:`max_ca_per_chain_for_cm` atoms per chain
    - All chains of the same chemical group have the same number of atoms
      (also in :attr:`ent_to_cm_2` according to :attr:`chem_mapping`)
    - All chains appearing in :attr:`chem_mapping` appear in this entity
      (so the two can be safely used together)

    This entity might be transformed (i.e. all positions rotated/translated by
    same transformation matrix) if this can speed up computations. So do not
    assume fixed global positions (but relative distances will remain fixed).

    :getter: Computed on first use (cached)
    :type: :class:`~ost.mol.EntityHandle`

    :raises: :class:`QSscoreError` if any chain ends up having less than 5 res.
    """
    if self._ent_to_cm_1 is None:
      self._ComputeAlignedEntities()
    return self._ent_to_cm_1

  @ent_to_cm_1.setter
  def ent_to_cm_1(self, ent_to_cm_1):
    self._ent_to_cm_1 = ent_to_cm_1

  @property
  def ent_to_cm_2(self):
    """Subset of :attr:`qs_ent_1` used to compute chain mapping and symmetries
    (see :attr:`ent_to_cm_1` for details).
    """
    if self._ent_to_cm_2 is None:
      self._ComputeAlignedEntities()
    return self._ent_to_cm_2

  @ent_to_cm_2.setter
  def ent_to_cm_2(self, ent_to_cm_2):
    self._ent_to_cm_2 = ent_to_cm_2

  @property
  def symm_1(self):
    """Symmetry groups for :attr:`qs_ent_1` used to speed up chain mapping.

    This is a list of chain-lists where each chain-list can be used reconstruct
    the others via cyclic C or dihedral D symmetry. The first chain-list is used
    as a representative symmetry group. For heteromers, the group-members must
    contain all different seqres in oligomer.

    Example: symm. groups [(A,B,C), (D,E,F), (G,H,I)] means that there are
    symmetry transformations to get (D,E,F) and (G,H,I) from (A,B,C).

    Properties:

    - All symmetry group tuples have the same length (num. of chains)
    - All chains in :attr:`ent_to_cm_1` appear (w/o duplicates)
    - For heteros: symmetry group tuples have all different chem. groups
    - Trivial symmetry group = one tuple with all chains (used if inconsistent
      data provided or if no symmetry is found)
    - Either compatible to :attr:`symm_2` or trivial symmetry groups used.
      Compatibility requires same lengths of symmetry group tuples and it must
      be possible to get an overlap (80% of residues covered within 6 A of a
      (chem. mapped) chain) of all chains in representative symmetry groups by
      superposing one pair of chains.

    :getter: Computed on first use (cached)
    :type: :class:`list` of :class:`tuple` of :class:`str` (chain names)
    """
    if self._symm_1 is None:
      self._ComputeSymmetry()
    return self._symm_1

  @property
  def symm_2(self):
    """Symmetry groups for :attr:`qs_ent_2` (see :attr:`symm_1` for details)."""
    if self._symm_2 is None:
      self._ComputeSymmetry()
    return self._symm_2

  def SetSymmetries(self, symm_1, symm_2):
    """Set user-provided symmetry groups.

    These groups are restricted to chain names appearing in :attr:`ent_to_cm_1`
    and :attr:`ent_to_cm_2` respectively. They are only valid if they cover all
    chains and both *symm_1* and *symm_2* have same lengths of symmetry group
    tuples. Otherwise trivial symmetry group used (see :attr:`symm_1`).

    :param symm_1: Value to set for :attr:`symm_1`.
    :param symm_2: Value to set for :attr:`symm_2`.
    """
    # restrict chain names
    self._symm_1 = _CleanUserSymmetry(symm_1, self.ent_to_cm_1)
    self._symm_2 = _CleanUserSymmetry(symm_2, self.ent_to_cm_2)
    # check that we have reasonable symmetries set (fallback: all chains)
    if not _AreValidSymmetries(self._symm_1, self._symm_2):
      self._symm_1 = [tuple(ch.name for ch in self.ent_to_cm_1.chains)]
      self._symm_2 = [tuple(ch.name for ch in self.ent_to_cm_2.chains)]

  @property
  def chain_mapping(self):
    """Mapping from :attr:`ent_to_cm_1` to :attr:`ent_to_cm_2`.

    Properties:

    - Mapping is between chains of same chem. group (see :attr:`chem_mapping`)
    - Each chain can appear only once in mapping
    - All chains of complex with less chains are mapped
    - Symmetry (:attr:`symm_1`, :attr:`symm_2`) is taken into account

    Details on algorithms used to find mapping:

    - We try all pairs of chem. mapped chains within symmetry group and get
      superpose-transformation for them
    - First option: check for "sufficient overlap" of other chain-pairs

      - For each chain-pair defined above: apply superposition to full oligomer
        and map chains based on structural overlap
      - Structural overlap = X% of residues in second oligomer covered within Y
        Angstrom of a (chem. mapped) chain in first oligomer. We successively
        try (X,Y) = (80,4), (40,6) and (20,8) to be less and less strict in
        mapping (warning shown for most permissive one).
      - If multiple possible mappings are found, we choose the one which leads
        to the lowest multi-chain-RMSD given the superposition

    - Fallback option: try all mappings to find minimal multi-chain-RMSD
      (warning shown)

      - For each chain-pair defined above: apply superposition, try all (!)
        possible chain mappings (within symmetry group) and keep mapping with
        lowest multi-chain-RMSD
      - Repeat procedure above to resolve symmetry. Within the symmetry group we
        can use the chain mapping computed before and we just need to find which
        symmetry group in first oligomer maps to which in the second one. We
        again try all possible combinations...
      - Limitations:
        
        - Trying all possible mappings is a combinatorial nightmare (factorial).
          We throw an exception if too many combinations (e.g. octomer vs
          octomer with no usable symmetry)
        - The mapping is forced: the "best" mapping will be chosen independently
          of how badly they fit in terms of multi-chain-RMSD
        - As a result, such a forced mapping can lead to a large range of
          resulting QS scores. An extreme example was observed between 1on3.1
          and 3u9r.1, where :attr:`global_score` can range from 0.12 to 0.43
          for mappings with very similar multi-chain-RMSD.

    :getter: Computed on first use (cached)
    :type: :class:`dict` with key / value = :class:`str` (chain names, key
           for :attr:`ent_to_cm_1`, value for :attr:`ent_to_cm_2`)
    :raises: :class:`QSscoreError` if there are too many combinations to check
             to find a chain mapping (see :attr:`max_mappings_extensive`).
    """
    if self._chain_mapping is None:
      self._chain_mapping, self._chain_mapping_scheme = \
        _GetChainMapping(self.ent_to_cm_1, self.ent_to_cm_2, self.symm_1,
                         self.symm_2, self.chem_mapping,
                         self.max_mappings_extensive)
      LogInfo('Mapping found: %s' % str(self._chain_mapping))
    return self._chain_mapping

  @chain_mapping.setter
  def chain_mapping(self, chain_mapping):
    self._chain_mapping = chain_mapping

  @property
  def chain_mapping_scheme(self):
    """Mapping scheme used to get :attr:`chain_mapping`.

    Possible values:

    - 'strict': 80% overlap needed within 4 Angstrom (overlap based mapping).
    - 'tolerant': 40% overlap needed within 6 Angstrom (overlap based mapping).
    - 'permissive': 20% overlap needed within 8 Angstrom (overlap based
      mapping). It's best if you check mapping manually!
    - 'extensive': Extensive search used for mapping detection (fallback). This
      approach has known limitations and may be removed in future versions.
      Mapping should be checked manually!
    - 'user': :attr:`chain_mapping` was set by user before first use of this
      attribute.

    :getter: Computed with :attr:`chain_mapping` on first use (cached)
    :type: :class:`str`
    :raises: :class:`QSscoreError` as in :attr:`chain_mapping`.
    """
    if self._chain_mapping_scheme is None:
      # default: user provided
      self._chain_mapping_scheme = 'user'
      # get chain mapping and make sure internal variable is set
      # -> will not compute and only update _chain_mapping if user provided
      # -> will compute and overwrite _chain_mapping_scheme else
      self._chain_mapping = self.chain_mapping
    return self._chain_mapping_scheme

  @property
  def alignments(self):
    """List of successful sequence alignments using :attr:`chain_mapping`.

    There will be one alignment for each mapped chain and they are ordered by
    their chain names in :attr:`qs_ent_1`.

    The first sequence of each alignment belongs to :attr:`qs_ent_1` and the
    second one to :attr:`qs_ent_2`. The sequences are named according to the
    mapped chain names and have views attached into :attr:`QSscoreEntity.ent`
    of :attr:`qs_ent_1` and :attr:`qs_ent_2`.

    If :attr:`res_num_alignment` is False, each alignment is performed using a
    global BLOSUM62-based alignment. Otherwise, the positions in the alignment
    sequences are simply given by the residue number so that residues with
    matching numbers are aligned.

    :getter: Computed on first use (cached)
    :type: :class:`list` of :class:`~ost.seq.AlignmentHandle`
    """
    if self._alignments is None:
      self._alignments = _GetMappedAlignments(self.qs_ent_1.ent,
                                              self.qs_ent_2.ent,
                                              self.chain_mapping,
                                              self.res_num_alignment)
    return self._alignments

  @alignments.setter
  def alignments(self, alignments):
    self._alignments = alignments

  @property
  def mapped_residues(self):
    """Mapping of shared residues in :attr:`alignments`.

    :getter: Computed on first use (cached)
    :type: :class:`dict` *mapped_residues[c1][r1] = r2* with:
           *c1* = Chain name in first entity (= first sequence in aln),
           *r1* = Residue number in first entity,
           *r2* = Residue number in second entity
    """
    if self._mapped_residues is None:
      self._mapped_residues = _GetMappedResidues(self.alignments)
    return self._mapped_residues

  @mapped_residues.setter
  def mapped_residues(self, mapped_residues):
    self._mapped_residues = mapped_residues

  @property
  def global_score(self):
    """QS-score with penalties.
    
    The range of the score is between 0 (i.e. no interface residues are shared
    between biounits) and 1 (i.e. the interfaces are identical).
    
    The global QS-score is computed applying penalties when interface residues
    or entire chains are missing (i.e. anything that is not mapped in
    :attr:`mapped_residues` / :attr:`chain_mapping`) in one of the biounits.

    :getter: Computed on first use (cached)
    :type: :class:`float`
    :raises: :class:`QSscoreError` if only one chain is mapped
    """
    if self._global_score is None:
      self._ComputeScores()
    return self._global_score

  @property
  def best_score(self):
    """QS-score without penalties.

    Like :attr:`global_score`, but neglecting additional residues or chains in
    one of the biounits (i.e. the score is calculated considering only mapped
    chains and residues).

    :getter: Computed on first use (cached)
    :type: :class:`float`
    :raises: :class:`QSscoreError` if only one chain is mapped
    """
    if self._best_score is None:
      self._ComputeScores()
    return self._best_score

  @property
  def superposition(self):
    """Superposition result based on shared CA atoms in :attr:`alignments`.

    The superposition can be used to map :attr:`QSscoreEntity.ent` of
    :attr:`qs_ent_1` onto the one of :attr:`qs_ent_2`. Use
    :func:`ost.geom.Invert` if you need the opposite transformation.

    :getter: Computed on first use (cached)
    :type: :class:`ost.mol.alg.SuperpositionResult`
    """
    if self._superposition is None:
      self._superposition = _GetQsSuperposition(self.alignments)
      # report it
      sup_rmsd = self._superposition.rmsd
      cmp_view = self._superposition.view1
      LogInfo('CA RMSD for %s aligned residues on %s chains: %.2f' \
              % (cmp_view.residue_count, cmp_view.chain_count, sup_rmsd))
    return self._superposition

  @property
  def clustalw_bin(self):
    """
    Full path to ``clustalw`` or ``clustalw2`` executable to use for multiple
    sequence alignments (unless :attr:`chain_mapping` is provided manually).

    :getter: Located in path on first use (cached)
    :type: :class:`str`
    """
    if self._clustalw_bin is None:
      self._clustalw_bin = settings.Locate(('clustalw', 'clustalw2'))
    return self._clustalw_bin

  @clustalw_bin.setter
  def clustalw_bin(self, clustalw_bin):
    self._clustalw_bin = clustalw_bin

  def GetOligoLDDTScorer(self, settings, penalize_extra_chains=True):
    """
    :return: :class:`OligoLDDTScorer` object, setup for this QS scoring problem.
             The scorer is set up with :attr:`qs_ent_1` as the reference and
             :attr:`qs_ent_2` as the model.
    :param settings: Passed to :class:`OligoLDDTScorer` constructor.
    :param penalize_extra_chains: Passed to :class:`OligoLDDTScorer` constructor.
    """
    if penalize_extra_chains:
      return OligoLDDTScorer(self.qs_ent_1.ent, self.qs_ent_2.ent,
                             self.alignments, self.calpha_only, settings,
                             True, self.chem_mapping)
    else:
      return OligoLDDTScorer(self.qs_ent_1.ent, self.qs_ent_2.ent,
                             self.alignments, self.calpha_only, settings, False)


  ##############################################################################
  # Class internal helpers (anything that doesnt easily work without this class)
  ##############################################################################

  def _ComputeAlignedEntities(self):
    """Fills cached ent_to_cm_1 and ent_to_cm_2."""
    # get aligned residues via MSA
    ev1, ev2 = _GetAlignedResidues(self.qs_ent_1, self.qs_ent_2,
                                   self.chem_mapping,
                                   self.max_ca_per_chain_for_cm,
                                   self.clustalw_bin)
    # extract new entities
    self._ent_to_cm_1 = mol.CreateEntityFromView(ev1, True)
    self._ent_to_cm_2 = mol.CreateEntityFromView(ev2, True)
    # name them
    self._ent_to_cm_1.SetName(self.qs_ent_1.GetName())
    self._ent_to_cm_2.SetName(self.qs_ent_2.GetName())

  def _ComputeSymmetry(self):
    """Fills cached symm_1 and symm_2."""
    # get them
    self._symm_1, self._symm_2 = \
      _FindSymmetry(self.qs_ent_1, self.qs_ent_2, self.ent_to_cm_1,
                    self.ent_to_cm_2, self.chem_mapping)
    # check that we have reasonable symmetries set (fallback: all chains)
    if not _AreValidSymmetries(self._symm_1, self._symm_2):
      self._symm_1 = [tuple(ch.name for ch in self.ent_to_cm_1.chains)]
      self._symm_2 = [tuple(ch.name for ch in self.ent_to_cm_2.chains)]

  def _ComputeScores(self):
    """Fills cached global_score and best_score."""
    if len(self.chain_mapping) < 2:
      raise QSscoreError("QS-score is not defined for monomers")
    # get contacts
    if self.calpha_only:
      contacts_1 = self.qs_ent_1.contacts_ca
      contacts_2 = self.qs_ent_2.contacts_ca
    else:
      contacts_1 = self.qs_ent_1.contacts
      contacts_2 = self.qs_ent_2.contacts
    # get scores
    scores = _GetScores(contacts_1, contacts_2, self.mapped_residues,
                        self.chain_mapping)
    self._best_score = scores[0]
    self._global_score = scores[1]
    # report scores
    LogInfo('QSscore %s, %s: best: %.2f, global: %.2f' \
            % (self.qs_ent_1.GetName(), self.qs_ent_2.GetName(),
               self._best_score, self._global_score))


###############################################################################
# Entity with cached entries for QS scoring
###############################################################################

class QSscoreEntity(object):
  """Entity with cached entries for QS scoring.

  Any known / precomputed information can be filled into the appropriate
  attribute here as long as they are labelled as read/write. Otherwise the
  quantities are computed on first access and cached (lazy evaluation). The
  heaviest load is expected when computing :attr:`contacts` and
  :attr:`contacts_ca`.

  :param ent: Entity to be used for QS scoring. A copy of it will be processed.
  :type ent:  :class:`~ost.mol.EntityHandle` or :class:`~ost.mol.EntityView`

  .. attribute:: is_valid

    True, if successfully initialized. False, if input structure has no protein
    chains with >= 20 residues.

    :type: :class:`bool`

  .. attribute:: original_name

    Name set for *ent* when object was created.

    :type: :class:`str`

  .. attribute:: ent

    Cleaned version of *ent* passed at construction. Hydrogens are removed, the
    entity is processed with a :class:`~ost.conop.RuleBasedProcessor` and chains
    listed in :attr:`removed_chains` have been removed. The name of this entity
    might change during scoring (see :func:`GetName`). Otherwise, this will be
    fixed.

    :type: :class:`~ost.mol.EntityHandle`

  .. attribute:: removed_chains

    Chains removed from *ent* passed at construction. These are ligand and water
    chains as well as small (< 20 res.) peptides or chains with no amino acids
    (determined by chem. type, which is set by rule based processor).

    :type: :class:`list` of :class:`str`

  .. attribute:: calpha_only

    Whether entity is CA-only (i.e. it has 0 CB atoms)

    :type: :class:`bool`
  """
  def __init__(self, ent):
    # copy entity and process/clean it
    self.original_name = ent.GetName()
    ent = mol.CreateEntityFromView(ent.Select('ele!=H and aname!=HN'), True)
    if not conop.GetDefaultLib():
      raise RuntimeError("QSscore computation requires a compound library!")
    pr = conop.RuleBasedProcessor(conop.GetDefaultLib())
    pr.Process(ent)
    self.ent, self.removed_chains, self.calpha_only = _CleanInputEntity(ent)
    # check if it's suitable for QS scoring
    if self.ent.chain_count == 0:
      LogError('Bad input file: ' + ent.GetName() + '. No chains left after '
               'removing water, ligands and small peptides.')
      self.is_valid = False
    elif self.ent.chain_count == 1:
      LogWarning('Structure ' + ent.GetName() + ' is a monomer.')
      self.is_valid = True
    else:
      self.is_valid = True
    # init cached stuff
    self._ca_entity = None
    self._ca_chains = None
    self._alignments = {}
    self._chem_groups = None
    self._angles = {}
    self._axis = {}
    self._contacts = None
    self._contacts_ca = None

  def GetName(self):
    """Wrapper to :func:`~ost.mol.EntityHandle.GetName` of :attr:`ent`.
    This is used to uniquely identify the entity while scoring. The name may
    therefore change while :attr:`original_name` remains fixed.
    """
    # for duck-typing and convenience
    return self.ent.GetName()

  def SetName(self, new_name):
    """Wrapper to :func:`~ost.mol.EntityHandle.SetName` of :attr:`ent`.
    Use this to change unique identifier while scoring (see :func:`GetName`).
    """
    # for duck-typing and convenience
    self.ent.SetName(new_name)

  @property
  def ca_entity(self):
    """
    Reduced representation of :attr:`ent` with only CA atoms.
    This guarantees that each included residue has exactly one atom.

    :getter: Computed on first use (cached)
    :type: :class:`~ost.mol.EntityHandle`
    """
    if self._ca_entity is None:
      self._ca_entity = _GetCAOnlyEntity(self.ent)
    return self._ca_entity

  @property
  def ca_chains(self):
    """
    Map of chain names in :attr:`ent` to sequences with attached view to CA-only
    chains (into :attr:`ca_entity`). Useful for alignments and superpositions.

    :getter: Computed on first use (cached)
    :type: :class:`dict` (key = :class:`str`,
           value = :class:`~ost.seq.SequenceHandle`)
    """
    if self._ca_chains is None:
      self._ca_chains = dict()
      ca_entity = self.ca_entity
      for ch in ca_entity.chains:
        self._ca_chains[ch.name] = seq.SequenceFromChain(ch.name, ch)
    return self._ca_chains

  def GetAlignment(self, c1, c2):
    """Get sequence alignment of chain *c1* with chain *c2*.
    Computed on first use based on :attr:`ca_chains` (cached).

    :param c1: Chain name for first chain to align.
    :type c1:  :class:`str`
    :param c2: Chain name for second chain to align.
    :type c2:  :class:`str`
    :rtype: :class:`~ost.seq.AlignmentHandle` or None if it failed.
    """
    if (c1,c2) not in self._alignments:
      ca_chains = self.ca_chains
      self._alignments[(c1,c2)] = _AlignAtomSeqs(ca_chains[c1], ca_chains[c2])
    return self._alignments[(c1,c2)]

  @property
  def chem_groups(self):
    """
    Intra-complex group of chemically identical (seq. id. > 95%) polypeptide
    chains as extracted from :attr:`ca_chains`. First chain in group is the one
    with the longest sequence.

    :getter: Computed on first use (cached)
    :type: :class:`list` of :class:`list` of :class:`str` (chain names)
    """
    if self._chem_groups is None:
      self._chem_groups = _GetChemGroups(self, 95)
      LogInfo('Chemically equivalent chain-groups in %s: %s' \
              % (self.GetName(), str(self._chem_groups)))
    return self._chem_groups

  def GetAngles(self, c1, c2):
    """Get Euler angles from superposition of chain *c1* with chain *c2*.
    Computed on first use based on :attr:`ca_chains` (cached).

    :param c1: Chain name for first chain to superpose.
    :type c1:  :class:`str`
    :param c2: Chain name for second chain to superpose.
    :type c2:  :class:`str`
    :return: 3 Euler angles (may contain nan if something fails).
    :rtype:  :class:`numpy.array`
    """
    if (c1,c2) not in self._angles:
      self._GetSuperposeData(c1, c2)
    return self._angles[(c1,c2)]

  def GetAxis(self, c1, c2):
    """Get axis of symmetry from superposition of chain *c1* with chain *c2*.
    Computed on first use based on :attr:`ca_chains` (cached).

    :param c1: Chain name for first chain to superpose.
    :type c1:  :class:`str`
    :param c2: Chain name for second chain to superpose.
    :type c2:  :class:`str`
    :return: Rotational axis (may contain nan if something fails).
    :rtype:  :class:`numpy.array`
    """
    if (c1,c2) not in self._axis:
      self._GetSuperposeData(c1, c2)
    return self._axis[(c1,c2)]

  @property
  def contacts(self):
    """
    Connectivity dictionary (**read/write**).
    As given by :func:`GetContacts` with *calpha_only* = False on :attr:`ent`.

    :getter: Computed on first use (cached)
    :setter: Uses :func:`FilterContacts` to ensure that we only keep contacts
             for chains in the cleaned entity.
    :type: See return type of :func:`GetContacts`
    """
    if self._contacts is None:
      self._contacts = GetContacts(self.ent, False)
    return self._contacts

  @contacts.setter
  def contacts(self, new_contacts):
    chain_names = set([ch.name for ch in self.ent.chains])
    self._contacts = FilterContacts(new_contacts, chain_names)
  
  @property
  def contacts_ca(self):
    """
    CA-only connectivity dictionary (**read/write**).
    Like :attr:`contacts` but with *calpha_only* = True in :func:`GetContacts`.
    """
    if self._contacts_ca is None:
      self._contacts_ca = GetContacts(self.ent, True)
    return self._contacts_ca
  
  @contacts_ca.setter
  def contacts_ca(self, new_contacts):
    chain_names = set([ch.name for ch in self.ent.chains])
    self._contacts_ca = FilterContacts(new_contacts, chain_names)

  ##############################################################################
  # Class internal helpers (anything that doesnt easily work without this class)
  ##############################################################################

  def _GetSuperposeData(self, c1, c2):
    """Fill _angles and _axis from superposition of CA chains of c1 and c2."""
    # get aligned views (must contain identical numbers of atoms!)
    aln = self.GetAlignment(c1, c2)
    if not aln:
      # fallback for non-aligned stuff (nan)
      self._angles[(c1,c2)] = np.empty(3) * np.nan
      self._axis[(c1,c2)] = np.empty(3) * np.nan
      return
    v1, v2 = seq.ViewsFromAlignment(aln)
    if v1.atom_count < 3:
      # fallback for non-aligned stuff (nan)
      self._angles[(c1,c2)] = np.empty(3) * np.nan
      self._axis[(c1,c2)] = np.empty(3) * np.nan
      return
    # get transformation
    sup_res = mol.alg.SuperposeSVD(v1, v2, apply_transform=False)
    Rt = sup_res.transformation
    # extract angles
    a,b,c = _GetAngles(Rt)
    self._angles[(c1,c2)] = np.asarray([a,b,c])
    # extract axis of symmetry
    R3 = geom.Rotation3(Rt.ExtractRotation())
    self._axis[(c1,c2)] = np.asarray(R3.GetRotationAxis().data)

###############################################################################
# Contacts computations
###############################################################################

def FilterContacts(contacts, chain_names):
  """Filter contacts to contain only contacts for chains in *chain_names*.

  :param contacts: Connectivity dictionary as produced by :func:`GetContacts`.
  :type contacts:  :class:`dict`
  :param chain_names: Chain names to keep.
  :type chain_names:  :class:`list` or (better) :class:`set`
  :return: New connectivity dictionary (format as in :func:`GetContacts`)
  :rtype:  :class:`dict`
  """
  # create new dict with subset
  filtered_contacts = dict()
  for c1 in contacts:
    if c1 in chain_names:
      new_contacts = dict()
      for c2 in contacts[c1]:
        if c2 in chain_names:
          new_contacts[c2] = contacts[c1][c2]
      # avoid adding empty dicts
      if new_contacts:
        filtered_contacts[c1] = new_contacts
  return filtered_contacts

def GetContacts(entity, calpha_only, dist_thr=12.0):
  """Get inter-chain contacts of a macromolecular entity.

  Contacts are pairs of residues within a given distance belonging to different
  chains. They are stored once per pair and include the CA/CB-CA/CB distance.

  :param entity: An entity to check connectivity for.
  :type entity:  :class:`~ost.mol.EntityHandle` or :class:`~ost.mol.EntityView`
  :param calpha_only: If True, we only consider CA-CA distances. Else, we use CB
                      unless the residue is a GLY.
  :type calpha_only:  :class:`bool`
  :param dist_thr: Maximal CA/CB-CA/CB distance to be considered in contact.
  :type dist_thr:  :class:`float`
  :return: A connectivity dictionary. A pair of residues with chain names
           *ch_name1* & *ch_name2* (*ch_name1* < *ch_name2*), residue numbers
           *res_num1* & *res_num2* and distance *dist* (<= *dist_thr*) are
           stored as *result[ch_name1][ch_name2][res_num1][res_num2]* = *dist*.
  :rtype:  :class:`dict`
  """
  # get ent copy to search on
  if calpha_only:
    ev = entity.Select("aname=CA")
  else:
    ev = entity.Select("(rname=GLY and aname=CA) or aname=CB")
  ent = mol.CreateEntityFromView(ev, True)
  # search all vs all
  contacts = dict()
  for atom in ent.atoms:
    ch_name1 = atom.chain.name
    res_num1 = atom.residue.number.num
    close_atoms = ent.FindWithin(atom.pos, dist_thr)
    for close_atom in close_atoms:
      ch_name2 = close_atom.chain.name
      if ch_name2 > ch_name1:
        res_num2 = close_atom.residue.number.num
        dist = geom.Distance(atom.pos, close_atom.pos)
        # add to contacts
        if ch_name1 not in contacts:
          contacts[ch_name1] = dict()
        if ch_name2 not in contacts[ch_name1]:
          contacts[ch_name1][ch_name2] = dict()
        if res_num1 not in contacts[ch_name1][ch_name2]:
          contacts[ch_name1][ch_name2][res_num1] = dict()
        contacts[ch_name1][ch_name2][res_num1][res_num2] = round(dist, 3)
  # DONE
  return contacts

###############################################################################
# Oligo-lDDT scores
###############################################################################

class OligoLDDTScorer(object):
  """Helper class to calculate oligomeric lDDT scores.

  This class can be used independently, but commonly it will be created by
  calling :func:`QSscorer.GetOligoLDDTScorer`.

  .. note::

    By construction, lDDT scores are not symmetric and hence it matters which
    structure is the reference (:attr:`ref`) and which one is the model
    (:attr:`mdl`). Extra residues in the model are generally not considered.
    Extra chains in both model and reference can be considered by setting the
    :attr:`penalize_extra_chains` flag to True.

  :param ref: Sets :attr:`ref`
  :param mdl: Sets :attr:`mdl`
  :param alignments: Sets :attr:`alignments`
  :param calpha_only: Sets :attr:`calpha_only`
  :param settings: Sets :attr:`settings`
  :param penalize_extra_chains: Sets :attr:`penalize_extra_chains`
  :param chem_mapping: Sets :attr:`chem_mapping`. Must be given if
                       *penalize_extra_chains* is True.
  
  .. attribute:: ref
                 mdl

    Full reference/model entity to be scored. The entity must contain all chains
    mapped in :attr:`alignments` and may also contain additional ones which are
    considered if :attr:`penalize_extra_chains` is True.

    :type: :class:`~ost.mol.EntityHandle`
  
  .. attribute:: alignments

    One alignment for each mapped chain of :attr:`ref`/:attr:`mdl` as defined in
    :attr:`QSscorer.alignments`. The first sequence of each alignment belongs to
    :attr:`ref` and the second one to :attr:`mdl`. Sequences must have sequence
    naming and attached views as defined in :attr:`QSscorer.alignments`.

    :type: :class:`list` of :class:`~ost.seq.AlignmentHandle`

  .. attribute:: calpha_only

    If True, restricts lDDT score to CA only.

    :type: :class:`bool`

  .. attribute:: settings

    Settings to use for lDDT scoring.

    :type: :class:`~ost.mol.alg.lDDTSettings`

  .. attribute:: penalize_extra_chains

    If True, extra chains in both :attr:`ref` and :attr:`mdl` will penalize the
    lDDT scores.

    :type: :class:`bool`

  .. attribute:: chem_mapping

    Inter-complex mapping of chemical groups as defined in
    :attr:`QSscorer.chem_mapping`. Used to find "chem-mapped" chains in
    :attr:`ref` for unmapped chains in :attr:`mdl` when penalizing scores.
    Each unmapped model chain can add extra reference-contacts according to the
    average total contacts of each single "chem-mapped" reference chain. If
    there is no "chem-mapped" reference chain, a warning is shown and the model
    chain is ignored.


    Only relevant if :attr:`penalize_extra_chains` is True.

    :type: :class:`dict` with key = :class:`tuple` of chain names in
           :attr:`ref` and value = :class:`tuple` of chain names in
           :attr:`mdl`.
  """

  # NOTE: one could also allow computation of both penalized and unpenalized
  #       in same object -> must regenerate lddt_ref / lddt_mdl though

  def __init__(self, ref, mdl, alignments, calpha_only, settings,
               penalize_extra_chains=False, chem_mapping=None):
    # sanity checks
    if chem_mapping is None and penalize_extra_chains:
      raise RuntimeError("Must provide chem_mapping when requesting penalty "
                         "for extra chains!")
    if not penalize_extra_chains:
      # warn for unmapped model chains
      unmapped_mdl_chains = self._GetUnmappedMdlChains(mdl, alignments)
      if unmapped_mdl_chains:
        LogWarning('MODEL contains chains unmapped to REFERENCE, '
                   'lDDT is not considering MODEL chains %s' \
                   % str(list(unmapped_mdl_chains)))
      # warn for unmapped reference chains
      ref_chains = set(ch.name for ch in ref.chains)
      mapped_ref_chains = set(aln.GetSequence(0).GetName() for aln in alignments)
      unmapped_ref_chains = (ref_chains - mapped_ref_chains)
      if unmapped_ref_chains:
        LogWarning('REFERENCE contains chains unmapped to MODEL, '
                   'lDDT is not considering REFERENCE chains %s' \
                   % str(list(unmapped_ref_chains)))
    # prepare fields
    self.ref = ref
    self.mdl = mdl
    self.alignments = alignments
    self.calpha_only = calpha_only
    self.settings = settings
    self.penalize_extra_chains = penalize_extra_chains
    self.chem_mapping = chem_mapping
    self._sc_lddt = None
    self._oligo_lddt = None
    self._weighted_lddt = None
    self._lddt_ref = None
    self._lddt_mdl = None
    self._oligo_lddt_scorer = None
    self._mapped_lddt_scorers = None
    self._ref_scorers = None
    self._model_penalty = None

  @property
  def oligo_lddt(self):
    """Oligomeric lDDT score.

    The score is computed as conserved contacts divided by the total contacts
    in the reference using the :attr:`oligo_lddt_scorer`, which uses the full
    complex as reference/model structure. If :attr:`penalize_extra_chains` is
    True, the reference/model complexes contain all chains (otherwise only the
    mapped ones) and additional contacts are added to the reference's total
    contacts for unmapped model chains according to the :attr:`chem_mapping`.

    The main difference with :attr:`weighted_lddt` is that the lDDT scorer
    "sees" the full complex here (incl. inter-chain contacts), while the
    weighted single chain score looks at each chain separately.

    :getter: Computed on first use (cached)
    :type: :class:`float`
    """
    if self._oligo_lddt is None:
      LogInfo('Reference %s has: %s chains' \
              % (self.ref.GetName(), self.ref.chain_count))
      LogInfo('Model %s has: %s chains' \
              % (self.mdl.GetName(), self.mdl.chain_count))

      # score with or w/o extra-chain penalty
      if self.penalize_extra_chains:
        denominator = self.oligo_lddt_scorer.total_contacts
        denominator += self._GetModelPenalty()
        if denominator > 0:
          oligo_lddt = self.oligo_lddt_scorer.conserved_contacts \
                     / float(denominator)
        else:
          oligo_lddt = 0.0
      else:
        oligo_lddt = self.oligo_lddt_scorer.global_score
      self._oligo_lddt = oligo_lddt
    return self._oligo_lddt

  @property
  def weighted_lddt(self):
    """Weighted average of single chain lDDT scores.

    The score is computed as a weighted average of single chain lDDT scores
    (see :attr:`sc_lddt_scorers`) using the total contacts of each single
    reference chain as weights. If :attr:`penalize_extra_chains` is True,
    unmapped chains are added with a 0 score and total contacts taken from
    the actual reference chains or (for unmapped model chains) using the
    :attr:`chem_mapping`.

    See :attr:`oligo_lddt` for a comparison of the two scores.

    :getter: Computed on first use (cached)
    :type: :class:`float`
    """
    if self._weighted_lddt is None:
      scores = [s.global_score for s in self.sc_lddt_scorers]
      weights = [s.total_contacts for s in self.sc_lddt_scorers]
      nominator = sum([s * w for s, w in zip(scores, weights)])
      if self.penalize_extra_chains:
        ref_scorers = self._GetRefScorers()
        denominator = sum(s.total_contacts for s in list(ref_scorers.values()))
        denominator += self._GetModelPenalty()
      else:
        denominator = sum(weights)
      if denominator > 0:
        self._weighted_lddt = nominator / float(denominator)
      else:
        self._weighted_lddt = 0.0
    return self._weighted_lddt

  @property
  def lddt_ref(self):
    """The reference entity used for oligomeric lDDT scoring
    (:attr:`oligo_lddt` / :attr:`oligo_lddt_scorer`).
    
    Since the lDDT computation requires a single chain with mapped residue
    numbering, all chains of :attr:`ref` are appended into a single chain X with
    unique residue numbers according to the column-index in the alignment. The
    alignments are in the same order as they appear in :attr:`alignments`.
    Additional residues are appended at the end of the chain with unique residue
    numbers. Unmapped chains are only added if :attr:`penalize_extra_chains` is
    True. Only CA atoms are considered if :attr:`calpha_only` is True.

    :getter: Computed on first use (cached)
    :type: :class:`~ost.mol.EntityHandle`
    """
    if self._lddt_ref is None:
      self._PrepareOligoEntities()
    return self._lddt_ref

  @property
  def lddt_mdl(self):
    """The model entity used for oligomeric lDDT scoring
    (:attr:`oligo_lddt` / :attr:`oligo_lddt_scorer`).

    Like :attr:`lddt_ref`, this is a single chain X containing all chains of
    :attr:`mdl`. The residue numbers match the ones in :attr:`lddt_ref` where
    aligned and have unique numbers for additional residues.

    :getter: Computed on first use (cached)
    :type: :class:`~ost.mol.EntityHandle`
    """
    if self._lddt_mdl is None:
      self._PrepareOligoEntities()
    return self._lddt_mdl

  @property
  def oligo_lddt_scorer(self):
    """lDDT Scorer object for :attr:`lddt_ref` and :attr:`lddt_mdl`.

    :getter: Computed on first use (cached)
    :type: :class:`~ost.mol.alg.lDDTScorer`
    """
    if self._oligo_lddt_scorer is None:
      self._oligo_lddt_scorer = lDDTScorer(
        references=[self.lddt_ref.Select("")],
        model=self.lddt_mdl.Select(""),
        settings=self.settings)
    return self._oligo_lddt_scorer

  @property
  def mapped_lddt_scorers(self):
    """List of scorer objects for each chain mapped in :attr:`alignments`.

    :getter: Computed on first use (cached)
    :type: :class:`list` of :class:`MappedLDDTScorer`
    """
    if self._mapped_lddt_scorers is None:
      self._mapped_lddt_scorers = list()
      for aln in self.alignments:
        mapped_lddt_scorer = MappedLDDTScorer(aln, self.calpha_only,
                                              self.settings)
        self._mapped_lddt_scorers.append(mapped_lddt_scorer)
    return self._mapped_lddt_scorers

  @property
  def sc_lddt_scorers(self):
    """List of lDDT scorer objects extracted from :attr:`mapped_lddt_scorers`.

    :type: :class:`list` of :class:`~ost.mol.alg.lDDTScorer`
    """
    return [mls.lddt_scorer for mls in self.mapped_lddt_scorers]

  @property
  def sc_lddt(self):
    """List of global scores extracted from :attr:`sc_lddt_scorers`.

    If scoring for a mapped chain fails, an error is displayed and a score of 0
    is assigned.

    :getter: Computed on first use (cached)
    :type: :class:`list` of :class:`float`
    """
    if self._sc_lddt is None:
      self._sc_lddt = list()
      for lddt_scorer in self.sc_lddt_scorers:
        try:
          self._sc_lddt.append(lddt_scorer.global_score)
        except Exception as ex:
          LogError('Single chain lDDT failed:', str(ex))
          self._sc_lddt.append(0.0)
    return self._sc_lddt

  ##############################################################################
  # Class internal helpers
  ##############################################################################

  def _PrepareOligoEntities(self):
    # simple wrapper to avoid code duplication
    self._lddt_ref, self._lddt_mdl = _MergeAlignedChains(
      self.alignments, self.ref, self.mdl, self.calpha_only,
      self.penalize_extra_chains)

  @staticmethod
  def _GetUnmappedMdlChains(mdl, alignments):
    # assume model is second sequence in alignment and is named by chain
    mdl_chains = set(ch.name for ch in mdl.chains)
    mapped_mdl_chains = set(aln.GetSequence(1).GetName() for aln in alignments)
    return (mdl_chains - mapped_mdl_chains)

  def _GetRefScorers(self):
    # single chain lddt scorers for each reference chain (key = chain name)
    if self._ref_scorers is None:
      # collect from mapped_lddt_scorers
      ref_scorers = dict()
      for mapped_lddt_scorer in self.mapped_lddt_scorers:
        ref_ch_name = mapped_lddt_scorer.reference_chain_name
        ref_scorers[ref_ch_name] = mapped_lddt_scorer.lddt_scorer
      # add new ones where needed
      for ch in self.ref.chains:
        if ch.name not in ref_scorers:
          if self.calpha_only:
            ref_chain = ch.Select('aname=CA')
          else:
            ref_chain = ch.Select('')
          ref_scorers[ch.name] = lDDTScorer(
            references=[ref_chain],
            model=ref_chain,
            settings=self.settings)
      # store in cache
      self._ref_scorers = ref_scorers
    # fetch from cache
    return self._ref_scorers

  def _GetModelPenalty(self):
    # extra value to add to total number of distances for extra model chains
    # -> estimated from chem-mapped reference chains
    if self._model_penalty is None:
      # sanity check
      if self.chem_mapping is None:
        raise RuntimeError("Must provide chem_mapping when requesting penalty "
                           "for extra model chains!")
      # get cached ref_scorers
      ref_scorers = self._GetRefScorers()
      # get unmapped model chains
      unmapped_mdl_chains = self._GetUnmappedMdlChains(self.mdl, self.alignments)
      # map extra chains to ref. chains
      model_penalty = 0
      for ch_name_mdl in sorted(unmapped_mdl_chains):
        # get penalty for chain
        cur_penalty = None
        for cm_ref, cm_mdl in self.chem_mapping.items():
          if ch_name_mdl in cm_mdl:
            # penalize by an average of the chem. mapped ref. chains
            cur_penalty = 0
            for ch_name_ref in cm_ref:
              # assumes that total_contacts is cached (for speed)
              cur_penalty += ref_scorers[ch_name_ref].total_contacts
            cur_penalty /= float(len(cm_ref))
            break
        # report penalty
        if cur_penalty is None:
          LogWarning('Extra MODEL chain %s could not be chemically mapped to '
                     'any chain in REFERENCE, lDDT cannot consider it!' \
                     % ch_name_mdl)
        else:
          LogScript('Extra MODEL chain %s added to lDDT score by considering '
                    'chemically mapped chains in REFERENCE.' % ch_name_mdl)
          model_penalty += cur_penalty
      # store in cache
      self._model_penalty = model_penalty
    # fetch from cache
    return self._model_penalty


class MappedLDDTScorer(object):
  """A simple class to calculate a single-chain lDDT score on a given chain to
  chain mapping as extracted from :class:`OligoLDDTScorer`.

  :param alignment: Sets :attr:`alignment`
  :param calpha_only: Sets :attr:`calpha_only`
  :param settings: Sets :attr:`settings`
  
  .. attribute:: alignment

    Alignment with two sequences named according to the mapped chains and with
    views attached to both sequences (e.g. one of the items of
    :attr:`QSscorer.alignments`).

    The first sequence is assumed to be the reference and the second one the
    model. Since the lDDT score is not symmetric (extra residues in model are
    ignored), the order is important.

    :type: :class:`~ost.seq.AlignmentHandle`
  
  .. attribute:: calpha_only

    If True, restricts lDDT score to CA only.

    :type: :class:`bool`

  .. attribute:: settings

    Settings to use for lDDT scoring.

    :type: :class:`~ost.mol.alg.lDDTSettings`

  .. attribute:: lddt_scorer

    lDDT Scorer object for the given chains.

    :type: :class:`~ost.mol.alg.lDDTScorer`

  .. attribute:: reference_chain_name

    Chain name of the reference.

    :type: :class:`str`

  .. attribute:: model_chain_name

    Chain name of the model.

    :type: :class:`str`
  """
  def __init__(self, alignment, calpha_only, settings):
    # prepare fields
    self.alignment = alignment
    self.calpha_only = calpha_only
    self.settings = settings
    self.lddt_scorer = None # set in _InitScorer
    self.reference_chain_name = alignment.sequences[0].name
    self.model_chain_name = alignment.sequences[1].name
    self._old_number_label = "old_num"
    self._extended_alignment = None  # set in _InitScorer
    # initialize lDDT scorer
    self._InitScorer()

  def GetPerResidueScores(self):
    """
    :return: Scores for each residue
    :rtype:  :class:`list` of :class:`dict` with one item for each residue
             existing in model and reference:

             - "residue_number": Residue number in reference chain
             - "residue_name": Residue name in reference chain
             - "lddt": local lDDT
             - "conserved_contacts": number of conserved contacts
             - "total_contacts": total number of contacts
    """
    scores = list()
    assigned_residues = list()
    # Make sure the score is calculated
    self.lddt_scorer.global_score
    for col in self._extended_alignment:
      if col[0] != "-" and col.GetResidue(3).IsValid():
        ref_res = col.GetResidue(0)
        mdl_res = col.GetResidue(1)
        ref_res_renum = col.GetResidue(2)
        mdl_res_renum = col.GetResidue(3)
        if ref_res.one_letter_code != ref_res_renum.one_letter_code:
          raise RuntimeError("Reference residue name mapping inconsistent: %s != %s" %
                             (ref_res.one_letter_code,
                              ref_res_renum.one_letter_code))
        if mdl_res.one_letter_code != mdl_res_renum.one_letter_code:
          raise RuntimeError("Model residue name mapping inconsistent: %s != %s" %
                             (mdl_res.one_letter_code,
                              mdl_res_renum.one_letter_code))
        if ref_res.GetNumber().num != ref_res_renum.GetIntProp(self._old_number_label):
          raise RuntimeError("Reference residue number mapping inconsistent: %s != %s" %
                             (ref_res.GetNumber().num,
                              ref_res_renum.GetIntProp(self._old_number_label)))
        if mdl_res.GetNumber().num != mdl_res_renum.GetIntProp(self._old_number_label):
          raise RuntimeError("Model residue number mapping inconsistent: %s != %s" %
                             (mdl_res.GetNumber().num,
                              mdl_res_renum.GetIntProp(self._old_number_label)))
        if ref_res.qualified_name in assigned_residues:
          raise RuntimeError("Duplicated residue in reference: " %
                             (ref_res.qualified_name))
        else:
          assigned_residues.append(ref_res.qualified_name)
        # check if property there (may be missing for CA-only)
        if mdl_res_renum.HasProp(self.settings.label):
          scores.append({
            "residue_number": ref_res.GetNumber().num,
            "residue_name": ref_res.name,
            "lddt": mdl_res_renum.GetFloatProp(self.settings.label),
            "conserved_contacts": mdl_res_renum.GetFloatProp(self.settings.label + "_conserved"),
            "total_contacts": mdl_res_renum.GetFloatProp(self.settings.label + "_total")})
    return scores

  ##############################################################################
  # Class internal helpers (anything that doesnt easily work without this class)
  ##############################################################################

  def _InitScorer(self):
    # Use copy of alignment (extended by 2 extra sequences for renumbering)
    aln = self.alignment.Copy()
    # Get chains and renumber according to alignment (for lDDT)
    reference = Renumber(
      aln.GetSequence(0),
      old_number_label=self._old_number_label).CreateFullView()
    refseq = seq.CreateSequence(
      "reference_renumbered",
      aln.GetSequence(0).GetString())
    refseq.AttachView(reference)
    aln.AddSequence(refseq)
    model = Renumber(
      aln.GetSequence(1),
      old_number_label=self._old_number_label).CreateFullView()
    modelseq = seq.CreateSequence(
      "model_renumbered",
      aln.GetSequence(1).GetString())
    modelseq.AttachView(model)
    aln.AddSequence(modelseq)
    # Filter to CA-only if desired (done after AttachView to not mess it up)
    if self.calpha_only:
      self.lddt_scorer = lDDTScorer(
        references=[reference.Select('aname=CA')],
        model=model.Select('aname=CA'),
        settings=self.settings)
    else:
      self.lddt_scorer = lDDTScorer(
        references=[reference],
        model=model,
        settings=self.settings)
    # Store alignment for later
    self._extended_alignment = aln

###############################################################################
# HELPERS
###############################################################################

# general

def _AlignAtomSeqs(seq_1, seq_2):
  """
  :type seq_1: :class:`ost.seq.SequenceHandle`
  :type seq_2: :class:`ost.seq.SequenceHandle`
  :return: Alignment of two sequences using a global alignment. Views attached
           to the input sequences will remain attached in the aln.
  :rtype:  :class:`~ost.seq.AlignmentHandle` or None if it failed.
  """
  # NOTE: If the two sequence have a greatly different length
  #       a local alignment could be a better choice...
  aln = None
  alns = seq.alg.GlobalAlign(seq_1, seq_2, seq.alg.BLOSUM62)
  if alns: aln = alns[0]
  if not aln:
    LogWarning('Failed to align %s to %s' % (seq_1.name, seq_2.name))
    LogWarning('%s:  %s' % (seq_1.name, seq_1.string))
    LogWarning('%s:  %s' % (seq_2.name, seq_2.string))
  return aln

def _FixSelectChainNames(ch_names):
  """
  :return: String to be used with Select(cname=<RETURN>). Takes care of joining
           and putting quotation marks where needed.
  :rtype:  :class:`str`
  :param ch_names: Some iterable list of chain names (:class:`str` items).
  """
  return ','.join(mol.QueryQuoteName(ch_name) for ch_name in ch_names)

# QS entity

def _CleanInputEntity(ent):
  """
  :param ent: The OST entity to be cleaned.
  :type ent:  :class:`EntityHandle` or :class:`EntityView`
  :return: A tuple of 3 items: :attr:`QSscoreEntity.ent`,
                               :attr:`QSscoreEntity.removed_chains`,
                               :attr:`QSscoreEntity.calpha_only`
  """
  # find chains to remove
  removed_chains = []
  for ch in ent.chains:
    # we remove chains if they are small-peptides or if the contain no aa
    # or if they contain only unknown or modified residues
    if    ch.name in ['-', '_'] \
       or ch.residue_count < 20 \
       or not any(r.chem_type.IsAminoAcid() for r in ch.residues) \
       or not (set(r.one_letter_code for r in ch.residues) - {'?', 'X'}):
      removed_chains.append(ch.name)

  # remove them from *ent*
  if removed_chains:
    view = ent.Select('cname!=%s' % _FixSelectChainNames(set(removed_chains)))
    ent_new = mol.CreateEntityFromView(view, True)
    ent_new.SetName(ent.GetName())
  else:
    ent_new = ent

  # check if CA only
  calpha_only = False
  if ent_new.atom_count > 0 and ent_new.Select('aname=CB').atom_count == 0:
    LogInfo('Structure %s is a CA only structure!' % ent_new.GetName())
    calpha_only = True

  # report and return
  if removed_chains:
    LogInfo('Chains removed from %s: %s' \
            % (ent_new.GetName(), ''.join(removed_chains)))
  LogInfo('Chains in %s: %s' \
          % (ent_new.GetName(), ''.join([c.name for c in ent_new.chains])))
  return ent_new, removed_chains, calpha_only

def _GetCAOnlyEntity(ent):
  """
  :param ent: Entity to process
  :type ent:  :class:`EntityHandle` or :class:`EntityView`
  :return: New entity with only CA and only one atom per residue
           (see :attr:`QSscoreEntity.ca_entity`)
  """
  # cook up CA only view (diff from Select = guaranteed 1 atom per residue)
  ca_view = ent.CreateEmptyView()
  # add chain by chain
  for res in ent.residues:
    ca_atom = res.FindAtom("CA")
    if ca_atom.IsValid():
      ca_view.AddAtom(ca_atom)
  # finalize
  return mol.CreateEntityFromView(ca_view, False)

def _GetChemGroups(qs_ent, seqid_thr=95.):
  """
  :return: Intra-complex group of chemically identical polypeptide chains
           (see :attr:`QSscoreEntity.chem_groups`)

  :param qs_ent: Entity to process
  :type qs_ent:  :class:`QSscoreEntity`
  :param seqid_thr: Threshold used to decide when two chains are identical.
                    95 percent tolerates the few mutations crystallographers
                    like to do.
  :type seqid_thr:  :class:`float`
  """
  # get data from qs_ent
  ca_chains = qs_ent.ca_chains
  chain_names = sorted(ca_chains.keys())
  # get pairs of identical chains
  # NOTE: this scales quadratically with number of chains and may be optimized
  #       -> one could merge it with "merge transitive pairs" below...
  id_seqs = []
  for ch_1, ch_2 in itertools.combinations(chain_names, 2):
    aln = qs_ent.GetAlignment(ch_1, ch_2)
    if aln and seq.alg.SequenceIdentity(aln) > seqid_thr:
      id_seqs.append((ch_1, ch_2))
  # trivial case: no matching pairs
  if not id_seqs:
    return [[name] for name in chain_names]

  # merge transitive pairs
  groups = []
  for ch_1, ch_2 in id_seqs:
    found = False
    for g in groups:
      if ch_1 in g or ch_2 in g:
        found = True
        g.add(ch_1)
        g.add(ch_2)
    if not found:
      groups.append(set([ch_1, ch_2]))
  # sort internally based on sequence length
  chem_groups = []
  for g in groups:
    ranked_g = sorted([(-ca_chains[ch].length, ch) for ch in g])
    chem_groups.append([ch for _,ch in ranked_g])
  # add other dissimilar chains
  for ch in chain_names:
    if not any(ch in g for g in chem_groups):
      chem_groups.append([ch])
  
  return chem_groups

def _GetAngles(Rt):
  """Computes the Euler angles given a transformation matrix.

  :param Rt: Rt operator.
  :type Rt:  :class:`ost.geom.Mat4`
  :return: A :class:`tuple` of angles for each axis (x,y,z)
  """
  rot = np.asarray(Rt.ExtractRotation().data).reshape(3,3)
  tx = np.arctan2(rot[2,1], rot[2,2])
  if tx < 0:
    tx += 2*np.pi
  ty = np.arctan2(rot[2,0], np.sqrt(rot[2,1]**2 + rot[2,2]**2))
  if ty < 0:
    ty += 2*np.pi
  tz = np.arctan2(rot[1,0], rot[0,0])
  if tz < 0:
    tz += 2*np.pi
  return tx,ty,tz

# QS scorer

def _GetChemGroupsMapping(qs_ent_1, qs_ent_2):
  """
  :return: Inter-complex mapping of chemical groups
           (see :attr:`QSscorer.chem_mapping`)

  :param qs_ent_1: See :attr:`QSscorer.qs_ent_1`
  :param qs_ent_2: See :attr:`QSscorer.qs_ent_2`
  """
  # get chem. groups and unique representative
  chem_groups_1 = qs_ent_1.chem_groups
  chem_groups_2 = qs_ent_2.chem_groups
  repr_chains_1 = {x[0]: tuple(x) for x in chem_groups_1}
  repr_chains_2 = {x[0]: tuple(x) for x in chem_groups_2}

  # if entities do not have different number of unique chains, we get the
  # mapping for the smaller set
  swapped = False
  if len(repr_chains_2) < len(repr_chains_1):
    repr_chains_1, repr_chains_2 = repr_chains_2, repr_chains_1
    qs_ent_1, qs_ent_2 = qs_ent_2, qs_ent_1
    swapped = True

  # find the closest to each chain between the two entities
  # NOTE: this may still be sensible to orthology problem
  # -> currently we use a global alignment and seq. id. to rank pairs
  # -> we also tried local alignments and weighting the seq. id. by the
  #    coverages of the alignments (gapless string in aln. / seq. length)
  #    but global aln performed better...
  chain_pairs = []
  ca_chains_1 = qs_ent_1.ca_chains
  ca_chains_2 = qs_ent_2.ca_chains
  for ch_1 in list(repr_chains_1.keys()):
    for ch_2 in list(repr_chains_2.keys()):
      aln = _AlignAtomSeqs(ca_chains_1[ch_1], ca_chains_2[ch_2])
      if aln:
        chains_seqid = seq.alg.SequenceIdentity(aln)
        LogVerbose('Sequence identity', ch_1, ch_2, 'seqid=%.2f' % chains_seqid)
        chain_pairs.append((chains_seqid, ch_1, ch_2))

  # get top matching groups first
  chain_pairs = sorted(chain_pairs, reverse=True)
  chem_mapping = {}
  for _, c1, c2 in chain_pairs:
    skip = False
    for a,b in chem_mapping.items():
      if repr_chains_1[c1] == a or repr_chains_2[c2] == b:
        skip = True
        break
    if not skip:
      chem_mapping[repr_chains_1[c1]] = repr_chains_2[c2]
  if swapped:
    chem_mapping = {y: x for x, y in chem_mapping.items()}
    qs_ent_1, qs_ent_2 = qs_ent_2, qs_ent_1

  # notify chains without partner
  mapped_1 = set([i for s in list(chem_mapping.keys()) for i in s])
  chains_1 = set([c.name for c in qs_ent_1.ent.chains])
  if chains_1 - mapped_1:
    LogWarning('Unmapped Chains in %s: %s'
               % (qs_ent_1.GetName(), ','.join(list(chains_1 - mapped_1))))

  mapped_2 = set([i for s in list(chem_mapping.values()) for i in s])
  chains_2 = set([c.name for c in qs_ent_2.ent.chains])
  if chains_2 - mapped_2:
    LogWarning('Unmapped Chains in %s: %s'
               % (qs_ent_2.GetName(), ','.join(list(chains_2 - mapped_2))))
  
  # check if we have any chains left
  LogInfo('Chemical chain-groups mapping: ' + str(chem_mapping))
  if len(mapped_1) < 1 or len(mapped_2) < 1:
    raise QSscoreError('Less than 1 chains left in chem_mapping.')
  return chem_mapping

def _SelectFew(l, max_elements):
  """Return l or copy of l with at most *max_elements* entries."""
  n_el = len(l)
  if n_el <= max_elements:
    return l
  else:
    # cheap integer ceiling (-1 to ensure that x*max_elements gets d_el = x)
    d_el = ((n_el - 1) // max_elements) + 1
    new_l = list()
    for i in range(0, n_el, d_el):
      new_l.append(l[i])
    return new_l

def _GetAlignedResidues(qs_ent_1, qs_ent_2, chem_mapping, max_ca_per_chain,
                        clustalw_bin):
  """
  :return: Tuple of two :class:`~ost.mol.EntityView` objects containing subsets
           of *qs_ent_1* and *qs_ent_2*. Two entities are later created from
           those views (see :attr:`QSscorer.ent_to_cm_1` and
           :attr:`QSscorer.ent_to_cm_2`)

  :param qs_ent_1: See :attr:`QSscorer.qs_ent_1`
  :param qs_ent_2: See :attr:`QSscorer.qs_ent_2`
  :param chem_mapping: See :attr:`QSscorer.chem_mapping`
  :param max_ca_per_chain: See :attr:`QSscorer.max_ca_per_chain_for_cm`
  """
  # make sure name doesn't contain spaces and is unique
  def _FixName(seq_name, seq_names):
    # get rid of spaces and make it unique
    seq_name = seq_name.replace(' ', '-')
    while seq_name in seq_names:
      seq_name += '-'
    return seq_name
  # resulting views into CA entities using CA chain sequences
  ent_view_1 = qs_ent_1.ca_entity.CreateEmptyView()
  ent_view_2 = qs_ent_2.ca_entity.CreateEmptyView()
  ca_chains_1 = qs_ent_1.ca_chains
  ca_chains_2 = qs_ent_2.ca_chains
  # go through all mapped chemical groups
  for group_1, group_2 in chem_mapping.items():
    # MSA with ClustalW
    seq_list = seq.CreateSequenceList()
    # keep sequence-name (must be unique) to view mapping
    seq_to_empty_view = dict()
    for ch in group_1:
      sequence = ca_chains_1[ch].Copy()
      sequence.name = _FixName(qs_ent_1.GetName() + '.' + ch, seq_to_empty_view)
      seq_to_empty_view[sequence.name] = ent_view_1
      seq_list.AddSequence(sequence)
    for ch in group_2:
      sequence = ca_chains_2[ch].Copy()
      sequence.name = _FixName(qs_ent_2.GetName() + '.' + ch, seq_to_empty_view)
      seq_to_empty_view[sequence.name] = ent_view_2
      seq_list.AddSequence(sequence)
    alnc = ClustalW(seq_list, clustalw=clustalw_bin)

    # collect aligned residues (list of lists of sequence_count valid residues)
    residue_lists = list()
    sequence_count = alnc.sequence_count
    for col in alnc:
      residues = [col.GetResidue(ir) for ir in range(sequence_count)]
      if all([r.IsValid() for r in residues]):
        residue_lists.append(residues)
    # check number of aligned residues
    if len(residue_lists) < 5:
      raise QSscoreError('Not enough aligned residues.')
    elif max_ca_per_chain > 0:
      residue_lists = _SelectFew(residue_lists, max_ca_per_chain)
    # check what view is for which residue
    res_views = [seq_to_empty_view[alnc.sequences[ir].name] \
                 for ir in range(sequence_count)]
    # add to view (note: only 1 CA atom per residue in here)
    for res_list in residue_lists:
      for res, view in zip(res_list, res_views):
        view.AddResidue(res, mol.ViewAddFlag.INCLUDE_ATOMS)
  # done
  return ent_view_1, ent_view_2


def _FindSymmetry(qs_ent_1, qs_ent_2, ent_to_cm_1, ent_to_cm_2, chem_mapping):
  """
  :return: A pair of comparable symmetry groups (for :attr:`QSscorer.symm_1`
           and :attr:`QSscorer.symm_2`) between the two structures.
           Empty lists if no symmetry identified.
  
  :param qs_ent_1: See :attr:`QSscorer.qs_ent_1`
  :param qs_ent_2: See :attr:`QSscorer.qs_ent_2`
  :param ent_to_cm_1: See :attr:`QSscorer.ent_to_cm_1`
  :param ent_to_cm_2: See :attr:`QSscorer.ent_to_cm_2`
  :param chem_mapping: See :attr:`QSscorer.chem_mapping`
  """
  LogInfo('Identifying Symmetry Groups...')

  # get possible symmetry groups
  symm_subg_1 = _GetSymmetrySubgroups(qs_ent_1, ent_to_cm_1,
                                      list(chem_mapping.keys()))
  symm_subg_2 = _GetSymmetrySubgroups(qs_ent_2, ent_to_cm_2,
                                      list(chem_mapping.values()))

  # check combination of groups
  LogInfo('Selecting Symmetry Groups...')
  # check how many mappings are to be checked for each pair of symmetry groups
  # -> lower number here will speed up worst-case runtime of _GetChainMapping
  # NOTE: this is tailored to speed up brute force all vs all mapping test
  #       for preferred _CheckClosedSymmetry this is suboptimal!
  best_symm = []
  for symm_1, symm_2 in itertools.product(symm_subg_1, symm_subg_2):
    s1 = symm_1[0]
    s2 = symm_2[0]
    if len(s1) != len(s2):
      LogDebug('Discarded', str(symm_1), str(symm_2),
               ': different size of symmetry groups')
      continue
    count = _CountSuperpositionsAndMappings(symm_1, symm_2, chem_mapping)
    nr_mapp = count['intra']['mappings'] + count['inter']['mappings']
    LogDebug('OK', str(symm_1), str(symm_2), ': %s mappings' % nr_mapp)
    best_symm.append((nr_mapp, symm_1, symm_2))

  for _, symm_1, symm_2 in sorted(best_symm):
    s1 = symm_1[0]
    s2 = symm_2[0]
    group_1 = ent_to_cm_1.Select('cname=%s' % _FixSelectChainNames(s1))
    group_2 = ent_to_cm_2.Select('cname=%s' % _FixSelectChainNames(s2))
    # check if by superposing a pair of chains within the symmetry group to
    # superpose all chains within the symmetry group
    # -> if successful, the symmetry groups are compatible
    closed_symm = _CheckClosedSymmetry(group_1, group_2, symm_1, symm_2,
                                       chem_mapping, 6, 0.8, find_best=False)

    if closed_symm:
      return symm_1, symm_2

  # nothing found
  LogInfo('No coherent symmetry identified between structures')
  return [], []


def _GetChainMapping(ent_1, ent_2, symm_1, symm_2, chem_mapping,
                     max_mappings_extensive):
  """
  :return: Tuple with mapping from *ent_1* to *ent_2* (see
           :attr:`QSscorer.chain_mapping`) and scheme used (see
           :attr:`QSscorer.chain_mapping_scheme`)

  :param ent_1: See :attr:`QSscorer.ent_to_cm_1`
  :param ent_2: See :attr:`QSscorer.ent_to_cm_2`
  :param symm_1: See :attr:`QSscorer.symm_1`
  :param symm_2: See :attr:`QSscorer.symm_2`
  :param chem_mapping: See :attr:`QSscorer.chem_mapping`
  :param max_mappings_extensive: See :attr:`QSscorer.max_mappings_extensive`
  """
  LogInfo('Symmetry-groups used in %s: %s' % (ent_1.GetName(), str(symm_1)))
  LogInfo('Symmetry-groups used in %s: %s' % (ent_2.GetName(), str(symm_2)))

  # quick check for closed symmetries
  thresholds = [(    'strict', 4, 0.8),
                (  'tolerant', 6, 0.4),
                ('permissive', 8, 0.2)]
  for scheme, sup_thr, sup_fract in thresholds:
    # check if by superposing a pair of chains within the symmetry group to
    # superpose all chains of the two oligomers
    # -> this also returns the resulting mapping of chains
    mapping = _CheckClosedSymmetry(ent_1, ent_2, symm_1, symm_2,
                                   chem_mapping, sup_thr, sup_fract)
    if mapping:
      LogInfo('Closed Symmetry with %s parameters' % scheme)
      if scheme == 'permissive':
        LogWarning('Permissive thresholds used for overlap based mapping ' + \
                   'detection: check mapping manually: %s' % mapping)
      return mapping, scheme
  
  # NOTE that what follows below is sub-optimal:
  # - if the two structures don't fit at all, we may map chains rather randomly
  # - we may need a lot of operations to find globally "best" mapping
  # -> COST = O(N^2) * O(N!)
  #    (first = all pairwise chain-superpose, latter = all chain mappings)
  # -> a greedy chain mapping choice algo (pick best RMSD for each one-by-one)
  #    could reduce this to O(N^2) * O(N^2) but it would need to be validated
  # - one could also try some other heuristic based on center of atoms of chains
  #   -> at this point we get a bad mapping anyways...

  # if we get here, we will need to check many combinations of mappings
  # -> check how many mappings must be checked and limit
  count = _CountSuperpositionsAndMappings(symm_1, symm_2, chem_mapping)
  LogInfo('Intra Symmetry-group mappings to check: %s' \
          % count['intra']['mappings'])
  LogInfo('Inter Symmetry-group mappings to check: %s' \
          % count['inter']['mappings'])
  nr_mapp = count['intra']['mappings'] + count['inter']['mappings']
  if nr_mapp > max_mappings_extensive:
    raise QSscoreError('Too many possible mappings: %s' % nr_mapp)

  # to speed up the computations we cache chain views and RMSDs
  cached_rmsd = _CachedRMSD(ent_1, ent_2)

  # get INTRA symmetry group chain mapping
  check = 0
  intra_mappings = [] # list of (RMSD, c1, c2, mapping) tuples (best superpose)
  # limit chem mapping to reference symmetry group
  ref_symm_1 = symm_1[0]
  ref_symm_2 = symm_2[0]
  asu_chem_mapping = _LimitChemMapping(chem_mapping, ref_symm_1, ref_symm_2)
  # for each chemically identical group
  for g1, g2 in asu_chem_mapping.items():
    # try to superpose all possible pairs
    for c1, c2 in itertools.product(g1, g2):
      # get superposition transformation
      LogVerbose('Superposing chains: %s, %s' % (c1, c2))
      res = cached_rmsd.GetSuperposition(c1, c2)
      # compute RMSD of possible mappings
      cur_mappings = [] # list of (RMSD, mapping) tuples
      for mapping in _ListPossibleMappings(c1, c2, asu_chem_mapping):
        check += 1
        chains_rmsd = cached_rmsd.GetMappedRMSD(mapping, res.transformation)
        cur_mappings.append((chains_rmsd, mapping))
        LogVerbose(str(mapping), str(chains_rmsd))
      best_rmsd, best_mapp = min(cur_mappings)
      intra_mappings.append((best_rmsd, c1, c2, best_mapp))
  # best chain-chain superposition defines the intra asu mapping
  _, intra_asu_c1, intra_asu_c2, intra_asu_mapp = min(intra_mappings)

  # if only one asu is present in any of the complexes, we're done
  if len(symm_1) == 1 or len(symm_2) == 1:
    mapping = intra_asu_mapp
  else:
    # the mapping is the element position within the asu chain list
    # -> this speed up a lot, assuming that the order of chain in asu groups
    #    is following the order of symmetry operations
    index_asu_c1 = ref_symm_1.index(intra_asu_c1)
    index_asu_c2 = ref_symm_2.index(intra_asu_c2)
    index_mapp = {ref_symm_1.index(s1): ref_symm_2.index(s2) \
                  for s1, s2 in intra_asu_mapp.items()}
    LogInfo('Intra symmetry-group mapping: %s' % str(intra_asu_mapp))

    # get INTER symmetry group chain mapping
    # we take the first symmetry group from the complex having the most
    if len(symm_1) < len(symm_2):
      symm_combinations = list(itertools.product(symm_1, [symm_2[0]]))
    else:
      symm_combinations = list(itertools.product([symm_1[0]], symm_2))

    full_asu_mapp = {tuple(symm_1): tuple(symm_2)}
    full_mappings = [] # list of (RMSD, mapping) tuples
    for s1, s2 in symm_combinations:
      # get superposition transformation (we take best chain-pair in asu)
      LogVerbose('Superposing symmetry-group: %s, %s in %s, %s' \
                 % (s1[index_asu_c1], s2[index_asu_c2], s1, s2))
      res = cached_rmsd.GetSuperposition(s1[index_asu_c1], s2[index_asu_c2])
      # compute RMSD of possible mappings
      for asu_mapp in _ListPossibleMappings(s1, s2, full_asu_mapp):
        check += 1
        # need to extract full chain mapping (use indexing)
        mapping = {}
        for a, b in asu_mapp.items():
          for id_a, id_b in index_mapp.items():
            mapping[a[id_a]] = b[id_b]
        chains_rmsd = cached_rmsd.GetMappedRMSD(mapping, res.transformation)
        full_mappings.append((chains_rmsd, mapping))
        LogVerbose(str(mapping), str(chains_rmsd))

    if check != nr_mapp:
      LogWarning('Mapping number estimation was wrong') # sanity check

    # return best (lowest RMSD) mapping
    mapping = min(full_mappings)[1]

  LogWarning('Extensive search used for mapping detection (fallback). This ' + \
             'approach has known limitations. Check mapping manually: %s' \
             % mapping)
  return mapping, 'extensive'


def _GetSymmetrySubgroups(qs_ent, ent, chem_groups):
  """
  Identify the symmetry (either cyclic C or dihedral D) of the protein and find
  all possible symmetry subgroups. This is done testing all combination of chain
  superposition and clustering them by the angles (D) or the axis (C) of the Rt
  operator.

  Groups of superposition which can fully reconstruct the structure are possible
  symmetry subgroups.

  :param qs_ent: Entity with cached angles and axis.
  :type qs_ent:  :class:`QSscoreEntity`
  :param ent: Entity from which to extract chains (perfect alignment assumed
              for superposition as in :attr:`QSscorer.ent_to_cm_1`).
  :type ent:  :class:`EntityHandle` or :class:`EntityView`
  :param chem_groups: List of tuples/lists of chain names in *ent*. Each list
                      contains all chains belonging to a chem. group (could be
                      from :attr:`QSscoreEntity.chem_groups` or from "keys()"
                      or "values()" of :attr:`QSscorer.chem_mapping`)

  :return: A list of possible symmetry subgroups (each in same format as
           :attr:`QSscorer.symm_1`). If no symmetry is found, we return a list
           with a single symmetry subgroup with a single group with all chains.
  """
  # get angles / axis for pairwise transformations within same chem. group
  angles = {}
  axis = {}
  for cnames in chem_groups:
    for c1, c2 in itertools.combinations(cnames, 2):
      cur_angles = qs_ent.GetAngles(c1, c2)
      if not np.any(np.isnan(cur_angles)):
        angles[(c1,c2)] = cur_angles
      cur_axis = qs_ent.GetAxis(c1, c2)
      if not np.any(np.isnan(cur_axis)):
        axis[(c1,c2)] = cur_axis

  # cluster superpositions angles at different thresholds
  symm_groups = []
  LogVerbose('Possible symmetry-groups for: %s' % ent.GetName())
  for angle_thr in np.arange(0.1, 1, 0.1):
    d_symm_groups = _GetDihedralSubgroups(ent, chem_groups, angles, angle_thr)
    if d_symm_groups:
      for group in d_symm_groups:
        if group not in symm_groups:
          symm_groups.append(group)
          LogVerbose('Dihedral: %s' % group)
      LogInfo('Symmetry threshold %.1f used for angles of %s' \
              % (angle_thr, ent.GetName()))
      break
  
  # cluster superpositions axis at different thresholds
  for axis_thr in np.arange(0.1, 1, 0.1):
    c_symm_groups = _GetCyclicSubgroups(ent, chem_groups, axis, axis_thr)
    if c_symm_groups:
      for group in c_symm_groups:
        if group not in symm_groups:
          symm_groups.append(group)
          LogVerbose('Cyclic: %s' % group)
      LogInfo('Symmetry threshold %.1f used for axis of %s' \
              % (axis_thr, ent.GetName()))
      break

  # fall back to single "group" containing all chains if none found
  if not symm_groups:
    LogInfo('No symmetry-group detected in %s' % ent.GetName())
    symm_groups = [[tuple([c for g in chem_groups for c in g])]]

  return symm_groups

def _GetDihedralSubgroups(ent, chem_groups, angles, angle_thr):
  """
  :return: Dihedral subgroups for :func:`_GetSymmetrySubgroups`
           (same return type as there). Empty list if fail.

  :param ent: See :func:`_GetSymmetrySubgroups`
  :param chem_groups: See :func:`_GetSymmetrySubgroups`
  :param angles: :class:`dict` (key = chain-pair-tuple, value = Euler angles)
  :param angle_thr: Euler angles distance threshold for clustering (float).
  """
  # cluster superpositions angles
  if len(angles) == 0:
    # nothing to be done here
    return []
  else:
    same_angles = _ClusterData(angles, angle_thr, _AngleArrayDistance)

  # get those which are non redundant and covering all chains
  symm_groups = []
  for clst in list(same_angles.values()):
    group = list(clst.keys())
    if _ValidChainGroup(group, ent):
      if len(chem_groups) > 1:
        # if hetero, we want to group toghether different chains only
        symm_groups.append(list(zip(*group)))
      else:
        # if homo, we also want pairs
        symm_groups.append(group)
        symm_groups.append(list(zip(*group)))
  return symm_groups

def _GetCyclicSubgroups(ent, chem_groups, axis, axis_thr):
  """
  :return: Cyclic subgroups for :func:`_GetSymmetrySubgroups`
           (same return type as there). Empty list if fail.

  :param ent: See :func:`_GetSymmetrySubgroups`
  :param chem_groups: See :func:`_GetSymmetrySubgroups`
  :param angles: :class:`dict` (key = chain-pair-tuple, value = rotation axis)
  :param angle_thr: Axis distance threshold for clustering (float).
  """
  # cluster superpositions axis
  if len(axis) == 0:
    # nothing to be done here
    return []
  else:
    same_axis = _ClusterData(axis, axis_thr, _AxisDistance)

  # use to get grouping
  symm_groups = []
  for clst in list(same_axis.values()):
    all_chain = [item for sublist in list(clst.keys()) for item in sublist]
    if len(set(all_chain)) == ent.chain_count:
      # if we have an hetero we can exploit cyclic symmetry for grouping
      if len(chem_groups) > 1:
        ref_chem_group = chem_groups[0]
        # try two criteria for grouping
        cyclic_group_closest = []
        cyclic_group_iface = []
        for c1 in ref_chem_group:
          group_closest = [c1]
          group_iface = [c1]
          for chains in chem_groups[1:]:
            # center of atoms distance
            closest = _GetClosestChain(ent, c1, chains)
            # chain with bigger interface
            closest_iface = _GetClosestChainInterface(ent, c1, chains)
            group_closest.append(closest)
            group_iface.append(closest_iface)
          cyclic_group_closest.append(tuple(group_closest))
          cyclic_group_iface.append(tuple(group_iface))
        if _ValidChainGroup(cyclic_group_closest, ent):
          symm_groups.append(cyclic_group_closest)
        if _ValidChainGroup(cyclic_group_iface, ent):
          symm_groups.append(cyclic_group_iface)
      # if it is a homo, there's not much we can group
      else:
        symm_groups.append(chem_groups)
  return symm_groups

def _ClusterData(data, thr, metric):
  """Wrapper for fclusterdata to get dict of clusters.
  
  :param data: :class:`dict` (keys for ID, values used for clustering)
  :return: :class:`dict` {cluster_idx: {data-key: data-value}}
  """
  # special case length 1
  if len(data) == 1:
    return {0: {list(data.keys())[0]: list(data.values())[0]} }
  # do the clustering
  cluster_indices = fclusterdata(np.asarray(list(data.values())), thr,
                                 method='complete', criterion='distance',
                                 metric=metric)
  # fclusterdata output is cluster ids -> put into dict of clusters
  res = {}
  for cluster_idx, data_key in zip(cluster_indices, list(data.keys())):
    if not cluster_idx in res:
      res[cluster_idx] = {}
    res[cluster_idx][data_key] = data[data_key]
  return res

def _AngleArrayDistance(u, v):
  """
  :return: Average angular distance of two arrays of angles.
  :param u: Euler angles.
  :param v: Euler angles.
  """
  dist = []
  for a,b in zip(u,v):
    delta = abs(a-b)
    if delta > np.pi:
      delta = abs(2*np.pi - delta)
    dist.append(delta)
  return np.mean(dist)

def _AxisDistance(u, v):
  """
  :return: Euclidean distance between two rotation axes. Axes can point in
           either direction so we ensure to use the closer one.
  :param u: Rotation axis.
  :param v: Rotation axis.
  """
  # get axes as arrays (for convenience)
  u = np.asarray(u)
  v = np.asarray(v)
  # get shorter of two possible distances (using v or -v)
  dv1 = u - v
  dv2 = u + v
  d1 = np.dot(dv1, dv1)
  d2 = np.dot(dv2, dv2)
  return np.sqrt(min(d1, d2))

def _GetClosestChain(ent, ref_chain, chains):
  """
  :return: Chain closest to *ref_chain* based on center of atoms distance.
  :rtype:  :class:`str`
  :param ent: See :func:`_GetSymmetrySubgroups`
  :param ref_chain: We look for chains closest to this one
  :type ref_chain:  :class:`str`
  :param chains: We only consider these chains
  :type chains:  :class:`list` of :class:`str`
  """
  contacts = []
  ref_center = ent.FindChain(ref_chain).center_of_atoms
  for ch in chains:
    center = ent.FindChain(ch).center_of_atoms
    contacts.append((geom.Distance(ref_center, center), ch))
  closest_chain = min(contacts)[1]
  return closest_chain

def _GetClosestChainInterface(ent, ref_chain, chains):
  """
  :return: Chain with biggest interface (within 10 A) to *ref_chain*.
  :rtype:  :class:`str`
  :param ent: See :func:`_GetSymmetrySubgroups`
  :param ref_chain: We look for chains closest to this one
  :type ref_chain:  :class:`str`
  :param chains: We only consider these chains
  :type chains:  :class:`list` of :class:`str`
  """
  # NOTE: this is computed on a subset of the full biounit and might be
  # inaccurate. Also it could be extracted from QSscoreEntity.contacts.
  closest = []
  for ch in chains:
    iface_view = ent.Select('cname="%s" and 10 <> [cname="%s"]' \
                            % (ref_chain, ch))
    nr_res = iface_view.residue_count
    closest.append((nr_res, ch))
  closest_chain = max(closest)[1]
  return closest_chain

def _ValidChainGroup(group, ent):
  """
  :return: True, if *group* has unique chain names and as many chains as *ent*
  :rtype:  :class:`bool`
  :param group: Symmetry groups to check
  :type group:  Same as :attr:`QSscorer.symm_1`
  :param ent: See :func:`_GetSymmetrySubgroups`
  """
  all_chain = [item for sublist in group for item in sublist]
  if len(all_chain) == len(set(all_chain)) and\
     len(all_chain) == ent.chain_count:
    return True
  return False

def _LimitChemMapping(chem_mapping, limit_1, limit_2):
  """
  :return: Chem. mapping containing only chains in *limit_1* and *limit_2*
  :rtype:  Same as :attr:`QSscorer.chem_mapping`
  :param chem_mapping: See :attr:`QSscorer.chem_mapping`
  :param limit_1: Limits chain names in chem_mapping.keys()
  :type limit_1:  List/tuple of strings
  :param limit_2: Limits chain names in chem_mapping.values()
  :type limit_2:  List/tuple of strings
  """
  # use set intersection to keep only mapping for chains in limit_X
  limit_1_set = set(limit_1)
  limit_2_set = set(limit_2)
  asu_chem_mapping = dict()
  for key, value in chem_mapping.items():
    new_key = tuple(set(key) & limit_1_set)
    if new_key:
      asu_chem_mapping[new_key] = tuple(set(value) & limit_2_set)
  return asu_chem_mapping


def _CountSuperpositionsAndMappings(symm_1, symm_2, chem_mapping):
  """
  :return: Dictionary of number of mappings and superpositions to be performed.
           Returned as *result[X][Y] = number* with X = "intra" or "inter" and
           Y = "mappings" or "superpositions". The idea is that for each
           pairwise superposition we check all possible mappings.
           We can check the combinations within (intra) a symmetry group and
           once established, we check the combinations between different (inter)
           symmetry groups.
  :param symm_1: See :attr:`QSscorer.symm_1`
  :param symm_2: See :attr:`QSscorer.symm_2`
  :param chem_mapping: See :attr:`QSscorer.chem_mapping`
  """
  # setup
  c = {}
  c['intra'] = {}
  c['inter'] = {}
  c['intra']['mappings'] = 0
  c['intra']['superpositions'] = 0
  c['inter']['mappings'] = 0
  c['inter']['superpositions'] = 0
  # intra ASU mappings within reference symmetry groups
  ref_symm_1 = symm_1[0]
  ref_symm_2 = symm_2[0]
  asu_chem_mapping = _LimitChemMapping(chem_mapping, ref_symm_1, ref_symm_2)
  for g1, g2 in asu_chem_mapping.items():
    chain_superpositions = len(g1) * len(g2)
    c['intra']['superpositions'] += chain_superpositions
    map_per_sup = _PermutationOrCombinations(g1[0], g2[0], asu_chem_mapping)
    c['intra']['mappings'] += chain_superpositions * map_per_sup
  if len(symm_1) == 1 or len(symm_2) == 1:
    return c
  # inter ASU mappings
  asu_superposition = min(len(symm_1), len(symm_2))
  c['inter']['superpositions'] = asu_superposition
  asu = {tuple(symm_1): tuple(symm_2)}
  map_per_sup = _PermutationOrCombinations(ref_symm_1, ref_symm_2, asu)
  c['inter']['mappings'] = asu_superposition * map_per_sup
  return c

def _PermutationOrCombinations(sup1, sup2, chem_mapping):
  """Should match len(_ListPossibleMappings(sup1, sup2, chem_mapping))."""
  # remove superposed elements, put smallest complex as key
  chem_map = {}
  for a,b in chem_mapping.items():
    new_a = tuple([x for x in a if x != sup1])
    new_b = tuple([x for x in b if x != sup2])
    chem_map[new_a] = new_b
  mapp_nr = 1.0
  for a, b in chem_map.items():
    if len(a) == len(b):
      mapp_nr *= factorial(len(a))
    elif len(a) < len(b):
      mapp_nr *= binom(len(b), len(a))
    elif len(a) > len(b):
      mapp_nr *= binom(len(a), len(b))
  return int(mapp_nr)

def _ListPossibleMappings(sup1, sup2, chem_mapping):
  """
  Return a flat list of all possible mappings given *chem_mapping* and keeping
  mapping of *sup1* and *sup2* fixed. For instance if elements are chain names
  this is all the mappings to check for a given superposition.

  Elements in first complex are defined by *chem_mapping.keys()* (list of list
  of elements) and elements in second complex by *chem_mapping.values()*. If
  complexes don't have same number of elements, we map only elements for the one
  with less. Also mapping is only between elements of mapped groups according to
  *chem_mapping*.
           
  :rtype:  :class:`list` of :class:`dict` (key = element in chem_mapping-key,
           value = element in chem_mapping-value)
  :param sup1: Element for which mapping is fixed.
  :type sup1:  Like element in chem_mapping-key
  :param sup2: Element for which mapping is fixed.
  :type sup2:  Like element in chem_mapping-value
  :param chem_mapping: Defines mapping between groups of elements (e.g. result
                       of :func:`_LimitChemMapping`).
  :type chem_mapping:  :class:`dict` with key / value = :class:`tuple`

  :raises: :class:`QSscoreError` if reference complex (first one or one with
           less elements) has more elements for any given mapped group.
  """
  # find smallest complex
  chain1 = [i for s in list(chem_mapping.keys()) for i in s]
  chain2 = [i for s in list(chem_mapping.values()) for i in s]
  swap = False
  if len(chain1) > len(chain2):
    swap = True
  # remove superposed elements, put smallest complex as key
  chem_map = {}
  for a, b in chem_mapping.items():
    new_a = tuple([x for x in a if x != sup1])
    new_b = tuple([x for x in b if x != sup2])
    if swap:
      chem_map[new_b] = new_a
    else:
      chem_map[new_a] = new_b
  # generate permutations or combinations of chemically
  # equivalent chains
  chem_perm = []
  chem_ref = []
  for a, b in chem_map.items():
    if len(a) == len(b):
      chem_perm.append(list(itertools.permutations(b)))
      chem_ref.append(a)
    elif len(a) < len(b):
      chem_perm.append(list(itertools.combinations(b, len(a))))
      chem_ref.append(a)
    else:
      raise QSscoreError('Impossible to define reference group: %s' \
                         % str(chem_map))
  # save the list of possible mappings
  mappings = []
  flat_ref = [i for s in chem_ref for i in s]
  for perm in itertools.product(*chem_perm):
    flat_perm = [i for s in perm for i in s]
    d = {c1: c2 for c1, c2 in zip(flat_ref, flat_perm)}
    if swap:
      d = {v: k for k, v in list(d.items())}
    d.update({sup1: sup2})
    mappings.append(d)
  return mappings


def _CheckClosedSymmetry(ent_1, ent_2, symm_1, symm_2, chem_mapping,
                         sup_thr=4, sup_fract=0.8, find_best=True):
  """
  Quick check if we can superpose two chains and get a mapping for all other
  chains using the same transformation. The mapping is defined by sufficient
  overlap of the transformed chain of *ent_1* onto another chain in *ent_2*.

  :param ent_1: Entity to map to *ent_2* (perfect alignment assumed between
                chains of same chem. group as in :attr:`QSscorer.ent_to_cm_1`).
                Views are ok but only to select full chains!
  :param ent_2: Entity to map to (perfect alignment assumed between
                chains of same chem. group as in :attr:`QSscorer.ent_to_cm_2`).
                Views are ok but only to select full chains!
  :param symm_1: Symmetry groups to use. We only superpose chains within
                 reference symmetry group of *symm_1* and *symm_2*.
                 See :attr:`QSscorer.symm_1`
  :param symm_2: See :attr:`QSscorer.symm_2`
  :param chem_mapping: See :attr:`QSscorer.chem_mapping`.
                       All chains in *ent_1* / *ent_2* must be listed here!
  :param sup_thr: Distance around transformed chain in *ent_1* to check for
                  overlap.
  :type sup_thr:  :class:`float`
  :param sup_fract: Fraction of atoms in chain of *ent_2* that must be
                    overlapped for overlap to be sufficient.
  :type sup_fract:  :class:`float`
  :param find_best: If True, we look for best mapping according to
                    :func:`_GetMappedRMSD`. Otherwise, we return first suitable
                    mapping.
  :type find_best:  :class:`bool`

  :return: Mapping from *ent_1* to *ent_2* or None if none found. Mapping, if
           found, is as in :attr:`QSscorer.chain_mapping`.
  :rtype:  :class:`dict` (:class:`str` / :class:`str`)
  """
  # limit chem mapping to reference symmetry group
  ref_symm_1 = symm_1[0]
  ref_symm_2 = symm_2[0]
  asu_chem_mapping = _LimitChemMapping(chem_mapping, ref_symm_1, ref_symm_2)
  # for each chemically identical group
  rmsd_mappings = []
  for g1, g2 in asu_chem_mapping.items():
    # try to superpose all possible pairs
    # -> note that some chain-chain combinations may work better than others
    #    to superpose the full oligomer (e.g. if some chains are open/closed)
    for c1, c2 in itertools.product(g1, g2):
      # get superposition transformation
      chain_1 = ent_1.Select('cname="%s"' % c1)
      chain_2 = ent_2.Select('cname="%s"' % c2)
      res = mol.alg.SuperposeSVD(chain_1, chain_2, apply_transform=False)
      # look for overlaps
      mapping = _GetSuperpositionMapping(ent_1, ent_2, chem_mapping,
                                         res.transformation, sup_thr,
                                         sup_fract)
      if not mapping:
        continue
      # early abort if we only want the first one
      if not find_best:
        return mapping
      else:
        # get RMSD for sorting
        rmsd = _GetMappedRMSD(ent_1, ent_2, mapping, res.transformation)
        rmsd_mappings.append((rmsd, mapping))
  # return best mapping
  if rmsd_mappings:
    return min(rmsd_mappings, key=lambda x: x[0])[1]
  else:
    return None

def _GetSuperpositionMapping(ent_1, ent_2, chem_mapping, transformation,
                             sup_thr, sup_fract):
  """
  :return: Dict with chain mapping from *ent_1* to *ent_2* or None if failed
           (see :func:`_CheckClosedSymmetry`).
  :param ent_1: See :func:`_CheckClosedSymmetry`
  :param ent_2: See :func:`_CheckClosedSymmetry`
  :param chem_mapping: See :func:`_CheckClosedSymmetry`
  :param transformation: Superposition transformation to be applied to *ent_1*.
  :param sup_thr: See :func:`_CheckClosedSymmetry`
  :param sup_fract: See :func:`_CheckClosedSymmetry`
  """
  # NOTE: this seems to be the comput. most expensive part in QS scoring
  # -> it could be moved to C++ for speed-up
  # -> the next expensive bits are ClustalW and GetContacts btw...

  # swap if needed (ent_1 must be shorter or same)
  if ent_1.chain_count > ent_2.chain_count:
    swap = True
    ent_1, ent_2 = ent_2, ent_1
    transformation = geom.Invert(transformation)
    chem_pairs = list(zip(list(chem_mapping.values()), list(chem_mapping.keys())))
  else:
    swap = False
    chem_pairs = iter(chem_mapping.items())
  # sanity check
  if ent_1.chain_count == 0:
    return None
  # extract valid partners for each chain in ent_1
  chem_partners = dict()
  for cg1, cg2 in chem_pairs:
    for ch in cg1:
      chem_partners[ch] = set(cg2)
  # get mapping for each chain in ent_1
  mapping = dict()
  mapped_chains = set()
  for ch_1 in ent_1.chains:
    # get all neighbors split by chain (NOTE: this muight be moved to C++)
    ch_atoms = {ch_2.name: set() for ch_2 in ent_2.chains}
    for a_1 in ch_1.handle.atoms:
      transformed_pos = geom.Vec3(transformation * geom.Vec4(a_1.pos))
      a_within = ent_2.FindWithin(transformed_pos, sup_thr)
      for a_2 in a_within:
        ch_atoms[a_2.chain.name].add(a_2.hash_code)
    # get one with most atoms in overlap
    max_num, max_name = max((len(atoms), name)
                            for name, atoms in ch_atoms.items())
    # early abort if overlap fraction not good enough
    ch_2 = ent_2.FindChain(max_name)
    if max_num == 0 or max_num / float(ch_2.handle.atom_count) < sup_fract:
      return None
    # early abort if mapped to chain of different chem. group
    ch_1_name = ch_1.name
    if ch_1_name not in chem_partners:
      raise RuntimeError("Chem. mapping doesn't contain all chains!")
    if max_name not in chem_partners[ch_1_name]:
      return None
    # early abort if multiple ones mapped to same chain
    if max_name in mapped_chains:
      return None
    # set mapping
    mapping[ch_1_name] = max_name
    mapped_chains.add(max_name)
  # unswap if needed and return
  if swap:
    mapping = {v: k for k, v in mapping.items()}
  return mapping

def _GetMappedRMSD(ent_1, ent_2, chain_mapping, transformation):
  """
  :return: RMSD between complexes considering chain mapping.
  :param ent_1: Entity mapped to *ent_2* (perfect alignment assumed between
                mapped chains as in :attr:`QSscorer.ent_to_cm_1`).
  :param ent_2: Entity which was mapped to (perfect alignment assumed between
                mapped chains as in :attr:`QSscorer.ent_to_cm_2`).
  :param chain_mapping: See :attr:`QSscorer.chain_mapping`
  :param transformation: Superposition transformation to be applied to *ent_1*.
  """
  # collect RMSDs and atom counts chain by chain and combine afterwards
  rmsds = []
  atoms = []
  for c1, c2 in chain_mapping.items():
    # get views and atom counts
    chain_1 = ent_1.Select('cname="%s"' % c1)
    chain_2 = ent_2.Select('cname="%s"' % c2)
    atom_count = chain_1.atom_count
    if atom_count != chain_2.atom_count:
      raise RuntimeError('Chains in _GetMappedRMSD must be perfectly aligned!')
    # get RMSD
    rmsd = mol.alg.CalculateRMSD(chain_1, chain_2, transformation)
    # update lists
    rmsds.append(rmsd)
    atoms.append(float(atom_count))
  # combine (reminder: RMSD = sqrt(sum(atom_distance^2)/atom_count))
  rmsd = np.sqrt( sum([a * r**2 for a, r in zip(atoms, rmsds)]) / sum(atoms) )
  return rmsd

class _CachedRMSD:
  """Helper object for repetitive RMSD calculations.
  Meant to speed up :func:`_GetChainMapping` but could also be used to replace
  :func:`_GetMappedRMSD` in :func:`_CheckClosedSymmetry`.

  :param ent_1: See :attr:`QSscorer.ent_to_cm_1`
  :param ent_2: See :attr:`QSscorer.ent_to_cm_2`
  """
  def __init__(self, ent_1, ent_2):
    # initialize caches and keep entities
    self.ent_1 = ent_1
    self.ent_2 = ent_2
    self._chain_views_1 = dict()
    self._chain_views_2 = dict()
    self._pair_rmsd = dict()  # key = (c1,c2), value = (atom_count,rmsd)

  def GetChainView1(self, cname):
    """Get cached view on chain *cname* for :attr:`ent_1`."""
    if cname not in self._chain_views_1:
      self._chain_views_1[cname] = self.ent_1.Select('cname="%s"' % cname)
    return self._chain_views_1[cname]

  def GetChainView2(self, cname):
    """Get cached view on chain *cname* for :attr:`ent_2`."""
    if cname not in self._chain_views_2:
      self._chain_views_2[cname] = self.ent_2.Select('cname="%s"' % cname)
    return self._chain_views_2[cname]

  def GetSuperposition(self, c1, c2):
    """Get superposition result (no change in entities!) for *c1* to *c2*.
    This invalidates cached RMSD results used in :func:`GetMappedRMSD`.

    :param c1: chain name for :attr:`ent_1`.
    :param c2: chain name for :attr:`ent_2`.
    """
    # invalidate _pair_rmsd
    self._pair_rmsd = dict()
    # superpose
    chain_1 = self.GetChainView1(c1)
    chain_2 = self.GetChainView2(c2)
    if chain_1.atom_count != chain_2.atom_count:
      raise RuntimeError('Chains in GetSuperposition not perfectly aligned!')
    return mol.alg.SuperposeSVD(chain_1, chain_2, apply_transform=False)

  def GetMappedRMSD(self, chain_mapping, transformation):
    """
    :return: RMSD between complexes considering chain mapping.
    :param chain_mapping: See :attr:`QSscorer.chain_mapping`.
    :param transformation: Superposition transformation (e.g. res.transformation
                           for res = :func:`GetSuperposition`).
    """
    # collect RMSDs and atom counts chain by chain and combine afterwards
    rmsds = []
    atoms = []
    for c1, c2 in chain_mapping.items():
      # cached?
      if (c1, c2) in self._pair_rmsd:
        atom_count, rmsd = self._pair_rmsd[(c1, c2)]
      else:
        # get views and atom counts
        chain_1 = self.GetChainView1(c1)
        chain_2 = self.GetChainView2(c2)
        atom_count = chain_1.atom_count
        if atom_count != chain_2.atom_count:
          raise RuntimeError('Chains in GetMappedRMSD not perfectly aligned!')
        # get RMSD
        rmsd = mol.alg.CalculateRMSD(chain_1, chain_2, transformation)
        self._pair_rmsd[(c1, c2)] = (atom_count, rmsd)
      # update lists
      rmsds.append(rmsd)
      atoms.append(float(atom_count))
    # combine (reminder: RMSD = sqrt(sum(atom_distance^2)/atom_count))
    rmsd = np.sqrt( sum([a * r**2 for a, r in zip(atoms, rmsds)]) / sum(atoms) )
    return rmsd


def _CleanUserSymmetry(symm, ent):
  """Clean-up user provided symmetry.

  :param symm: Symmetry group as in :attr:`QSscorer.symm_1`
  :param ent: Entity from which to extract chain names
  :type ent:  :class:`~ost.mol.EntityHandle` or :class:`~ost.mol.EntityView`
  :return: New symmetry group limited to chains in sub-structure *ent*
           (see :attr:`QSscorer.symm_1`). Empty list if invalid symmetry.
  """
  # restrict symm to only contain chains in ent
  chain_names = set([ch.name for ch in ent.chains])
  new_symm = list()
  for symm_group in symm:
    new_group = tuple(ch for ch in symm_group if ch in chain_names)
    if new_group:
      new_symm.append(new_group)
  # report it
  if new_symm != symm:
    LogInfo("Cleaned user symmetry for %s: %s" % (ent.GetName(), new_symm))
  # do all groups have same length?
  lengths = [len(symm_group) for symm_group in new_symm]
  if lengths[1:] != lengths[:-1]:
    LogWarning('User symmetry group of different sizes for %s. Ignoring it!' \
               % (ent.GetName()))
    return []
  # do we cover all chains?
  s_chain_names = set([ch for symm_group in new_symm for ch in symm_group])
  if len(s_chain_names) != sum(lengths):
    LogWarning('User symmetry group for %s has duplicate chains. Ignoring it!' \
               % (ent.GetName()))
    return []
  if s_chain_names != chain_names:
    LogWarning('User symmetry group for %s is missing chains. Ignoring it!' \
               % (ent.GetName()))
    return []
  # ok all good
  return new_symm

def _AreValidSymmetries(symm_1, symm_2):
  """Check symmetry pair for major problems.

  :return: False if any of the two is empty or if they're compatible in size.
  :rtype:  :class:`bool`
  """
  if not symm_1 or not symm_2:
    return False
  if len(symm_1) != 1 or len(symm_2) != 1:
    if not len(symm_1[0]) == len(symm_2[0]):
      LogWarning('Symmetry groups of different sizes. Ignoring them!')
      return False
  return True

def _GetMappedAlignments(ent_1, ent_2, chain_mapping, res_num_alignment):
  """
  :return: Alignments of 2 structures given chain mapping
           (see :attr:`QSscorer.alignments`).
  :param ent_1: Entity containing all chains in *chain_mapping.keys()*.
                Views to this entity attached to first sequence of each aln.
  :param ent_2: Entity containing all chains in *chain_mapping.values()*.
                Views to this entity attached to second sequence of each aln.
  :param chain_mapping: See :attr:`QSscorer.chain_mapping`
  :param res_num_alignment: See :attr:`QSscorer.res_num_alignment`
  """
  alns = list()
  for ch_1_name in sorted(chain_mapping):
    # get both sequences incl. attached view
    ch_1 = ent_1.FindChain(ch_1_name)
    ch_2 = ent_2.FindChain(chain_mapping[ch_1_name])
    if res_num_alignment:
      max_res_num = max([r.number.GetNum() for r in ch_1.residues] +
                        [r.number.GetNum() for r in ch_2.residues])
      ch1_aln = ["-"] * max_res_num
      ch2_aln = ["-"] * max_res_num
      for res in ch_1.residues:
        ch1_aln[res.number.GetNum() - 1] = res.GetOneLetterCode()
      ch1_aln = "".join(ch1_aln)
      seq_1 = seq.CreateSequence(ch_1.name, str(ch1_aln))
      seq_1.AttachView(ch_1.Select(""))
      for res in ch_2.residues:
        ch2_aln[res.number.GetNum() - 1] = res.GetOneLetterCode()
      ch2_aln = "".join(ch2_aln)
      seq_2 = seq.CreateSequence(ch_2.name, str(ch2_aln))
      seq_2.AttachView(ch_2.Select(""))
      # Create alignment
      aln = seq.CreateAlignment()
      aln.AddSequence(seq_1)
      aln.AddSequence(seq_2)
    else:
      seq_1 = seq.SequenceFromChain(ch_1.name, ch_1)
      seq_2 = seq.SequenceFromChain(ch_2.name, ch_2)
      # align them
      aln = _AlignAtomSeqs(seq_1, seq_2)
    if aln:
      alns.append(aln)
  return alns

def _GetMappedResidues(alns):
  """
  :return: Mapping of shared residues in *alns* (with views attached)
           (see :attr:`QSscorer.mapped_residues`).
  :param alns: See :attr:`QSscorer.alignments`
  """
  mapped_residues = dict()
  for aln in alns:
    # prepare dict for c1
    c1 = aln.GetSequence(0).name
    mapped_residues[c1] = dict()
    # get shared residues
    v1, v2 = seq.ViewsFromAlignment(aln)
    for res_1, res_2 in zip(v1.residues, v2.residues):
      r1 = res_1.number.num
      r2 = res_2.number.num
      mapped_residues[c1][r1] = r2

  return mapped_residues

def _GetExtraWeights(contacts, done, mapped_residues):
  """Return sum of extra weights for contacts of chains in set and not in done.
  :return: Tuple (weight_extra_mapped, weight_extra_all).
           weight_extra_mapped only sums if both cX,rX in mapped_residues
           weight_extra_all sums all
  :param contacts: See :func:`GetContacts` for first entity
  :param done: List of (c1, c2, r1, r2) tuples to ignore
  :param mapped_residues: See :func:`_GetMappedResidues`
  """
  weight_extra_mapped = 0
  weight_extra_non_mapped = 0
  for c1 in contacts:
    for c2 in contacts[c1]:
      for r1 in contacts[c1][c2]:
        for r2 in contacts[c1][c2][r1]:
          if (c1, c2, r1, r2) not in done:
            weight = _weight(contacts[c1][c2][r1][r2])
            if     c1 in mapped_residues and r1 in mapped_residues[c1] \
               and c2 in mapped_residues and r2 in mapped_residues[c2]:
              weight_extra_mapped += weight
            else:
              weight_extra_non_mapped += weight
  return weight_extra_mapped, weight_extra_mapped + weight_extra_non_mapped

def _GetScores(contacts_1, contacts_2, mapped_residues, chain_mapping):
  """Get QS scores (see :class:`QSscorer`).

  Note that if some chains are not to be considered at all, they must be removed
  from *contacts_1* / *contacts_2* prior to calling this.

  :param contacts_1: See :func:`GetContacts` for first entity
  :param contacts_2: See :func:`GetContacts` for second entity
  :param mapped_residues: See :func:`_GetMappedResidues`
  :param chain_mapping: Maps any chain name in *mapped_residues* to chain name
                        for *contacts_2*.
  :type chain_mapping:  :class:`dict` (:class:`str` / :class:`str`)
  :return: Tuple (QS_best, QS_global) (see :attr:`QSscorer.best_score` and
           see :attr:`QSscorer.global_score`)
  """
  # keep track of considered contacts as set of (c1,c2,r1,r2) tuples
  done_1 = set()
  done_2 = set()
  weighted_scores = 0
  weight_sum = 0
  # naming cXY, rXY: X = 1/2 for contact, Y = 1/2/T for mapping (T = tmp)
  # -> d1 = contacts_1[c11][c21][r11][r21], d2 = contacts_2[c12][c22][r12][r22]
  for c11 in contacts_1:
    # map to other one
    if c11 not in mapped_residues: continue
    c1T = chain_mapping[c11]
    # second chain
    for c21 in contacts_1[c11]:
      # map to other one
      if c21 not in mapped_residues: continue
      c2T = chain_mapping[c21]
      # flip if needed
      flipped_chains = (c1T > c2T)
      if flipped_chains:
        c12, c22 = c2T, c1T
      else:
        c12, c22 = c1T, c2T
      # skip early if no contacts there
      if c12 not in contacts_2 or c22 not in contacts_2[c12]: continue
      # loop over res.num. in c11
      for r11 in contacts_1[c11][c21]:
        # map to other one and skip if no contacts there
        if r11 not in mapped_residues[c11]: continue
        r1T = mapped_residues[c11][r11]
        # loop over res.num. in c21
        for r21 in contacts_1[c11][c21][r11]:
          # map to other one and skip if no contacts there
          if r21 not in mapped_residues[c21]: continue
          r2T = mapped_residues[c21][r21]
          # flip if needed
          if flipped_chains:
            r12, r22 = r2T, r1T
          else:
            r12, r22 = r1T, r2T
          # skip early if no contacts there
          if r12 not in contacts_2[c12][c22]: continue
          if r22 not in contacts_2[c12][c22][r12]: continue
          # ok now we can finally do the scoring
          d1 = contacts_1[c11][c21][r11][r21]
          d2 = contacts_2[c12][c22][r12][r22]
          weight = _weight(min(d1, d2))
          weighted_scores += weight * (1 - abs(d1 - d2) / 12)
          weight_sum += weight
          # keep track of done ones
          done_1.add((c11, c21, r11, r21))
          done_2.add((c12, c22, r12, r22))

  LogVerbose("Shared contacts: %d" % len(done_1))

  # add the extra weights
  weights_extra_1 = _GetExtraWeights(contacts_1, done_1, mapped_residues)
  mapped_residues_2 = dict()
  for c1 in mapped_residues:
    c2 = chain_mapping[c1]
    mapped_residues_2[c2] = set()
    for r1 in mapped_residues[c1]:
      mapped_residues_2[c2].add(mapped_residues[c1][r1])
  weights_extra_2 = _GetExtraWeights(contacts_2, done_2, mapped_residues_2)
  weight_extra_mapped = weights_extra_1[0] + weights_extra_2[0]
  weight_extra_all    = weights_extra_1[1] + weights_extra_2[1]
  
  # get scores
  denom_best = weight_sum + weight_extra_mapped
  denom_all  = weight_sum + weight_extra_all
  if denom_best == 0:
    QS_best = 0
  else:
    QS_best = weighted_scores / denom_best
  if denom_all == 0:
    QS_global = 0
  else:
    QS_global = weighted_scores / denom_all
  return QS_best, QS_global

def _weight(dist):
  """
  This weight expresses the probability of two residues to interact given the CB
  distance (from Xu et al. http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2573399/)
  """
  if dist <= 5.0:
    return 1
  x = (dist-5.0)/4.28;
  return np.exp(-2*x*x)


def _GetQsSuperposition(alns):
  """
  :return: Superposition result based on shared CA atoms in *alns*
           (with views attached) (see :attr:`QSscorer.superposition`).
  :param alns: See :attr:`QSscorer.alignments`
  """
  # check input
  if not alns:
    raise QSscoreError('No successful alignments for superposition!')
  # prepare views
  view_1 = alns[0].GetSequence(0).attached_view.CreateEmptyView()
  view_2 = alns[0].GetSequence(1).attached_view.CreateEmptyView()
  # collect CA from alignments in proper order
  for aln in alns:
    for col in aln:
      res_1 = col.GetResidue(0)
      res_2 = col.GetResidue(1)
      if res_1.IsValid() and res_2.IsValid():
        ca_1 = res_1.FindAtom('CA')
        ca_2 = res_2.FindAtom('CA')
        if ca_1.IsValid() and ca_2.IsValid():
          view_1.AddAtom(ca_1)
          view_2.AddAtom(ca_2)
  # superpose them without chainging entities
  res = mol.alg.SuperposeSVD(view_1, view_2, apply_transform=False)
  return res


def _AddResidue(edi, res, rnum, chain, calpha_only):
  """
  Add residue *res* with res. num. *run* to given *chain* using editor *edi*.
  Either all atoms added or (if *calpha_only*) only CA.
  """
  if calpha_only:
    ca_atom = res.FindAtom('CA')
    if ca_atom.IsValid():
      new_res = edi.AppendResidue(chain, res.name, rnum)
      edi.InsertAtom(new_res, ca_atom.name, ca_atom.pos)
  else:
    new_res = edi.AppendResidue(chain, res.name, rnum)
    for atom in res.atoms:
      edi.InsertAtom(new_res, atom.name, atom.pos)

def _MergeAlignedChains(alns, ent_1, ent_2, calpha_only, penalize_extra_chains):
  """
  Create two new entities (based on the alignments attached views) where all
  residues have same numbering (when they're aligned) and they are all pushed to
  a single chain X. Also append extra chains contained in *ent_1* and *ent_2*
  but not contained in *alns*.

  Used for :attr:`QSscorer.lddt_ref` and :attr:`QSscorer.lddt_mdl`

  :param alns: List of alignments with attached views (first sequence: *ent_1*,
               second: *ent_2*). Residue number in single chain is column index
               of current alignment + sum of lengths of all previous alignments
               (order of alignments as in input list).
  :type alns:  See :attr:`QSscorer.alignments`
  :param ent_1: First entity to process.
  :type ent_1:  :class:`~ost.mol.EntityHandle`
  :param ent_2: Second entity to process.
  :type ent_2:  :class:`~ost.mol.EntityHandle`
  :param calpha_only: If True, we only include CA atoms instead of all.
  :type calpha_only:  :class:`bool`
  :param penalize_extra_chains: If True, extra chains are added to model and
                                reference. Otherwise, only mapped ones.
  :type penalize_extra_chains:  :class:`bool`

  :return: Tuple of two single chain entities (from *ent_1* and from *ent_2*)
  :rtype:  :class:`tuple` of :class:`~ost.mol.EntityHandle`
  """
  # first new entity
  ent_ren_1 = mol.CreateEntity()
  ed_1 = ent_ren_1.EditXCS()
  new_chain_1 = ed_1.InsertChain('X')
  # second one
  ent_ren_2 = mol.CreateEntity()
  ed_2 = ent_ren_2.EditXCS()
  new_chain_2 = ed_2.InsertChain('X')
  # the alignment already contains sorted chains
  rnum = 0
  chain_done_1 = set()
  chain_done_2 = set()
  for aln in alns:
    chain_done_1.add(aln.GetSequence(0).name)
    chain_done_2.add(aln.GetSequence(1).name)
    for col in aln:
      rnum += 1
      # add valid residues to single chain entities
      res_1 = col.GetResidue(0)
      if res_1.IsValid():
        _AddResidue(ed_1, res_1, rnum, new_chain_1, calpha_only)
      res_2 = col.GetResidue(1)
      if res_2.IsValid():
        _AddResidue(ed_2, res_2, rnum, new_chain_2, calpha_only)
  # extra chains?
  if penalize_extra_chains:
    for chain in ent_1.chains:
      if chain.name in chain_done_1:
        continue
      for res in chain.residues:
        rnum += 1
        _AddResidue(ed_1, res, rnum, new_chain_1, calpha_only)
    for chain in ent_2.chains:
      if chain.name in chain_done_2:
        continue
      for res in chain.residues:
        rnum += 1
        _AddResidue(ed_2, res, rnum, new_chain_2, calpha_only)
  # get entity names
  ent_ren_1.SetName(aln.GetSequence(0).GetAttachedView().GetName())
  ent_ren_2.SetName(aln.GetSequence(1).GetAttachedView().GetName())
  # connect atoms
  if not conop.GetDefaultLib():
    raise RuntimeError("QSscore computation requires a compound library!")
  pr = conop.RuleBasedProcessor(conop.GetDefaultLib())
  pr.Process(ent_ren_1)  
  ed_1.UpdateICS()
  pr.Process(ent_ren_2)  
  ed_2.UpdateICS()
  return ent_ren_1, ent_ren_2


# specify public interface
__all__ = ('QSscoreError', 'QSscorer', 'QSscoreEntity', 'FilterContacts',
           'GetContacts', 'OligoLDDTScorer', 'MappedLDDTScorer')
