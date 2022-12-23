import os
from ost import mol
from ost.mol.alg import chain_mapping
import numpy as np


class LigandScorer:
    """ Helper class to access the various small molecule ligand (non polymer)
    scores available from ost.mol.alg.

    Mostly expects cleaned up structures (you can use the
    `~ost.mol.alg.scoring.Scorer` outputs for that).

    :param model: Model structure - a deep copy is available as :attr:`model`.
                  No additional processing (ie. Molck), checks,
                  stereochemistry checks or sanitization is performed on the
                  input.
    :type model: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param target: Target structure - a deep copy is available as :attr:`target`.
                  No additional processing (ie. Molck), checks or sanitization
                  is performed on the input.
    :type target: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param model_ligands: Model ligands, either a :class:list of
                  :class:`ost.mol.ResidueHandle`/:class:`ost.mol.ResidueView`
                  or of :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
                  containing a single residue each. If `None`, ligands will be
                  extracted from the `model` entity, from chains with
                  :class:`~ost.mol.ChainType` `CHAINTYPE_NON_POLY` (this is
                  normally set properly in entities loaded from mmCIF).
    :type model_ligands: :class:`list`
    :param target_ligands: Target ligands, either a :class:list of
                  :class:`ost.mol.ResidueHandle`/:class:`ost.mol.ResidueView`
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
    :param chain_mapper: a chain mapper initialized for the target structure.
                         If None (default), a chain mapper will be initialized
                         lazily as required.
    :type chain_mapper:  :class:`ost.mol.alg.chain_mapping.ChainMapper`


    """
    def __init__(self, model, target, model_ligands=None, target_ligands=None,
                 resnum_alignments=False, chain_mapper=None):

        if isinstance(model, mol.EntityView):
            self._model = mol.CreateEntityFromView(model, False)
        elif isinstance(model, mol.EntityHandle):
            self._model = model.Copy()
        else:
            raise RuntimeError("model must be of type EntityView/EntityHandle")

        if isinstance(target, mol.EntityView):
            self._target = mol.CreateEntityFromView(target, False)
        elif isinstance(target, mol.EntityHandle):
            self._target = target.Copy()
        else:
            raise RuntimeError("model must be of type EntityView/EntityHandle")

        # Extract ligands from target
        if target_ligands is None:
            self.target_ligands = self._extract_ligands(self._target)
        else:
            # TODO: sanitize given ligands
            self.target_ligands = target_ligands

        # Extract ligands from model
        if model_ligands is None:
            self.model_ligands = self._extract_ligands(self._model)
        else:
            # TODO: sanitize given ligands
            self.model_ligands = target_ligands

        self._chain_mapper = chain_mapper
        self.resnum_alignments = resnum_alignments

        # lazily computed scores
        self._lddt_pli = None
        self._rmsd = None
        self._lddt_bs = None

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
                    extracted_ligands.append(residue)
        return extracted_ligands

