from ost import mol
from ost import seq
from ost import io
from ost.mol.alg import lddt
from ost.mol.alg import chain_mapping

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
