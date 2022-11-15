import os
import json

import numpy as np

from ost import geom
from ost import mol


def _PotentialDisulfid(a_one, a_two):
    """ Returns whether two atoms can potentially build a disulfid bond

    Assumes that they're from two distinct residues
    """
    if a_one.GetName() == "SG" and a_two.GetName() == "SG":
        if a_one.GetResidue().GetName() == "CYS":
            if a_two.GetResidue().GetName() == "CYS":
                return True
    return False


def _GetAngles(bonds):
    """ Returns list of angles based on bonds

    Returns list of tuples, each tuple has three atom handles
    representing angles
    """
    angles = list()
    done = set()
    for bond in bonds:
        h1 = bond.first.GetHashCode()
        h2 = bond.second.GetHashCode()
        for a in bond.first.GetBondPartners():
            h0 = a.GetHashCode()
            if h0 != h2:
                if ((h0, h1, h2)) not in done and (h2, h1, h0) not in done:
                    angles.append((a, bond.first, bond.second))
                    done.add((h0, h1, h2))
        for a in bond.second.GetBondPartners():
            h3 = a.GetHashCode()
            if h3 != h1:
                if ((h1, h2, h3)) not in done and (h3, h2, h1) not in done:
                    angles.append((bond.first, bond.second, a))
                    done.add((h1, h2, h3))
    return angles


def _ParseBondData(doc):
    """ Parse stereochemistry data for bonds

    That is expected distances and standard deviations from a
    :class:`gemmi.Document`. Concatenates results form all loops with tags:
    _chem_comp_bond.comp_id, _chem_comp_bond.atom_id_1,
    _chem_comp_bond.atom_id_2, _chem_comp_bond.value_dist,
    _chem_comp_bond.value_dist_esd

    :param doc: Gemmi doc representing cif file opened with
                gemmi.cif.read_file(filepath)
    :type doc: :class:`gemmi.Document`
    :returns: :class:`dict` with one key per compound, the respective value
              is again a dict with key f"{at_1}_{at_2}" and value
              [dist, dist_std].
    """
    data = dict()
    for block in doc:
        comp_id = block.find_values("_chem_comp_bond.comp_id")
        at_1 = block.find_values("_chem_comp_bond.atom_id_1")
        at_2 = block.find_values("_chem_comp_bond.atom_id_2")
        dist = block.find_values("_chem_comp_bond.value_dist")
        dist_std = block.find_values("_chem_comp_bond.value_dist_esd")
        if None not in [comp_id, at_1, at_2, dist, dist_std]:
            for a, b, c, d, e in zip(comp_id, at_1, at_2, dist, dist_std):
                if a not in data:
                    data[a] = dict()
                data[a][f"{b}_{c}"] = [float(d), float(e)]
    return data


def _ParseAngleData(doc):
    """ Parse stereochemistry data for angles

    That is expected distances and standard deviations from a
    :class:`gemmi.Document`. Concatenates results form all loops with tags:
    _chem_comp_angle.comp_id, _chem_comp_angle.atom_id_1,
    _chem_comp_angle.atom_id_2, _chem_comp_angle.atom_id_2,
    _chem_comp_angle.value_angle, _chem_comp_angle.value_angle_esd

    :param doc: Gemmi doc representing cif file opened with
                gemmi.cif.read_file(filepath)
    :type doc: :class:`gemmi.Document`
    :returns: :class:`dict` with one key per compound, the respective value
              is again a dict with key f"{at_1}_{at_2}_{at_3}" and value
              [angle, angle_std].
    """
    data = dict()
    for block in doc:
        comp_id = block.find_values("_chem_comp_angle.comp_id")
        at_1 = block.find_values("_chem_comp_angle.atom_id_1")
        at_2 = block.find_values("_chem_comp_angle.atom_id_2")
        at_3 = block.find_values("_chem_comp_angle.atom_id_3")
        angle = block.find_values("_chem_comp_angle.value_angle")
        angle_std = block.find_values("_chem_comp_angle.value_angle_esd")
        if None not in [comp_id, at_1, at_2, at_3, angle, angle_std]:
            for a, b, c, d, e, f in zip(comp_id, at_1, at_2, at_3, angle,
                                        angle_std):
                if a not in data:
                    data[a] = dict()
                data[a][f"{b}_{c}_{d}"] = [float(e), float(f)]
    return data


def StereoDataFromMON_LIB(mon_lib_path, compounds=None):
    """ Parses stereochemistry parameters from CCP4 MON_LIB

    CCP4 `MON_LIB <https://www.ccp4.ac.uk/html/mon_lib.html>`_ contains
    data on ideal bond lengths/angles for compounds.

    Original data (several updates in the meantime) come from:

    * Amino acid bond lengths and angles: Engh and Huber, Acta Cryst.
      A47, 392-400 (1991).
    * Purine and pyrimidine bond lengths and angles: O. Kennard & R. Taylor
      (1982), J. Am. Soc. Chem. vol. 104, pp. 3209-3212.
    * Sugar-phosphate backbone bond lengths and bond angles: W. Saenger’s
      Principles of Nucleic Acid Structure (1983), Springer-Verlag, pp. 70,86.

    This function adds a dependency to the
    `gemmi <https://github.com/project-gemmi/gemmi/>`_ library to read cif
    files.

    :param mon_lib_path: Path to CCP4 MON_LIB
    :type mon_lib_path: :class:`str`
    :param compounds: Compounds to parse - parses proteinogenic amino acids
                      and nucleotides if not given.
    :type compounds: :class:`list`
    :returns: :class:`dict` with stereochemistry parameters
    """
    if compounds is None:
        compounds = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                     'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'MSE', 'PHE', 'PRO',
                     'SER', 'THR', 'TRP', 'TYR', 'VAL', 'DA', 'A', 'DC', 'C',
                     'DG', 'G', 'DU', 'U', 'DT', 'DI', 'I']

    cif_paths = list()
    for c in compounds:
        p = os.path.join(mon_lib_path, c[0].lower(), c + ".cif")
        if not os.path.exists(p):
            raise RuntimeError(f"Tried to find cif file for compound {c} "
                               f"in specified MON_LIB ({mon_lib_path})."
                               f"Expected file ({p}) does not exist.")
        cif_paths.append(p)

    # hide import to avoid it as dependency for the whole module
    from gemmi import cif
    # construct return dict from first element and subsequently 
    # add the remainder
    doc = cif.read_file(cif_paths[0])
    data = {"bond_data": _ParseBondData(doc),
            "angle_data": _ParseAngleData(doc)}
    for cp in cif_paths[1:]:
        doc = cif.read_file(cp)
        bond_data = _ParseBondData(doc)
        angle_data = _ParseAngleData(doc)
        data["bond_data"].update(bond_data)
        data["angle_data"].update(angle_data)
    return data


def GetBondParam(a1, a2, stereo_data):
    """ Returns mean and standard deviation for bond

    :param a1: First atom that defines bond
    :type a1: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param a2: Second atom that defines bond
    :type a2: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param stereo_data: Stereochemistry data
    :type stereo_data: :class:`dict`
    :returns: :class:`tuple` with mean and standard deviation. Values are None
              if respective bond is not found in *stereo_data*
    """
    if a1.GetResidue().GetHashCode() == a2.GetResidue().GetHashCode():
        # intra residue case, inter-residue case not yet implemented
        rname = a1.GetResidue().GetName()
        a1name = a1.GetName()
        a2name = a2.GetName()
        if rname in stereo_data["bond_data"]:
            key = a1name + "_" + a2name
            if key in stereo_data["bond_data"][rname]:
                mean = stereo_data["bond_data"][rname][key][0]
                std = stereo_data["bond_data"][rname][key][1]
                return (mean, std)
            key = a2name + "_" + a1name
            if key in stereo_data["bond_data"][rname]:
                mean = stereo_data["bond_data"][rname][key][0]
                std = stereo_data["bond_data"][rname][key][1]
                return (mean, std)
    return (None, None)


def GetAngleParam(a1, a2, a3, stereo_data):
    """ Returns mean and standard deviation for angle

    :param a1: First atom that defines angle
    :type a1: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param a2: Second atom that defines angle
    :type a2: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param a3: Third atom that defines angle
    :type a3: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param stereo_data: Stereochemistry data
    :type stereo_data: :class:`dict`
    :returns: :class:`tuple` with mean and standard deviation. Values are None
              if respective angle is not found in *stereo_data*
    """
    h1 = a1.GetResidue().handle.GetHashCode()
    h2 = a2.GetResidue().handle.GetHashCode()
    h3 = a3.GetResidue().handle.GetHashCode()
    if h1 == h2 and h2 == h3:
        # intra residue case, inter-residue case not yet implemented
        rname = a1.GetResidue().GetName()
        a1name = a1.GetName()
        a2name = a2.GetName()
        a3name = a3.GetName()
        if rname in stereo_data["angle_data"]:
            key = a1name + "_" + a2name + "_" + a3name
            if key in stereo_data["angle_data"][rname]:
                mean = stereo_data["angle_data"][rname][key][0]
                std = stereo_data["angle_data"][rname][key][1]
                return (mean, std)
            key = a3name + "_" + a2name + "_" + a1name
            if key in stereo_data["angle_data"][rname]:
                mean = stereo_data["angle_data"][rname][key][0]
                std = stereo_data["angle_data"][rname][key][1]
                return (mean, std)
    return (None, None)


def GetClashes(ent, vdw_radii = None, tolerance = 1.5, disulfid_dist = 2.03,
               disulfid_tolerance = 1.0):
    """ Identifies clashing atoms

    A clash between two non-bonded atoms is defined as their distance d being
    below the sum of their vdw radii with some subtracted tolerance value.

    The default values are not very sensitive.

    :param ent: Entity for which you want to identify clashing atoms
    :type ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param vdw_radii: Element based van der Waals radii. Only atoms of these
                      elements will be considered. If not given, default values
                      for all elements occuring in proteins/nucleotides are
                      used. Must be provided as :class:`dict`, where they key
                      are elements (capitalized) and value the respective radii
                      in Angstrom.
    :type vdw_radii: :class:`dict`
    :param tolerance: Tolerance value
    :param disulfid_dist: Summed vdw radius that is used if two Sulfurs that can
                          potentially build a disulfid bond interact
    :type disulfid_dist: :class:`float`
    :param disulfid_tolerance: The respective tolerance
    :type disulfid_dist: :class:`float`
    :returns: A :class:`list` of pairs. Each pair consists of two
              :class:`ost.mol.AtomView` from *ent* that are clashing.
    """

    if vdw_radii is None:
        vdw_radii = {"C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80}

    for ele in vdw_radii.keys():
        if ele.upper() != ele:
            raise RuntimeError(f"Elements in vdw_radii must be upper case. "
                               f"Got {ele}")

    # it would be elegant to just do a selection by the ele property. However,
    # thats case sensitive so someone could define a vdw radius for Cl but
    # the element of the atom is CL.
    elements = set([ele.upper() for ele in vdw_radii.keys()])
    for a in ent.atoms:
        if a.GetElement().upper() in elements:
            a.SetIntProp("clash_check", 1)
    clash_ent = ent.Select("gaclash_check:0=1")

    max_radius = max(vdw_radii.values())
    max_radius = max(max_radius, 0.5*disulfid_dist)
    min_tolerance = min(tolerance, disulfid_tolerance)
    radius = 2*max_radius-min_tolerance

    done = set()
    return_list = list()
    for a in clash_ent.atoms:
        a_hash = a.handle.GetHashCode()
        close_atoms = clash_ent.FindWithin(a.GetPos(), radius)
        for ca in close_atoms:
            ca_hash = ca.handle.GetHashCode()
            if a_hash != ca_hash and not mol.BondExists(a.handle, ca.handle):
                d = geom.Distance(a.GetPos(), ca.GetPos())
                if _PotentialDisulfid(a, ca):
                    thresh = disulfid_dist - disulfid_tolerance
                else:
                    thresh = vdw_radii[a.GetElement().upper()]
                    thresh += vdw_radii[ca.GetElement().upper()]
                    thresh -= tolerance
                if d < thresh:
                    # check if already there, add if not
                    hash_pair = (min(a_hash, ca_hash), max(a_hash, ca_hash))
                    if hash_pair not in done:
                        done.add(hash_pair)
                        return_list.append((a, ca))
    return return_list


def GetBadBonds(ent, stereo_data = None, tolerance=12):
    """ Identify unrealistic bonds

    :param ent: Entity for which you want to identify unrealistic bonds
    :type ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param stereo_data: Stereochemistry data
    :type stereo_data: :class:`dict`
    :param tolerance: Bonds that devaiate more than *tolerance* times standard
                      deviation from expected mean are considered bad
    :type tolerance: :class:`int`
    :returns: :class:`list` of pairs. Each pair consists of two
              :class:`ost.mol.AtomHandle` from *ent* that represent bad bonds.

    """
    assert("bond_data" in stereo_data)
    return_list = list()
    for b in ent.bonds:
        a1 = b.first
        a2 = b.second
        mean, std = GetBondParam(a1, a2, stereo_data)
        if None not in [mean, std]:
            diff = abs(mean-b.length)
            if diff > tolerance*std:
                return_list.append((a1, a2))
    return return_list


def GetBadAngles(ent, stereo_data = None, tolerance=12):
    """ Identify unrealistic angles

    :param ent: Entity for which you want to identify unrealistic angles
    :type ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param stereo_data: Stereochemistry data
    :type stereo_data: :class:`dict`
    :param tolerance: Angles that devaiate more than *tolerance* times standard
                      deviation from expected mean are considered bad
    :type tolerance: :class:`int`
    :returns: :class:`list` of tuples. Each tuple consists of three
              :class:`ost.mol.AtomHandle` from *ent* that represent bad angles.
    """
    assert("angle_data" in stereo_data)
    return_list = list()
    for a in _GetAngles(ent.bonds):
        mean, std = GetAngleParam(a[0], a[1], a[2], stereo_data)
        if None not in [mean, std]:
            angle = geom.Angle(a[0].GetPos() - a[1].GetPos(),
                               a[2].GetPos() - a[1].GetPos())
            angle = angle/np.pi*180 # stereo params are in degrees
            diff = abs(mean-angle)
            if diff > tolerance*std:
                return_list.append(a)
    return return_list
