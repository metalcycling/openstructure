import os
import json
import datetime

import numpy as np

import ost
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


def _GetResidueType(atoms):
    """ Identifies type in StereoLinkData

    :param atoms: Atoms that define a bond or angle
    :type atoms: :class:`list` of :class:`AtomHandle`
    :returns: :class:`str` with which the respective parameters can be
              accessed in default stereo link data, None if no match is found
    """
    residues = [a.GetResidue().handle for a in atoms]
    chem_types = list(set([str(r.GetChemType()) for r in residues]))

    if len(chem_types) == 1 and chem_types[0] == 'N':
        return "NA"
    elif len(chem_types) == 1 and chem_types[0] == 'A':
        # in both cases, bond or angle, there should be exactly two residues
        # involved
        tmp = list()
        r_hashes = set()
        for r in residues:
            h = r.GetHashCode()
            if h not in r_hashes:
                r_hashes.add(h)
                tmp.append(r)
        residues = tmp
        if len(residues) != 2:
            return None

        # need to be sorted
        if residues[0].GetNumber() > residues[1].GetNumber():
            r0 = residues[1]
            r1 = residues[0]
        else:
            r0 = residues[0]
            r1 = residues[1]

        if r1.GetName() == "GLY":
            return "GLY"
        elif r1.GetName() == "PRO":
            a = r0.FindAtom("CA")
            b = r0.FindAtom("C")
            c = r1.FindAtom("N")
            d = r1.FindAtom("CA")
            if a.IsValid() and b.IsValid() and c.IsValid() and d.IsValid():
                omega = geom.DihedralAngle(a.GetPos(), b.GetPos(),
                                           c.GetPos(), d.GetPos())
                if abs(omega) < 1.57:
                    return "PRO_CIS"
                else:
                    return "PRO_TRANS"
        else:
            return "PEPTIDE"

    return None


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
                key = '_'.join([b.strip('\"'), c.strip('\"')])
                data[a][key] = [float(d), float(e)]
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
                key = '_'.join([b.strip('\"'), c.strip('\"'), d.strip('\"')])
                data[a][key] = [float(e), float(f)]
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

    # add license info
    copying_str = f"This data has been derived from the CCP4 MON_LIB on "
    copying_str += f"{datetime.datetime.now()}. MON_LIB is licensed under "
    copying_str += f"GNU LESSER GENERAL PUBLIC LICENSE Version 3. Consult the "
    copying_str += f"latest CCP4 for the full license text."
    data["COPYING"] = copying_str

    return data


def GetBondParam(a1, a2, stereo_data = None, stereo_link_data = None):
    """ Returns mean and standard deviation for bond

    :param a1: First atom that defines bond
    :type a1: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param a2: Second atom that defines bond
    :type a2: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param stereo_data: Stereochemistry data, use return value of
                        :func:`GetDefaultStereoData` if not given.
                        If you call this function repeatedly, you
                        really should provide *stereo_data*!
    :type stereo_data: :class:`dict`
    :param stereo_link_data: Stereochemistry data, use return value of
                             :func:`GetDefaultStereoLinkData` if not given.
                             If you call this function repeatedly, you
                             really should provide *stereo_link_data*!
    :type stereo_link_data: :class:`dict`
    :returns: :class:`tuple` with mean and standard deviation. Values are None
              if respective bond is not found in *stereo_data*
    """
    if stereo_data is None:
        stereo_data = GetDefaultStereoData()
    if stereo_link_data is None:
        stereo_link_data = GetDefaultStereoLinkData()

    residue_data = None
    if a1.GetResidue().GetHashCode() == a2.GetResidue().GetHashCode():
        # intra residue case
        rname = a1.GetResidue().GetName()
        if rname in stereo_data["bond_data"]:
            residue_data = stereo_data["bond_data"][rname]
    else:
        # inter residue case
        residue_type = _GetResidueType([a1, a2])
        if residue_type is not None:
            residue_data = stereo_link_data["bond_data"][residue_type]

    if residue_data is not None:
        a1name = a1.GetName()
        a2name = a2.GetName()
        key = a1name + "_" + a2name
        if key in residue_data:
            return (residue_data[key][0], residue_data[key][1])
        key = a2name + "_" + a1name
        if key in residue_data:
            return (residue_data[key][0], residue_data[key][1])

    return (None, None)


def GetAngleParam(a1, a2, a3, stereo_data = None, stereo_link_data = None):
    """ Returns mean and standard deviation for angle

    :param a1: First atom that defines angle
    :type a1: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param a2: Second atom that defines angle
    :type a2: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param a3: Third atom that defines angle
    :type a3: :class:`ost.mol.AtomView`/:class:`ost.mol.AtomHandle`
    :param stereo_data: Stereochemistry data, use return value of
                        :func:`GetDefaultStereoData` if not given.
                        If you call this function repeatedly, you
                        really should provide *stereo_data*!
    :type stereo_data: :class:`dict`
    :param stereo_link_data: Stereochemistry data, use return value of
                             :func:`GetDefaultStereoLinkData` if not given.
                             If you call this function repeatedly, you
                             really should provide *stereo_link_data*!
    :type stereo_link_data: :class:`dict`
    :returns: :class:`tuple` with mean and standard deviation. Values are None
              if respective angle is not found in *stereo_data*
    """
    if stereo_data is None:
        stereo_data = GetDefaultStereoData()
    if stereo_link_data is None:
        stereo_link_data = GetDefaultStereoLinkData()
    h1 = a1.GetResidue().handle.GetHashCode()
    h2 = a2.GetResidue().handle.GetHashCode()
    h3 = a3.GetResidue().handle.GetHashCode()
    residue_data = None
    if h1 == h2 and h2 == h3:
        # intra residue case
        rname = a1.GetResidue().GetName()
        if rname in stereo_data["angle_data"]:
            residue_data = stereo_data["angle_data"][rname]
    else:
        # inter residue case
        residue_type = _GetResidueType([a1, a2, a3])
        if residue_type in stereo_link_data["angle_data"]:
            residue_data = stereo_link_data["angle_data"][residue_type]

    if residue_data is not None:
        a1name = a1.GetName()
        a2name = a2.GetName()
        a3name = a3.GetName()
        key = a1name + "_" + a2name + "_" + a3name
        if key in residue_data:
            return (residue_data[key][0], residue_data[key][1])
        key = a3name + "_" + a2name + "_" + a1name
        if key in residue_data:
            return (residue_data[key][0], residue_data[key][1])
    return (None, None)


class ClashInfo:
    """ Object to hold info on clash

    Constructor arguments are available as attributes:

    * a1
    * a2
    * dist
    * tolerated_dist
    """
    def __init__(self, a1, a2, dist, tolerated_dist):
        self.a1 = a1
        self.a2 = a2
        self.dist = dist
        self.tolerated_dist = tolerated_dist

    def ToJSON(self):
        """ Return JSON serializable dict
        """
        return {"a1": self.a1.GetQualifiedName(),
                "a2": self.a2.GetQualifiedName(),
                "dist": self.dist,
                "tolerated_dist": self.tolerated_dist}


class BondViolationInfo:
    """ Object to hold info on bond violation

    Constructor arguments are available as attributes:

    * a1
    * a2
    * length
    * exp_length
    * std
    """
    def __init__(self, a1, a2, length, exp_length, std):
        self.a1 = a1
        self.a2 = a2
        self.length = length
        self.exp_length = exp_length
        self.std = std

    def ToJSON(self):
        """ Return JSON serializable dict
        """
        return {"a1": self.a1.GetQualifiedName(),
                "a2": self.a2.GetQualifiedName(),
                "length": self.length,
                "exp_length": self.exp_length,
                "std": self.std}


class AngleViolationInfo:
    """ Object to hold info on angle violation

    Constructor arguments are available as attributes:

    * a1
    * a2
    * a3
    * angle
    * exp_angle
    * std
    """
    def __init__(self, a1, a2, a3, angle, exp_angle, std):
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.angle = angle
        self.exp_angle = exp_angle
        self.std = std

    def ToJSON(self):
        """ Return JSON serializable dict
        """
        return {"a1": self.a1.GetQualifiedName(),
                "a2": self.a2.GetQualifiedName(),
                "a3": self.a3.GetQualifiedName(),
                "angle": self.angle,
                "exp_angle": self.exp_angle,
                "std": self.std}


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
    :returns: A :class:`list` of :class:`ClashInfo`
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
                        return_list.append(ClashInfo(a.handle, ca.handle, d,
                                                     thresh))
    return return_list


def GetBadBonds(ent, stereo_data = None, stereo_link_data = None, tolerance=12):
    """ Identify unrealistic bonds

    :param ent: Entity for which you want to identify unrealistic bonds
    :type ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param stereo_data: Stereochemistry data, use return value of
                        :func:`GetDefaultStereoData` if not given.
    :type stereo_data: :class:`dict`
    :param stereo_link_data: Stereochemistry data, use return value of
                             :func:`GetDefaultStereoLinkData` if not given.
    :type stereo_link_data: :class:`dict`
    :param tolerance: Bonds that devaiate more than *tolerance* times standard
                      deviation from expected mean are considered bad
    :type tolerance: :class:`int`
    :returns: :class:`list` :class:`BondViolationInfo`

    """
    if stereo_data is None:
        stereo_data = GetDefaultStereoData()
    if stereo_link_data is None:
        stereo_link_data = GetDefaultStereoLinkData()
    return_list = list()
    for b in ent.bonds:
        a1 = b.first
        a2 = b.second
        mean, std = GetBondParam(a1, a2, stereo_data = stereo_data,
                                 stereo_link_data = stereo_link_data)
        if None not in [mean, std]:
            l = b.length
            if abs(mean-l) > tolerance*std:
                return_list.append(BondViolationInfo(a1, a2, l, mean, std))
    return return_list


def GetBadAngles(ent, stereo_data = None, stereo_link_data = None,
                 tolerance = 12):
    """ Identify unrealistic angles

    :param ent: Entity for which you want to identify unrealistic angles
    :type ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param stereo_data: Stereochemistry data, use return value of
                        :func:`GetDefaultStereoData` if not given.
    :type stereo_data: :class:`dict`
    :param stereo_link_data: Stereochemistry data, use return value of
                             :func:`GetDefaultStereoLinkData` if not given.
    :type stereo_link_data: :class:`dict`
    :param tolerance: Angles that devaiate more than *tolerance* times standard
                      deviation from expected mean are considered bad
    :type tolerance: :class:`int`
    :returns: :class:`list` of :class:`AngleViolationInfo`
    """
    if stereo_data is None:
        stereo_data = GetDefaultStereoData()
    if stereo_link_data is None:
        stereo_link_data = GetDefaultStereoLinkData()
    return_list = list()
    for a in _GetAngles(ent.bonds):
        mean, std = GetAngleParam(a[0], a[1], a[2], stereo_data = stereo_data,
                                  stereo_link_data = stereo_link_data)
        if None not in [mean, std]:
            angle = geom.Angle(a[0].GetPos() - a[1].GetPos(),
                               a[2].GetPos() - a[1].GetPos())
            angle = angle/np.pi*180 # stereo params are in degrees
            diff = abs(mean-angle)
            if diff > tolerance*std:
                return_list.append(AngleViolationInfo(a[0], a[1], a[2], angle,
                                                      mean, std))
    return return_list


def StereoCheck(ent, stereo_data = None, stereo_link_data = None):
    """ Remove atoms with stereochemical problems

    Selects for peptide/nucleotides and calls :func:`GetClashes`,
    :func:`GetBadBonds` and :func:`GetBadAngles` with default
    parameters.

    * Amino acids: Remove full residue if backbone atom is involved in
      stereochemistry issue ("N", "CA", "C", "O"). Remove sidechain if any of
      the sidechain atoms is involved in stereochemistry issues.
    * Nucleotides: Remove full residue if backbone atom is involved in
      stereochemistry issue ("P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'",
      "C3'", "C2'", "C1'", "O4'", "O3'", "O2'"). Remove sidechain (base) if any
      of the sidechain atoms is involved in stereochemistry issues.

    :param ent: Entity to be stereochecked
    :type ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param stereo_data: Stereochemistry data, use return value of
                        :func:`GetDefaultStereoData` if not given.
    :type stereo_data: :class:`dict`
    :param stereo_link_data: Stereochemistry data, use return value of
                             :func:`GetDefaultStereoLinkData` if not given.
    :type stereo_link_data: :class:`dict`
    :returns: Tuple with four elements: 1) :class:`ost.mol.EntityView` of
              *ent* processed as described above 2) Return value of
              :func:`GetClashes` 3) return value of :func:`GetBadBonds`
              4) return value of :func:`GetBadAngles`
    """
    if stereo_data is None:
        stereo_data = GetDefaultStereoData()

    sel = ent.Select("peptide=true or nucleotide=true")
    clashes = GetClashes(sel)
    bad_bonds = GetBadBonds(sel, stereo_data = stereo_data)
    bad_angles = GetBadAngles(sel, stereo_data = stereo_data)

    # set stereo problems as properties on an atom level
    for clash in clashes:
        clash.a1.SetIntProp("stereo_problem", 1)
        clash.a2.SetIntProp("stereo_problem", 1)

    for bond in bad_bonds:
        bond.a1.SetIntProp("stereo_problem", 1)
        bond.a2.SetIntProp("stereo_problem", 1)

    for angle in bad_angles:
        angle.a1.SetIntProp("stereo_problem", 1)
        angle.a2.SetIntProp("stereo_problem", 1)
        angle.a3.SetIntProp("stereo_problem", 1)

    # set stereo problems as properties on a residue level
    bad_ent = ent.Select("gastereo_problem:0=1")
    if len(bad_ent.residues) > 0:
        pep_bb = set(["N", "CA", "C", "O"])
        nuc_bb = set(["P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "C3'",
                      "C2'", "C1'", "O4'", "O3'", "O2'"])

        for r in bad_ent.residues:
            bad_atoms = set([a.GetName() for a in r.atoms])
            r.SetIntProp("stereo_problem", 1)
            if r.GetChemType() == mol.ChemType.NUCLEOTIDES:
                if len(nuc_bb.intersection(bad_atoms)) > 0:
                    r.SetIntProp("stereo_problem_bb", 1)
            elif r.GetChemType() == mol.ChemType.AMINOACIDS:
                if len(pep_bb.intersection(bad_atoms)) > 0:
                    r.SetIntProp("stereo_problem_bb", 1)

        # explicitely add " as OpenStructure query language would not
        # understand ' otherwise
        nuc_bb = [f"\"{name}\"" for name in nuc_bb]

        pep_query = f"(peptide=true and grstereo_problem:0=0) or "
        pep_query += f"(peptide=true and grstereo_problem_bb:0=0 and "
        pep_query += f"aname={','.join(pep_bb)})"
        nuc_query = f"(nucleotide=true and grstereo_problem:0=0) or "
        nuc_query += f"(nucleotide=true and grstereo_problem_bb:0=0 and "
        nuc_query += f"aname={','.join(nuc_bb)})"
        query = pep_query + " or " + nuc_query
        return_view = sel.Select(query)
    else:
        return_view = sel

    return return_view, clashes, bad_bonds, bad_angles


def GetDefaultStereoData():
    """ Get default stereo data derived from CCP4 MON_LIB

    Used as default if not provided in :func:`GetBadBonds`, :func:`GetBadAngles`
    and :func:`StereoCheck`.

    MON_LIB is licensed under GNU LESSER GENERAL PUBLIC LICENSE Version 3.
    Consult the latest CCP4 for the full license text.
    """
    data_path = os.path.join(ost.GetSharedDataPath(), "stereo_data.json")
    with open(data_path, 'r') as fh:
        return json.load(fh)


def GetDefaultStereoLinkData():
    """ Get default stereo data for links between compounds

    Hardcoded from arbitrary sources, see comments in the code.
    
    :returns: Data for peptide bonds, nucleotide links and disulfid bonds that
              are used as default if not provided in :func:`GetBadBonds`,
              :func:`GetBadAngles` and :func:`StereoCheck`.
    """
    data = {"bond_data": dict(),
            "angle_data": dict()}

    # data for nucleotides - deliberately stolen from
    # geostd (https://github.com/phenix-project/geostd) which is basically
    # the Phenix equivalent for MON_LIB
    # used file: $GEOSTD_DIR/rna_dna/chain_link_rna2p.cif
    # Reason to not use the same data origin as peptides is that in CCP4
    # there is a bit a more fine grained differentiation of NA types
    # which makes things more complicated.
    data["bond_data"]["NA"] = dict()
    data["bond_data"]["NA"]["O3'_P"] = [1.607, 0.015]

    data["angle_data"]["NA"] = dict()
    data["angle_data"]["NA"]["O3'_P_O5'"] = [104.000, 1.500]
    data["angle_data"]["NA"]["O3'_P_OP1"] = [108.000, 3.000]
    data["angle_data"]["NA"]["O3'_P_OP2"] = [108.000, 3.000]
    data["angle_data"]["NA"]["C3'_O3'_P"] = [120.200, 1.500]

    # data for peptides - deliberately stolen from standard_geometry.cif file
    # which is shipped with CCP4
    # (_standard_geometry.version "Fri Feb 22 17:25:15 GMT 2013").
    data["bond_data"]["PEPTIDE"] = dict()
    data["bond_data"]["PEPTIDE"]["C_N"] = [1.336, 0.023]
    data["bond_data"]["PEPTIDE"]["SG_SG"] = [2.033, 0.016]

    data["bond_data"]["GLY"] = dict()
    data["bond_data"]["GLY"]["C_N"] = [1.326, 0.018]

    data["bond_data"]["PRO_CIS"] = dict()
    data["bond_data"]["PRO_CIS"]["C_N"] = [1.338, 0.019]
    data["bond_data"]["PRO_TRANS"] = dict()
    data["bond_data"]["PRO_TRANS"]["C_N"] = [1.338, 0.019]

    data["angle_data"]["PEPTIDE"] = dict()
    data["angle_data"]["PEPTIDE"]["CA_C_N"] = [117.2, 2.2]
    data["angle_data"]["PEPTIDE"]["O_C_N"] = [122.7, 1.6]
    data["angle_data"]["PEPTIDE"]["C_N_CA"] = [121.7, 2.5]

    data["angle_data"]["GLY"] = dict()
    data["angle_data"]["GLY"]["CA_C_N"] = [116.2, 2.0]
    data["angle_data"]["GLY"]["O_C_N"] = [123.2, 1.7]
    data["angle_data"]["GLY"]["C_N_CA"] = [122.3, 2.1]

    data["angle_data"]["PRO_TRANS"] = dict()
    data["angle_data"]["PRO_TRANS"]["CA_C_N"] = [117.1, 2.8]
    data["angle_data"]["PRO_TRANS"]["O_C_N"] = [121.1, 1.9]
    data["angle_data"]["PRO_TRANS"]["C_N_CA"] = [119.3, 1.5]
    data["angle_data"]["PRO_TRANS"]["C_N_CD"] = [128.4, 2.1]

    data["angle_data"]["PRO_CIS"] = dict()
    data["angle_data"]["PRO_CIS"]["CA_C_N"] = [117.1, 2.8]
    data["angle_data"]["PRO_CIS"]["O_C_N"] = [121.1, 1.9]
    data["angle_data"]["PRO_CIS"]["C_N_CA"] = [127.0, 2.4]
    data["angle_data"]["PRO_CIS"]["C_N_CD"] = [120.6, 2.2]

    return data
