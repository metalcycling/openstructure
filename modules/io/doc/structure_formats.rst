Supported Structure File Formats
================================================================================

.. currentmodule:: ost.io

The following file formats are supported by :func:`LoadEntity`:


CRD - CARD format file used by CHARMM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This trajectory file format is used by the CHARMM program suite (Molecular Modelling).

*Recognized File Extensions*
  .crd

PDB - Brookhaven PDB File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Fine grained control over PDB file import is available via the :func:`LoadPDB`
function. The PDB importer supports loading gzipped PDB files, which are auto-
detected by the .gz file extension.

*Recognized File Extensions*
  .ent, .pdb, .ent.gz, .pdb.gz

*Format Name*
  pdb

mmCIF - macromolecular Crystallographic Information File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Fine grained control over mmCIFile import is available via the :func:`LoadMMCIF`
function. Most notably, this gives you access to the :class:`MMCifInfo` class.
The mmCIF importer supports loading gzipped files, which are auto-detected by
the .gz file extension.

*Recognized File Extensions*
  .cif, .cif.gz

*Format Name*
  mmcif

PQR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A variant of the PDB format that contains data related to atom charges and
radii.

*Recognized File Extensions*
  .pqr

*Format Name*
  pqr
  
SDF - Structured Data File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Chemical table (Ctab) file format (V2000; read-only V3000 experimental support),
aka MDL Molfile.
The SDF format does not support residues, chains or atom names natively.
The SDF importer supports loading gzipped files, which are auto-detected by the
.gz file extension.

The reader assigns 1-based atom indices as atom names.
SDF files containing several molecules are loaded into distinct chains,
named after the molecule name in the MOLfile header with a numerical prefix.
Residues are named after the name in the MOLfile header as well.

Chains are written as separate molecules. If a chain contains more than one
residue, they will be merged into a single molecule.

*Recognized File Extensions*
  .sdf, .sdf.gz
  
