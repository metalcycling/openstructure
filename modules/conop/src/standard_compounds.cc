CompoundPtr MakeALACompound() {
  CompoundPtr c( new Compound("ALA"));
  c->SetOneLetterCode('A');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C3 H7 N O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "OXT", "OXT", "O", true, true),
    AtomSpec(6, "H", "H", "H", false, true),
    AtomSpec(7, "H2", "HN2", "H", true, true),
    AtomSpec(8, "HA", "HA", "H", false, true),
    AtomSpec(9, "HB1", "1HB", "H", false, true),
    AtomSpec(10, "HB2", "2HB", "H", false, true),
    AtomSpec(11, "HB3", "3HB", "H", false, true),
    AtomSpec(12, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<13; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 6, 1),
    BondSpec(0, 7, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 8, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 5, 1),
    BondSpec(4, 9, 1),
    BondSpec(4, 10, 1),
    BondSpec(4, 11, 1),
    BondSpec(5, 12, 1)
  };
  for (int i=0; i<12; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeCYSCompound() {
  CompoundPtr c( new Compound("CYS"));
  c->SetOneLetterCode('C');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C3 H7 N O2 S");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "SG", "SG", "S", false, true),
    AtomSpec(6, "OXT", "OXT", "O", true, true),
    AtomSpec(7, "H", "H", "H", false, true),
    AtomSpec(8, "H2", "HN2", "H", true, true),
    AtomSpec(9, "HA", "HA", "H", false, true),
    AtomSpec(10, "HB2", "1HB", "H", false, true),
    AtomSpec(11, "HB3", "2HB", "H", false, true),
    AtomSpec(12, "HG", "HG", "H", false, true),
    AtomSpec(13, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<14; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 7, 1),
    BondSpec(0, 8, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 9, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 6, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 10, 1),
    BondSpec(4, 11, 1),
    BondSpec(5, 12, 1),
    BondSpec(6, 13, 1)
  };
  for (int i=0; i<13; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeASPCompound() {
  CompoundPtr c( new Compound("ASP"));
  c->SetOneLetterCode('D');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C4 H7 N O4");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "OD1", "OD1", "O", false, true),
    AtomSpec(7, "OD2", "OD2", "O", false, true),
    AtomSpec(8, "OXT", "OXT", "O", true, true),
    AtomSpec(9, "H", "H", "H", false, true),
    AtomSpec(10, "H2", "HN2", "H", true, true),
    AtomSpec(11, "HA", "HA", "H", false, true),
    AtomSpec(12, "HB2", "HB1", "H", false, true),
    AtomSpec(13, "HB3", "HB2", "H", false, true),
    AtomSpec(14, "HD2", "HD2", "H", false, true),
    AtomSpec(15, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<16; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 9, 1),
    BondSpec(0, 10, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 11, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 8, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 12, 1),
    BondSpec(4, 13, 1),
    BondSpec(5, 6, 2),
    BondSpec(5, 7, 1),
    BondSpec(7, 14, 1),
    BondSpec(8, 15, 1)
  };
  for (int i=0; i<15; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeGLUCompound() {
  CompoundPtr c( new Compound("GLU"));
  c->SetOneLetterCode('E');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C5 H9 N O4");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD", "CD", "C", false, true),
    AtomSpec(7, "OE1", "OE1", "O", false, true),
    AtomSpec(8, "OE2", "OE2", "O", false, true),
    AtomSpec(9, "OXT", "OXT", "O", true, true),
    AtomSpec(10, "H", "H", "H", false, true),
    AtomSpec(11, "H2", "HN2", "H", true, true),
    AtomSpec(12, "HA", "HA", "H", false, true),
    AtomSpec(13, "HB2", "HB1", "H", false, true),
    AtomSpec(14, "HB3", "HB2", "H", false, true),
    AtomSpec(15, "HG2", "HG1", "H", false, true),
    AtomSpec(16, "HG3", "HG2", "H", false, true),
    AtomSpec(17, "HE2", "HE2", "H", false, true),
    AtomSpec(18, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<19; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 10, 1),
    BondSpec(0, 11, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 12, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 9, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 13, 1),
    BondSpec(4, 14, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 15, 1),
    BondSpec(5, 16, 1),
    BondSpec(6, 7, 2),
    BondSpec(6, 8, 1),
    BondSpec(8, 17, 1),
    BondSpec(9, 18, 1)
  };
  for (int i=0; i<18; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakePHECompound() {
  CompoundPtr c( new Compound("PHE"));
  c->SetOneLetterCode('F');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C9 H11 N O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD1", "CD1", "C", false, true),
    AtomSpec(7, "CD2", "CD2", "C", false, true),
    AtomSpec(8, "CE1", "CE1", "C", false, true),
    AtomSpec(9, "CE2", "CE2", "C", false, true),
    AtomSpec(10, "CZ", "CZ", "C", false, true),
    AtomSpec(11, "OXT", "OXT", "O", true, true),
    AtomSpec(12, "H", "H", "H", false, true),
    AtomSpec(13, "H2", "HN2", "H", true, true),
    AtomSpec(14, "HA", "HA", "H", false, true),
    AtomSpec(15, "HB2", "1HB", "H", false, true),
    AtomSpec(16, "HB3", "2HB", "H", false, true),
    AtomSpec(17, "HD1", "HD1", "H", false, true),
    AtomSpec(18, "HD2", "HD2", "H", false, true),
    AtomSpec(19, "HE1", "HE1", "H", false, true),
    AtomSpec(20, "HE2", "HE2", "H", false, true),
    AtomSpec(21, "HZ", "HZ", "H", false, true),
    AtomSpec(22, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<23; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 12, 1),
    BondSpec(0, 13, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 14, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 11, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 15, 1),
    BondSpec(4, 16, 1),
    BondSpec(5, 6, 2),
    BondSpec(5, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 17, 1),
    BondSpec(7, 9, 2),
    BondSpec(7, 18, 1),
    BondSpec(8, 10, 2),
    BondSpec(8, 19, 1),
    BondSpec(9, 10, 1),
    BondSpec(9, 20, 1),
    BondSpec(10, 21, 1),
    BondSpec(11, 22, 1)
  };
  for (int i=0; i<23; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeGLYCompound() {
  CompoundPtr c( new Compound("GLY"));
  c->SetOneLetterCode('G');
  c->SetChemClass(mol::ChemClass('P'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C2 H5 N O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "OXT", "OXT", "O", true, true),
    AtomSpec(5, "H", "H", "H", false, true),
    AtomSpec(6, "H2", "HN2", "H", true, true),
    AtomSpec(7, "HA2", "HA1", "H", false, true),
    AtomSpec(8, "HA3", "HA2", "H", false, true),
    AtomSpec(9, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<10; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 5, 1),
    BondSpec(0, 6, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 7, 1),
    BondSpec(1, 8, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 4, 1),
    BondSpec(4, 9, 1)
  };
  for (int i=0; i<9; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeHISCompound() {
  CompoundPtr c( new Compound("HIS"));
  c->SetOneLetterCode('H');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C6 H10 N3 O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "ND1", "ND1", "N", false, true),
    AtomSpec(7, "CD2", "CD2", "C", false, true),
    AtomSpec(8, "CE1", "CE1", "C", false, true),
    AtomSpec(9, "NE2", "NE2", "N", false, true),
    AtomSpec(10, "OXT", "OXT", "O", true, true),
    AtomSpec(11, "H", "H", "H", false, true),
    AtomSpec(12, "H2", "HN2", "H", true, true),
    AtomSpec(13, "HA", "HA", "H", false, true),
    AtomSpec(14, "HB2", "1HB", "H", false, true),
    AtomSpec(15, "HB3", "2HB", "H", false, true),
    AtomSpec(16, "HD1", "HD1", "H", false, true),
    AtomSpec(17, "HD2", "HD2", "H", false, true),
    AtomSpec(18, "HE1", "HE1", "H", false, true),
    AtomSpec(19, "HE2", "HE2", "H", false, true),
    AtomSpec(20, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<21; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 11, 1),
    BondSpec(0, 12, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 13, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 10, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 14, 1),
    BondSpec(4, 15, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 7, 2),
    BondSpec(6, 8, 2),
    BondSpec(6, 16, 1),
    BondSpec(7, 9, 1),
    BondSpec(7, 17, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 18, 1),
    BondSpec(9, 19, 1),
    BondSpec(10, 20, 1)
  };
  for (int i=0; i<21; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeILECompound() {
  CompoundPtr c( new Compound("ILE"));
  c->SetOneLetterCode('I');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C6 H13 N O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG1", "CG1", "C", false, true),
    AtomSpec(6, "CG2", "CG2", "C", false, true),
    AtomSpec(7, "CD1", "CD1", "C", false, true),
    AtomSpec(8, "OXT", "OXT", "O", true, true),
    AtomSpec(9, "H", "H", "H", false, true),
    AtomSpec(10, "H2", "HN2", "H", true, true),
    AtomSpec(11, "HA", "HA", "H", false, true),
    AtomSpec(12, "HB", "HB", "H", false, true),
    AtomSpec(13, "HG12", "1HG1", "H", false, true),
    AtomSpec(14, "HG13", "2HG1", "H", false, true),
    AtomSpec(15, "HG21", "1HG2", "H", false, true),
    AtomSpec(16, "HG22", "2HG2", "H", false, true),
    AtomSpec(17, "HG23", "3HG2", "H", false, true),
    AtomSpec(18, "HD11", "1HD1", "H", false, true),
    AtomSpec(19, "HD12", "2HD1", "H", false, true),
    AtomSpec(20, "HD13", "3HD1", "H", false, true),
    AtomSpec(21, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<22; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 9, 1),
    BondSpec(0, 10, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 11, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 8, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 6, 1),
    BondSpec(4, 12, 1),
    BondSpec(5, 7, 1),
    BondSpec(5, 13, 1),
    BondSpec(5, 14, 1),
    BondSpec(6, 15, 1),
    BondSpec(6, 16, 1),
    BondSpec(6, 17, 1),
    BondSpec(7, 18, 1),
    BondSpec(7, 19, 1),
    BondSpec(7, 20, 1),
    BondSpec(8, 21, 1)
  };
  for (int i=0; i<21; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeLYSCompound() {
  CompoundPtr c( new Compound("LYS"));
  c->SetOneLetterCode('K');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C6 H15 N2 O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD", "CD", "C", false, true),
    AtomSpec(7, "CE", "CE", "C", false, true),
    AtomSpec(8, "NZ", "NZ", "N", false, true),
    AtomSpec(9, "OXT", "OXT", "O", true, true),
    AtomSpec(10, "H", "H", "H", false, true),
    AtomSpec(11, "H2", "HN2", "H", true, true),
    AtomSpec(12, "HA", "HA", "H", false, true),
    AtomSpec(13, "HB2", "1HB", "H", false, true),
    AtomSpec(14, "HB3", "2HB", "H", false, true),
    AtomSpec(15, "HG2", "1HG", "H", false, true),
    AtomSpec(16, "HG3", "2HG", "H", false, true),
    AtomSpec(17, "HD2", "1HD", "H", false, true),
    AtomSpec(18, "HD3", "2HD", "H", false, true),
    AtomSpec(19, "HE2", "1HE", "H", false, true),
    AtomSpec(20, "HE3", "2HE", "H", false, true),
    AtomSpec(21, "HZ1", "1HZ", "H", false, true),
    AtomSpec(22, "HZ2", "2HZ", "H", false, true),
    AtomSpec(23, "HZ3", "3HZ", "H", false, true),
    AtomSpec(24, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<25; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 10, 1),
    BondSpec(0, 11, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 12, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 9, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 13, 1),
    BondSpec(4, 14, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 15, 1),
    BondSpec(5, 16, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 17, 1),
    BondSpec(6, 18, 1),
    BondSpec(7, 8, 1),
    BondSpec(7, 19, 1),
    BondSpec(7, 20, 1),
    BondSpec(8, 21, 1),
    BondSpec(8, 22, 1),
    BondSpec(8, 23, 1),
    BondSpec(9, 24, 1)
  };
  for (int i=0; i<24; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeLEUCompound() {
  CompoundPtr c( new Compound("LEU"));
  c->SetOneLetterCode('L');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C6 H13 N O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD1", "CD1", "C", false, true),
    AtomSpec(7, "CD2", "CD2", "C", false, true),
    AtomSpec(8, "OXT", "OXT", "O", true, true),
    AtomSpec(9, "H", "H", "H", false, true),
    AtomSpec(10, "H2", "HN2", "H", true, true),
    AtomSpec(11, "HA", "HA", "H", false, true),
    AtomSpec(12, "HB2", "1HB", "H", false, true),
    AtomSpec(13, "HB3", "2HB", "H", false, true),
    AtomSpec(14, "HG", "HG", "H", false, true),
    AtomSpec(15, "HD11", "1HD1", "H", false, true),
    AtomSpec(16, "HD12", "2HD1", "H", false, true),
    AtomSpec(17, "HD13", "3HD1", "H", false, true),
    AtomSpec(18, "HD21", "1HD2", "H", false, true),
    AtomSpec(19, "HD22", "2HD2", "H", false, true),
    AtomSpec(20, "HD23", "3HD2", "H", false, true),
    AtomSpec(21, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<22; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 9, 1),
    BondSpec(0, 10, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 11, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 8, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 12, 1),
    BondSpec(4, 13, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 7, 1),
    BondSpec(5, 14, 1),
    BondSpec(6, 15, 1),
    BondSpec(6, 16, 1),
    BondSpec(6, 17, 1),
    BondSpec(7, 18, 1),
    BondSpec(7, 19, 1),
    BondSpec(7, 20, 1),
    BondSpec(8, 21, 1)
  };
  for (int i=0; i<21; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeMETCompound() {
  CompoundPtr c( new Compound("MET"));
  c->SetOneLetterCode('M');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C5 H11 N O2 S");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "SD", "SD", "S", false, true),
    AtomSpec(7, "CE", "CE", "C", false, true),
    AtomSpec(8, "OXT", "OXT", "O", true, true),
    AtomSpec(9, "H", "H", "H", false, true),
    AtomSpec(10, "H2", "HN2", "H", true, true),
    AtomSpec(11, "HA", "HA", "H", false, true),
    AtomSpec(12, "HB2", "1HB", "H", false, true),
    AtomSpec(13, "HB3", "2HB", "H", false, true),
    AtomSpec(14, "HG2", "1HG", "H", false, true),
    AtomSpec(15, "HG3", "2HG", "H", false, true),
    AtomSpec(16, "HE1", "1HE", "H", false, true),
    AtomSpec(17, "HE2", "2HE", "H", false, true),
    AtomSpec(18, "HE3", "3HE", "H", false, true),
    AtomSpec(19, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<20; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 9, 1),
    BondSpec(0, 10, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 11, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 8, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 12, 1),
    BondSpec(4, 13, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 14, 1),
    BondSpec(5, 15, 1),
    BondSpec(6, 7, 1),
    BondSpec(7, 16, 1),
    BondSpec(7, 17, 1),
    BondSpec(7, 18, 1),
    BondSpec(8, 19, 1)
  };
  for (int i=0; i<19; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeASNCompound() {
  CompoundPtr c( new Compound("ASN"));
  c->SetOneLetterCode('N');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C4 H8 N2 O3");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "OD1", "OD1", "O", false, true),
    AtomSpec(7, "ND2", "ND2", "N", true, true),
    AtomSpec(8, "OXT", "OXT", "O", true, true),
    AtomSpec(9, "H", "H", "H", false, true),
    AtomSpec(10, "H2", "HN2", "H", true, true),
    AtomSpec(11, "HA", "HA", "H", false, true),
    AtomSpec(12, "HB2", "HB1", "H", false, true),
    AtomSpec(13, "HB3", "HB2", "H", false, true),
    AtomSpec(14, "HD21", "HD21", "H", true, true),
    AtomSpec(15, "HD22", "HD22", "H", true, true),
    AtomSpec(16, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<17; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 9, 1),
    BondSpec(0, 10, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 11, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 8, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 12, 1),
    BondSpec(4, 13, 1),
    BondSpec(5, 6, 2),
    BondSpec(5, 7, 1),
    BondSpec(7, 14, 1),
    BondSpec(7, 15, 1),
    BondSpec(8, 16, 1)
  };
  for (int i=0; i<16; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakePROCompound() {
  CompoundPtr c( new Compound("PRO"));
  c->SetOneLetterCode('P');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C5 H9 N O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD", "CD", "C", false, true),
    AtomSpec(7, "OXT", "OXT", "O", true, true),
    AtomSpec(8, "H", "HT1", "H", true, true),
    AtomSpec(9, "HA", "HA", "H", false, true),
    AtomSpec(10, "HB2", "1HB", "H", false, true),
    AtomSpec(11, "HB3", "2HB", "H", false, true),
    AtomSpec(12, "HG2", "1HG", "H", false, true),
    AtomSpec(13, "HG3", "2HG", "H", false, true),
    AtomSpec(14, "HD2", "1HD", "H", false, true),
    AtomSpec(15, "HD3", "2HD", "H", false, true),
    AtomSpec(16, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<17; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 6, 1),
    BondSpec(0, 8, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 9, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 7, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 10, 1),
    BondSpec(4, 11, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 12, 1),
    BondSpec(5, 13, 1),
    BondSpec(6, 14, 1),
    BondSpec(6, 15, 1),
    BondSpec(7, 16, 1)
  };
  for (int i=0; i<17; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeGLNCompound() {
  CompoundPtr c( new Compound("GLN"));
  c->SetOneLetterCode('Q');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C5 H10 N2 O3");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD", "CD", "C", false, true),
    AtomSpec(7, "OE1", "OE1", "O", false, true),
    AtomSpec(8, "NE2", "NE2", "N", false, true),
    AtomSpec(9, "OXT", "OXT", "O", true, true),
    AtomSpec(10, "H", "H", "H", false, true),
    AtomSpec(11, "H2", "HN2", "H", true, true),
    AtomSpec(12, "HA", "HA", "H", false, true),
    AtomSpec(13, "HB2", "1HB", "H", false, true),
    AtomSpec(14, "HB3", "2HB", "H", false, true),
    AtomSpec(15, "HG2", "1HG", "H", false, true),
    AtomSpec(16, "HG3", "2HG", "H", false, true),
    AtomSpec(17, "HE21", "1HE2", "H", false, true),
    AtomSpec(18, "HE22", "2HE2", "H", false, true),
    AtomSpec(19, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<20; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 10, 1),
    BondSpec(0, 11, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 12, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 9, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 13, 1),
    BondSpec(4, 14, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 15, 1),
    BondSpec(5, 16, 1),
    BondSpec(6, 7, 2),
    BondSpec(6, 8, 1),
    BondSpec(8, 17, 1),
    BondSpec(8, 18, 1),
    BondSpec(9, 19, 1)
  };
  for (int i=0; i<19; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeARGCompound() {
  CompoundPtr c( new Compound("ARG"));
  c->SetOneLetterCode('R');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C6 H15 N4 O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD", "CD", "C", false, true),
    AtomSpec(7, "NE", "NE", "N", false, true),
    AtomSpec(8, "CZ", "CZ", "C", false, true),
    AtomSpec(9, "NH1", "NH1", "N", false, true),
    AtomSpec(10, "NH2", "NH2", "N", false, true),
    AtomSpec(11, "OXT", "OXT", "O", true, true),
    AtomSpec(12, "H", "H", "H", false, true),
    AtomSpec(13, "H2", "HN2", "H", true, true),
    AtomSpec(14, "HA", "HA", "H", false, true),
    AtomSpec(15, "HB2", "1HB", "H", false, true),
    AtomSpec(16, "HB3", "2HB", "H", false, true),
    AtomSpec(17, "HG2", "1HG", "H", false, true),
    AtomSpec(18, "HG3", "2HG", "H", false, true),
    AtomSpec(19, "HD2", "1HD", "H", false, true),
    AtomSpec(20, "HD3", "2HD", "H", false, true),
    AtomSpec(21, "HE", "HE", "H", false, true),
    AtomSpec(22, "HH11", "1HH1", "H", false, true),
    AtomSpec(23, "HH12", "2HH1", "H", false, true),
    AtomSpec(24, "HH21", "1HH2", "H", false, true),
    AtomSpec(25, "HH22", "2HH2", "H", false, true),
    AtomSpec(26, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<27; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 12, 1),
    BondSpec(0, 13, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 14, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 11, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 15, 1),
    BondSpec(4, 16, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 17, 1),
    BondSpec(5, 18, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 19, 1),
    BondSpec(6, 20, 1),
    BondSpec(7, 8, 1),
    BondSpec(7, 21, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 2),
    BondSpec(9, 22, 1),
    BondSpec(9, 23, 1),
    BondSpec(10, 24, 1),
    BondSpec(10, 25, 1),
    BondSpec(11, 26, 1)
  };
  for (int i=0; i<26; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeSERCompound() {
  CompoundPtr c( new Compound("SER"));
  c->SetOneLetterCode('S');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C3 H7 N O3");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "OG", "OG", "O", false, true),
    AtomSpec(6, "OXT", "OXT", "O", true, true),
    AtomSpec(7, "H", "H", "H", false, true),
    AtomSpec(8, "H2", "HN2", "H", true, true),
    AtomSpec(9, "HA", "HA", "H", false, true),
    AtomSpec(10, "HB2", "1HB", "H", false, true),
    AtomSpec(11, "HB3", "2HB", "H", false, true),
    AtomSpec(12, "HG", "HG", "H", false, true),
    AtomSpec(13, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<14; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 7, 1),
    BondSpec(0, 8, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 9, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 6, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 10, 1),
    BondSpec(4, 11, 1),
    BondSpec(5, 12, 1),
    BondSpec(6, 13, 1)
  };
  for (int i=0; i<13; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeTHRCompound() {
  CompoundPtr c( new Compound("THR"));
  c->SetOneLetterCode('T');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C4 H9 N O3");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "OG1", "OG1", "O", false, true),
    AtomSpec(6, "CG2", "CG2", "C", false, true),
    AtomSpec(7, "OXT", "OXT", "O", true, true),
    AtomSpec(8, "H", "H", "H", false, true),
    AtomSpec(9, "H2", "HN2", "H", true, true),
    AtomSpec(10, "HA", "HA", "H", false, true),
    AtomSpec(11, "HB", "HB", "H", false, true),
    AtomSpec(12, "HG1", "HG1", "H", false, true),
    AtomSpec(13, "HG21", "1HG2", "H", false, true),
    AtomSpec(14, "HG22", "2HG2", "H", false, true),
    AtomSpec(15, "HG23", "3HG2", "H", false, true),
    AtomSpec(16, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<17; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 8, 1),
    BondSpec(0, 9, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 10, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 7, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 6, 1),
    BondSpec(4, 11, 1),
    BondSpec(5, 12, 1),
    BondSpec(6, 13, 1),
    BondSpec(6, 14, 1),
    BondSpec(6, 15, 1),
    BondSpec(7, 16, 1)
  };
  for (int i=0; i<16; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeVALCompound() {
  CompoundPtr c( new Compound("VAL"));
  c->SetOneLetterCode('V');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C5 H11 N O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG1", "CG1", "C", false, true),
    AtomSpec(6, "CG2", "CG2", "C", false, true),
    AtomSpec(7, "OXT", "OXT", "O", true, true),
    AtomSpec(8, "H", "H", "H", false, true),
    AtomSpec(9, "H2", "HN2", "H", true, true),
    AtomSpec(10, "HA", "HA", "H", false, true),
    AtomSpec(11, "HB", "HB", "H", false, true),
    AtomSpec(12, "HG11", "1HG1", "H", false, true),
    AtomSpec(13, "HG12", "2HG1", "H", false, true),
    AtomSpec(14, "HG13", "3HG1", "H", false, true),
    AtomSpec(15, "HG21", "1HG2", "H", false, true),
    AtomSpec(16, "HG22", "2HG2", "H", false, true),
    AtomSpec(17, "HG23", "3HG2", "H", false, true),
    AtomSpec(18, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<19; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 8, 1),
    BondSpec(0, 9, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 10, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 7, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 6, 1),
    BondSpec(4, 11, 1),
    BondSpec(5, 12, 1),
    BondSpec(5, 13, 1),
    BondSpec(5, 14, 1),
    BondSpec(6, 15, 1),
    BondSpec(6, 16, 1),
    BondSpec(6, 17, 1),
    BondSpec(7, 18, 1)
  };
  for (int i=0; i<18; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeTRPCompound() {
  CompoundPtr c( new Compound("TRP"));
  c->SetOneLetterCode('W');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C11 H12 N2 O2");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD1", "CD1", "C", false, true),
    AtomSpec(7, "CD2", "CD2", "C", false, true),
    AtomSpec(8, "NE1", "NE1", "N", false, true),
    AtomSpec(9, "CE2", "CE2", "C", false, true),
    AtomSpec(10, "CE3", "CE3", "C", false, true),
    AtomSpec(11, "CZ2", "CZ2", "C", false, true),
    AtomSpec(12, "CZ3", "CZ3", "C", false, true),
    AtomSpec(13, "CH2", "CH2", "C", false, true),
    AtomSpec(14, "OXT", "OXT", "O", true, true),
    AtomSpec(15, "H", "H", "H", false, true),
    AtomSpec(16, "H2", "HN2", "H", true, true),
    AtomSpec(17, "HA", "HA", "H", false, true),
    AtomSpec(18, "HB2", "1HB", "H", false, true),
    AtomSpec(19, "HB3", "2HB", "H", false, true),
    AtomSpec(20, "HD1", "HD1", "H", false, true),
    AtomSpec(21, "HE1", "HE1", "H", false, true),
    AtomSpec(22, "HE3", "HE3", "H", false, true),
    AtomSpec(23, "HZ2", "HZ2", "H", false, true),
    AtomSpec(24, "HZ3", "HZ3", "H", false, true),
    AtomSpec(25, "HH2", "HH2", "H", false, true),
    AtomSpec(26, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<27; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 15, 1),
    BondSpec(0, 16, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 17, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 14, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 18, 1),
    BondSpec(4, 19, 1),
    BondSpec(5, 6, 2),
    BondSpec(5, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 20, 1),
    BondSpec(7, 9, 2),
    BondSpec(7, 10, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 21, 1),
    BondSpec(9, 11, 1),
    BondSpec(10, 12, 2),
    BondSpec(10, 22, 1),
    BondSpec(11, 13, 2),
    BondSpec(11, 23, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 24, 1),
    BondSpec(13, 25, 1),
    BondSpec(14, 26, 1)
  };
  for (int i=0; i<28; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeTYRCompound() {
  CompoundPtr c( new Compound("TYR"));
  c->SetOneLetterCode('Y');
  c->SetChemClass(mol::ChemClass('L'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C9 H11 N O3");
  AtomSpec atoms[] = {
    AtomSpec(0, "N", "N", "N", false, true),
    AtomSpec(1, "CA", "CA", "C", false, true),
    AtomSpec(2, "C", "C", "C", false, true),
    AtomSpec(3, "O", "O", "O", false, true),
    AtomSpec(4, "CB", "CB", "C", false, true),
    AtomSpec(5, "CG", "CG", "C", false, true),
    AtomSpec(6, "CD1", "CD1", "C", false, true),
    AtomSpec(7, "CD2", "CD2", "C", false, true),
    AtomSpec(8, "CE1", "CE1", "C", false, true),
    AtomSpec(9, "CE2", "CE2", "C", false, true),
    AtomSpec(10, "CZ", "CZ", "C", false, true),
    AtomSpec(11, "OH", "OH", "O", false, true),
    AtomSpec(12, "OXT", "OXT", "O", true, true),
    AtomSpec(13, "H", "H", "H", false, true),
    AtomSpec(14, "H2", "HN2", "H", true, true),
    AtomSpec(15, "HA", "HA", "H", false, true),
    AtomSpec(16, "HB2", "1HB", "H", false, true),
    AtomSpec(17, "HB3", "2HB", "H", false, true),
    AtomSpec(18, "HD1", "HD1", "H", false, true),
    AtomSpec(19, "HD2", "HD2", "H", false, true),
    AtomSpec(20, "HE1", "HE1", "H", false, true),
    AtomSpec(21, "HE2", "HE2", "H", false, true),
    AtomSpec(22, "HH", "HH", "H", false, true),
    AtomSpec(23, "HXT", "HXT", "H", true, true)
  };
  for (int i=0; i<24; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 13, 1),
    BondSpec(0, 14, 1),
    BondSpec(1, 2, 1),
    BondSpec(1, 4, 1),
    BondSpec(1, 15, 1),
    BondSpec(2, 3, 2),
    BondSpec(2, 12, 1),
    BondSpec(4, 5, 1),
    BondSpec(4, 16, 1),
    BondSpec(4, 17, 1),
    BondSpec(5, 6, 2),
    BondSpec(5, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 18, 1),
    BondSpec(7, 9, 2),
    BondSpec(7, 19, 1),
    BondSpec(8, 10, 2),
    BondSpec(8, 20, 1),
    BondSpec(9, 10, 1),
    BondSpec(9, 21, 1),
    BondSpec(10, 11, 1),
    BondSpec(11, 22, 1),
    BondSpec(12, 23, 1)
  };
  for (int i=0; i<24; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeGCompound() {
  CompoundPtr c( new Compound("G"));
  c->SetOneLetterCode('G');
  c->SetChemClass(mol::ChemClass('R'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C10 H14 N5 O8 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "O2'", "O2*", "O", false, true),
    AtomSpec(12, "C1'", "C1*", "C", false, true),
    AtomSpec(13, "N9", "N9", "N", false, true),
    AtomSpec(14, "C8", "C8", "C", false, true),
    AtomSpec(15, "N7", "N7", "N", false, true),
    AtomSpec(16, "C5", "C5", "C", false, true),
    AtomSpec(17, "C6", "C6", "C", false, true),
    AtomSpec(18, "O6", "O6", "O", false, true),
    AtomSpec(19, "N1", "N1", "N", false, true),
    AtomSpec(20, "C2", "C2", "C", false, true),
    AtomSpec(21, "N2", "N2", "N", false, true),
    AtomSpec(22, "N3", "N3", "N", false, true),
    AtomSpec(23, "C4", "C4", "C", false, true),
    AtomSpec(24, "HOP3", "3HOP", "H", false, true),
    AtomSpec(25, "HOP2", "2HOP", "H", false, true),
    AtomSpec(26, "H5'", "1H5*", "H", false, true),
    AtomSpec(27, "H5''", "2H5*", "H", false, true),
    AtomSpec(28, "H4'", "H4*", "H", false, true),
    AtomSpec(29, "H3'", "H3*", "H", false, true),
    AtomSpec(30, "HO3'", "H3T", "H", true, true),
    AtomSpec(31, "H2'", "H2*", "H", false, true),
    AtomSpec(32, "HO2'", "2HO*", "H", false, true),
    AtomSpec(33, "H1'", "H1*", "H", false, true),
    AtomSpec(34, "H8", "H8", "H", false, true),
    AtomSpec(35, "H1", "H1", "H", false, true),
    AtomSpec(36, "H21", "1H2", "H", false, true),
    AtomSpec(37, "H22", "2H2", "H", false, true)
  };
  for (int i=0; i<38; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 24, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 25, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 26, 1),
    BondSpec(5, 27, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 28, 1),
    BondSpec(7, 12, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 29, 1),
    BondSpec(9, 30, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 12, 1),
    BondSpec(10, 31, 1),
    BondSpec(11, 32, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 33, 1),
    BondSpec(13, 14, 1),
    BondSpec(13, 23, 1),
    BondSpec(14, 15, 2),
    BondSpec(14, 34, 1),
    BondSpec(15, 16, 1),
    BondSpec(16, 17, 1),
    BondSpec(16, 23, 2),
    BondSpec(17, 18, 2),
    BondSpec(17, 19, 1),
    BondSpec(19, 20, 1),
    BondSpec(19, 35, 1),
    BondSpec(20, 21, 1),
    BondSpec(20, 22, 2),
    BondSpec(21, 36, 1),
    BondSpec(21, 37, 1),
    BondSpec(22, 23, 1)
  };
  for (int i=0; i<40; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeTCompound() {
  CompoundPtr c( new Compound("T"));
  c->SetOneLetterCode('T');
  c->SetChemClass(mol::ChemClass('R'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C10 H15 N2 O8 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "C1'", "C1*", "C", false, true),
    AtomSpec(12, "N1", "N1", "N", false, true),
    AtomSpec(13, "C2", "C2", "C", false, true),
    AtomSpec(14, "O2", "O2", "O", false, true),
    AtomSpec(15, "N3", "N3", "N", false, true),
    AtomSpec(16, "C4", "C4", "C", false, true),
    AtomSpec(17, "O4", "O4", "O", false, true),
    AtomSpec(18, "C5", "C5", "C", false, true),
    AtomSpec(19, "C7", "C5M", "C", false, true),
    AtomSpec(20, "C6", "C6", "C", false, true),
    AtomSpec(21, "HOP3", "3HOP", "H", false, true),
    AtomSpec(22, "HOP2", "2HOP", "H", false, true),
    AtomSpec(23, "H5'", "1H5*", "H", false, true),
    AtomSpec(24, "H5''", "2H5*", "H", false, true),
    AtomSpec(25, "H4'", "H4*", "H", false, true),
    AtomSpec(26, "H3'", "H3*", "H", false, true),
    AtomSpec(27, "HO3'", "H3T", "H", true, true),
    AtomSpec(28, "H2'", "1H2*", "H", false, true),
    AtomSpec(29, "H2''", "2H2*", "H", false, true),
    AtomSpec(30, "H1'", "H1*", "H", false, true),
    AtomSpec(31, "H3", "H3", "H", false, true),
    AtomSpec(32, "H71", "1H5M", "H", false, true),
    AtomSpec(33, "H72", "2H5M", "H", false, true),
    AtomSpec(34, "H73", "3H5M", "H", false, true),
    AtomSpec(35, "H6", "H6", "H", false, true)
  };
  for (int i=0; i<36; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 21, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 22, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 23, 1),
    BondSpec(5, 24, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 25, 1),
    BondSpec(7, 11, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 26, 1),
    BondSpec(9, 27, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 28, 1),
    BondSpec(10, 29, 1),
    BondSpec(11, 12, 1),
    BondSpec(11, 30, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 20, 1),
    BondSpec(13, 14, 2),
    BondSpec(13, 15, 1),
    BondSpec(15, 16, 1),
    BondSpec(15, 31, 1),
    BondSpec(16, 17, 2),
    BondSpec(16, 18, 1),
    BondSpec(18, 19, 1),
    BondSpec(18, 20, 2),
    BondSpec(19, 32, 1),
    BondSpec(19, 33, 1),
    BondSpec(19, 34, 1),
    BondSpec(20, 35, 1)
  };
  for (int i=0; i<37; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeACompound() {
  CompoundPtr c( new Compound("A"));
  c->SetOneLetterCode('A');
  c->SetChemClass(mol::ChemClass('R'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C10 H14 N5 O7 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "O2'", "O2*", "O", false, true),
    AtomSpec(12, "C1'", "C1*", "C", false, true),
    AtomSpec(13, "N9", "N9", "N", false, true),
    AtomSpec(14, "C8", "C8", "C", false, true),
    AtomSpec(15, "N7", "N7", "N", false, true),
    AtomSpec(16, "C5", "C5", "C", false, true),
    AtomSpec(17, "C6", "C6", "C", false, true),
    AtomSpec(18, "N6", "N6", "N", false, true),
    AtomSpec(19, "N1", "N1", "N", false, true),
    AtomSpec(20, "C2", "C2", "C", false, true),
    AtomSpec(21, "N3", "N3", "N", false, true),
    AtomSpec(22, "C4", "C4", "C", false, true),
    AtomSpec(23, "HOP3", "3HOP", "H", false, true),
    AtomSpec(24, "HOP2", "2HOP", "H", false, true),
    AtomSpec(25, "H5'", "1H5*", "H", false, true),
    AtomSpec(26, "H5''", "2H5*", "H", false, true),
    AtomSpec(27, "H4'", "H4*", "H", false, true),
    AtomSpec(28, "H3'", "H3*", "H", false, true),
    AtomSpec(29, "HO3'", "H3T", "H", true, true),
    AtomSpec(30, "H2'", "H2*", "H", false, true),
    AtomSpec(31, "HO2'", "2HO*", "H", false, true),
    AtomSpec(32, "H1'", "H1*", "H", false, true),
    AtomSpec(33, "H8", "H8", "H", false, true),
    AtomSpec(34, "H61", "1H6", "H", false, true),
    AtomSpec(35, "H62", "2H6", "H", false, true),
    AtomSpec(36, "H2", "H2", "H", false, true)
  };
  for (int i=0; i<37; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 23, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 24, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 25, 1),
    BondSpec(5, 26, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 27, 1),
    BondSpec(7, 12, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 28, 1),
    BondSpec(9, 29, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 12, 1),
    BondSpec(10, 30, 1),
    BondSpec(11, 31, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 32, 1),
    BondSpec(13, 14, 1),
    BondSpec(13, 22, 1),
    BondSpec(14, 15, 2),
    BondSpec(14, 33, 1),
    BondSpec(15, 16, 1),
    BondSpec(16, 17, 1),
    BondSpec(16, 22, 2),
    BondSpec(17, 18, 1),
    BondSpec(17, 19, 2),
    BondSpec(18, 34, 1),
    BondSpec(18, 35, 1),
    BondSpec(19, 20, 1),
    BondSpec(20, 21, 2),
    BondSpec(20, 36, 1),
    BondSpec(21, 22, 1)
  };
  for (int i=0; i<39; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeCCompound() {
  CompoundPtr c( new Compound("C"));
  c->SetOneLetterCode('C');
  c->SetChemClass(mol::ChemClass('R'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C9 H14 N3 O8 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "O2'", "O2*", "O", false, true),
    AtomSpec(12, "C1'", "C1*", "C", false, true),
    AtomSpec(13, "N1", "N1", "N", false, true),
    AtomSpec(14, "C2", "C2", "C", false, true),
    AtomSpec(15, "O2", "O2", "O", false, true),
    AtomSpec(16, "N3", "N3", "N", false, true),
    AtomSpec(17, "C4", "C4", "C", false, true),
    AtomSpec(18, "N4", "N4", "N", false, true),
    AtomSpec(19, "C5", "C5", "C", false, true),
    AtomSpec(20, "C6", "C6", "C", false, true),
    AtomSpec(21, "HOP3", "3HOP", "H", false, true),
    AtomSpec(22, "HOP2", "2HOP", "H", false, true),
    AtomSpec(23, "H5'", "1H5*", "H", false, true),
    AtomSpec(24, "H5''", "2H5*", "H", false, true),
    AtomSpec(25, "H4'", "H4*", "H", false, true),
    AtomSpec(26, "H3'", "H3*", "H", false, true),
    AtomSpec(27, "HO3'", "H3T", "H", true, true),
    AtomSpec(28, "H2'", "H2*", "H", false, true),
    AtomSpec(29, "HO2'", "2HO*", "H", false, true),
    AtomSpec(30, "H1'", "H1*", "H", false, true),
    AtomSpec(31, "H41", "1H4", "H", false, true),
    AtomSpec(32, "H42", "2H4", "H", false, true),
    AtomSpec(33, "H5", "H5", "H", false, true),
    AtomSpec(34, "H6", "H6", "H", false, true)
  };
  for (int i=0; i<35; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 21, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 22, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 23, 1),
    BondSpec(5, 24, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 25, 1),
    BondSpec(7, 12, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 26, 1),
    BondSpec(9, 27, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 12, 1),
    BondSpec(10, 28, 1),
    BondSpec(11, 29, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 30, 1),
    BondSpec(13, 14, 1),
    BondSpec(13, 20, 1),
    BondSpec(14, 15, 2),
    BondSpec(14, 16, 1),
    BondSpec(16, 17, 2),
    BondSpec(17, 18, 1),
    BondSpec(17, 19, 1),
    BondSpec(18, 31, 1),
    BondSpec(18, 32, 1),
    BondSpec(19, 20, 2),
    BondSpec(19, 33, 1),
    BondSpec(20, 34, 1)
  };
  for (int i=0; i<36; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeUCompound() {
  CompoundPtr c( new Compound("U"));
  c->SetOneLetterCode('U');
  c->SetChemClass(mol::ChemClass('R'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C9 H13 N2 O9 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "O2'", "O2*", "O", false, true),
    AtomSpec(12, "C1'", "C1*", "C", false, true),
    AtomSpec(13, "N1", "N1", "N", false, true),
    AtomSpec(14, "C2", "C2", "C", false, true),
    AtomSpec(15, "O2", "O2", "O", false, true),
    AtomSpec(16, "N3", "N3", "N", false, true),
    AtomSpec(17, "C4", "C4", "C", false, true),
    AtomSpec(18, "O4", "O4", "O", false, true),
    AtomSpec(19, "C5", "C5", "C", false, true),
    AtomSpec(20, "C6", "C6", "C", false, true),
    AtomSpec(21, "HOP3", "3HOP", "H", false, true),
    AtomSpec(22, "HOP2", "2HOP", "H", false, true),
    AtomSpec(23, "H5'", "1H5*", "H", false, true),
    AtomSpec(24, "H5''", "2H5*", "H", false, true),
    AtomSpec(25, "H4'", "H4*", "H", false, true),
    AtomSpec(26, "H3'", "H3*", "H", false, true),
    AtomSpec(27, "HO3'", "H3T", "H", true, true),
    AtomSpec(28, "H2'", "H2*", "H", false, true),
    AtomSpec(29, "HO2'", "2HO*", "H", false, true),
    AtomSpec(30, "H1'", "H1*", "H", false, true),
    AtomSpec(31, "H3", "H3", "H", false, true),
    AtomSpec(32, "H5", "H5", "H", false, true),
    AtomSpec(33, "H6", "H6", "H", false, true)
  };
  for (int i=0; i<34; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 21, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 22, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 23, 1),
    BondSpec(5, 24, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 25, 1),
    BondSpec(7, 12, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 26, 1),
    BondSpec(9, 27, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 12, 1),
    BondSpec(10, 28, 1),
    BondSpec(11, 29, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 30, 1),
    BondSpec(13, 14, 1),
    BondSpec(13, 20, 1),
    BondSpec(14, 15, 2),
    BondSpec(14, 16, 1),
    BondSpec(16, 17, 1),
    BondSpec(16, 31, 1),
    BondSpec(17, 18, 2),
    BondSpec(17, 19, 1),
    BondSpec(19, 20, 2),
    BondSpec(19, 32, 1),
    BondSpec(20, 33, 1)
  };
  for (int i=0; i<35; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeDGCompound() {
  CompoundPtr c( new Compound("DG"));
  c->SetOneLetterCode('G');
  c->SetChemClass(mol::ChemClass('S'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C10 H14 N5 O7 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "C1'", "C1*", "C", false, true),
    AtomSpec(12, "N9", "N9", "N", false, true),
    AtomSpec(13, "C8", "C8", "C", false, true),
    AtomSpec(14, "N7", "N7", "N", false, true),
    AtomSpec(15, "C5", "C5", "C", false, true),
    AtomSpec(16, "C6", "C6", "C", false, true),
    AtomSpec(17, "O6", "O6", "O", false, true),
    AtomSpec(18, "N1", "N1", "N", false, true),
    AtomSpec(19, "C2", "C2", "C", false, true),
    AtomSpec(20, "N2", "N2", "N", false, true),
    AtomSpec(21, "N3", "N3", "N", false, true),
    AtomSpec(22, "C4", "C4", "C", false, true),
    AtomSpec(23, "HOP3", "3HOP", "H", false, true),
    AtomSpec(24, "HOP2", "2HOP", "H", false, true),
    AtomSpec(25, "H5'", "1H5*", "H", false, true),
    AtomSpec(26, "H5''", "2H5*", "H", false, true),
    AtomSpec(27, "H4'", "H4*", "H", false, true),
    AtomSpec(28, "H3'", "H3*", "H", false, true),
    AtomSpec(29, "HO3'", "H3T", "H", true, true),
    AtomSpec(30, "H2'", "1H2*", "H", false, true),
    AtomSpec(31, "H2''", "2H2*", "H", false, true),
    AtomSpec(32, "H1'", "H1*", "H", false, true),
    AtomSpec(33, "H8", "H8", "H", false, true),
    AtomSpec(34, "H1", "H1", "H", false, true),
    AtomSpec(35, "H21", "1H2", "H", false, true),
    AtomSpec(36, "H22", "2H2", "H", false, true)
  };
  for (int i=0; i<37; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 23, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 24, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 25, 1),
    BondSpec(5, 26, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 27, 1),
    BondSpec(7, 11, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 28, 1),
    BondSpec(9, 29, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 30, 1),
    BondSpec(10, 31, 1),
    BondSpec(11, 12, 1),
    BondSpec(11, 32, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 22, 1),
    BondSpec(13, 14, 2),
    BondSpec(13, 33, 1),
    BondSpec(14, 15, 1),
    BondSpec(15, 16, 1),
    BondSpec(15, 22, 2),
    BondSpec(16, 17, 2),
    BondSpec(16, 18, 1),
    BondSpec(18, 19, 1),
    BondSpec(18, 34, 1),
    BondSpec(19, 20, 1),
    BondSpec(19, 21, 2),
    BondSpec(20, 35, 1),
    BondSpec(20, 36, 1),
    BondSpec(21, 22, 1)
  };
  for (int i=0; i<39; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeDTCompound() {
  CompoundPtr c( new Compound("DT"));
  c->SetOneLetterCode('T');
  c->SetChemClass(mol::ChemClass('S'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C10 H15 N2 O8 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "C1'", "C1*", "C", false, true),
    AtomSpec(12, "N1", "N1", "N", false, true),
    AtomSpec(13, "C2", "C2", "C", false, true),
    AtomSpec(14, "O2", "O2", "O", false, true),
    AtomSpec(15, "N3", "N3", "N", false, true),
    AtomSpec(16, "C4", "C4", "C", false, true),
    AtomSpec(17, "O4", "O4", "O", false, true),
    AtomSpec(18, "C5", "C5", "C", false, true),
    AtomSpec(19, "C7", "C5M", "C", false, true),
    AtomSpec(20, "C6", "C6", "C", false, true),
    AtomSpec(21, "HOP3", "3HOP", "H", false, true),
    AtomSpec(22, "HOP2", "2HOP", "H", false, true),
    AtomSpec(23, "H5'", "1H5*", "H", false, true),
    AtomSpec(24, "H5''", "2H5*", "H", false, true),
    AtomSpec(25, "H4'", "H4*", "H", false, true),
    AtomSpec(26, "H3'", "H3*", "H", false, true),
    AtomSpec(27, "HO3'", "H3T", "H", true, true),
    AtomSpec(28, "H2'", "1H2*", "H", false, true),
    AtomSpec(29, "H2''", "2H2*", "H", false, true),
    AtomSpec(30, "H1'", "H1*", "H", false, true),
    AtomSpec(31, "H3", "H3", "H", false, true),
    AtomSpec(32, "H71", "1H5M", "H", false, true),
    AtomSpec(33, "H72", "2H5M", "H", false, true),
    AtomSpec(34, "H73", "3H5M", "H", false, true),
    AtomSpec(35, "H6", "H6", "H", false, true)
  };
  for (int i=0; i<36; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 21, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 22, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 23, 1),
    BondSpec(5, 24, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 25, 1),
    BondSpec(7, 11, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 26, 1),
    BondSpec(9, 27, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 28, 1),
    BondSpec(10, 29, 1),
    BondSpec(11, 12, 1),
    BondSpec(11, 30, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 20, 1),
    BondSpec(13, 14, 2),
    BondSpec(13, 15, 1),
    BondSpec(15, 16, 1),
    BondSpec(15, 31, 1),
    BondSpec(16, 17, 2),
    BondSpec(16, 18, 1),
    BondSpec(18, 19, 1),
    BondSpec(18, 20, 2),
    BondSpec(19, 32, 1),
    BondSpec(19, 33, 1),
    BondSpec(19, 34, 1),
    BondSpec(20, 35, 1)
  };
  for (int i=0; i<37; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeDUCompound() {
  CompoundPtr c( new Compound("DU"));
  c->SetOneLetterCode('U');
  c->SetChemClass(mol::ChemClass('S'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C9 H13 N2 O8 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "C1'", "C1*", "C", false, true),
    AtomSpec(12, "N1", "N1", "N", false, true),
    AtomSpec(13, "C2", "C2", "C", false, true),
    AtomSpec(14, "O2", "O2", "O", false, true),
    AtomSpec(15, "N3", "N3", "N", false, true),
    AtomSpec(16, "C4", "C4", "C", false, true),
    AtomSpec(17, "O4", "O4", "O", false, true),
    AtomSpec(18, "C5", "C5", "C", false, true),
    AtomSpec(19, "C6", "C6", "C", false, true),
    AtomSpec(20, "HOP3", "3HOP", "H", false, true),
    AtomSpec(21, "HOP2", "2HOP", "H", false, true),
    AtomSpec(22, "H5'", "1H5*", "H", false, true),
    AtomSpec(23, "H5''", "2H5*", "H", false, true),
    AtomSpec(24, "H4'", "H4*", "H", false, true),
    AtomSpec(25, "H3'", "H3*", "H", false, true),
    AtomSpec(26, "HO3'", "H3T", "H", true, true),
    AtomSpec(27, "H2'", "1H2*", "H", false, true),
    AtomSpec(28, "H2''", "2H2*", "H", false, true),
    AtomSpec(29, "H1'", "H1*", "H", false, true),
    AtomSpec(30, "H3", "H3", "H", false, true),
    AtomSpec(31, "H5", "H5", "H", false, true),
    AtomSpec(32, "H6", "H6", "H", false, true)
  };
  for (int i=0; i<33; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 20, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 21, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 22, 1),
    BondSpec(5, 23, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 24, 1),
    BondSpec(7, 11, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 25, 1),
    BondSpec(9, 26, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 27, 1),
    BondSpec(10, 28, 1),
    BondSpec(11, 12, 1),
    BondSpec(11, 29, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 19, 1),
    BondSpec(13, 14, 2),
    BondSpec(13, 15, 1),
    BondSpec(15, 16, 1),
    BondSpec(15, 30, 1),
    BondSpec(16, 17, 2),
    BondSpec(16, 18, 1),
    BondSpec(18, 19, 2),
    BondSpec(18, 31, 1),
    BondSpec(19, 32, 1)
  };
  for (int i=0; i<34; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeDCCompound() {
  CompoundPtr c( new Compound("DC"));
  c->SetOneLetterCode('C');
  c->SetChemClass(mol::ChemClass('S'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C9 H14 N3 O7 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "C1'", "C1*", "C", false, true),
    AtomSpec(12, "N1", "N1", "N", false, true),
    AtomSpec(13, "C2", "C2", "C", false, true),
    AtomSpec(14, "O2", "O2", "O", false, true),
    AtomSpec(15, "N3", "N3", "N", false, true),
    AtomSpec(16, "C4", "C4", "C", false, true),
    AtomSpec(17, "N4", "N4", "N", false, true),
    AtomSpec(18, "C5", "C5", "C", false, true),
    AtomSpec(19, "C6", "C6", "C", false, true),
    AtomSpec(20, "HOP3", "3HOP", "H", false, true),
    AtomSpec(21, "HOP2", "2HOP", "H", false, true),
    AtomSpec(22, "H5'", "1H5*", "H", false, true),
    AtomSpec(23, "H5''", "2H5*", "H", false, true),
    AtomSpec(24, "H4'", "H4*", "H", false, true),
    AtomSpec(25, "H3'", "H3*", "H", false, true),
    AtomSpec(26, "HO3'", "H3T", "H", true, true),
    AtomSpec(27, "H2'", "1H2*", "H", false, true),
    AtomSpec(28, "H2''", "2H2*", "H", false, true),
    AtomSpec(29, "H1'", "H1*", "H", false, true),
    AtomSpec(30, "H41", "1H4", "H", false, true),
    AtomSpec(31, "H42", "2H4", "H", false, true),
    AtomSpec(32, "H5", "H5", "H", false, true),
    AtomSpec(33, "H6", "H6", "H", false, true)
  };
  for (int i=0; i<34; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 20, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 21, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 22, 1),
    BondSpec(5, 23, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 24, 1),
    BondSpec(7, 11, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 25, 1),
    BondSpec(9, 26, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 27, 1),
    BondSpec(10, 28, 1),
    BondSpec(11, 12, 1),
    BondSpec(11, 29, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 19, 1),
    BondSpec(13, 14, 2),
    BondSpec(13, 15, 1),
    BondSpec(15, 16, 2),
    BondSpec(16, 17, 1),
    BondSpec(16, 18, 1),
    BondSpec(17, 30, 1),
    BondSpec(17, 31, 1),
    BondSpec(18, 19, 2),
    BondSpec(18, 32, 1),
    BondSpec(19, 33, 1)
  };
  for (int i=0; i<35; ++i) { c->AddBond(bonds[i]); }
  return c;
}
CompoundPtr MakeDACompound() {
  CompoundPtr c( new Compound("DA"));
  c->SetOneLetterCode('A');
  c->SetChemClass(mol::ChemClass('S'));
  c->SetChemType(mol::ChemType('U'));
  c->SetFormula("C10 H14 N5 O6 P");
  AtomSpec atoms[] = {
    AtomSpec(0, "OP3", "O3P", "O", true, true),
    AtomSpec(1, "P", "P", "P", false, true),
    AtomSpec(2, "OP1", "O1P", "O", false, true),
    AtomSpec(3, "OP2", "O2P", "O", false, true),
    AtomSpec(4, "O5'", "O5*", "O", false, true),
    AtomSpec(5, "C5'", "C5*", "C", false, true),
    AtomSpec(6, "C4'", "C4*", "C", false, true),
    AtomSpec(7, "O4'", "O4*", "O", false, true),
    AtomSpec(8, "C3'", "C3*", "C", false, true),
    AtomSpec(9, "O3'", "O3*", "O", false, true),
    AtomSpec(10, "C2'", "C2*", "C", false, true),
    AtomSpec(11, "C1'", "C1*", "C", false, true),
    AtomSpec(12, "N9", "N9", "N", false, true),
    AtomSpec(13, "C8", "C8", "C", false, true),
    AtomSpec(14, "N7", "N7", "N", false, true),
    AtomSpec(15, "C5", "C5", "C", false, true),
    AtomSpec(16, "C6", "C6", "C", false, true),
    AtomSpec(17, "N6", "N6", "N", false, true),
    AtomSpec(18, "N1", "N1", "N", false, true),
    AtomSpec(19, "C2", "C2", "C", false, true),
    AtomSpec(20, "N3", "N3", "N", false, true),
    AtomSpec(21, "C4", "C4", "C", false, true),
    AtomSpec(22, "HOP3", "3HOP", "H", false, true),
    AtomSpec(23, "HOP2", "2HOP", "H", false, true),
    AtomSpec(24, "H5'", "1H5*", "H", false, true),
    AtomSpec(25, "H5''", "2H5*", "H", false, true),
    AtomSpec(26, "H4'", "H4*", "H", false, true),
    AtomSpec(27, "H3'", "H3*", "H", false, true),
    AtomSpec(28, "HO3'", "H3T", "H", true, true),
    AtomSpec(29, "H2'", "1H2*", "H", false, true),
    AtomSpec(30, "H2''", "2H2*", "H", false, true),
    AtomSpec(31, "H1'", "H1*", "H", false, true),
    AtomSpec(32, "H8", "H8", "H", false, true),
    AtomSpec(33, "H61", "1H6", "H", false, true),
    AtomSpec(34, "H62", "2H6", "H", false, true),
    AtomSpec(35, "H2", "H2", "H", false, true)
  };
  for (int i=0; i<36; ++i) { c->AddAtom(atoms[i]); }
  BondSpec bonds[] = {
    BondSpec(0, 1, 1),
    BondSpec(0, 22, 1),
    BondSpec(1, 2, 2),
    BondSpec(1, 3, 1),
    BondSpec(1, 4, 1),
    BondSpec(3, 23, 1),
    BondSpec(4, 5, 1),
    BondSpec(5, 6, 1),
    BondSpec(5, 24, 1),
    BondSpec(5, 25, 1),
    BondSpec(6, 7, 1),
    BondSpec(6, 8, 1),
    BondSpec(6, 26, 1),
    BondSpec(7, 11, 1),
    BondSpec(8, 9, 1),
    BondSpec(8, 10, 1),
    BondSpec(8, 27, 1),
    BondSpec(9, 28, 1),
    BondSpec(10, 11, 1),
    BondSpec(10, 29, 1),
    BondSpec(10, 30, 1),
    BondSpec(11, 12, 1),
    BondSpec(11, 31, 1),
    BondSpec(12, 13, 1),
    BondSpec(12, 21, 1),
    BondSpec(13, 14, 2),
    BondSpec(13, 32, 1),
    BondSpec(14, 15, 1),
    BondSpec(15, 16, 1),
    BondSpec(15, 21, 2),
    BondSpec(16, 17, 1),
    BondSpec(16, 18, 2),
    BondSpec(17, 33, 1),
    BondSpec(17, 34, 1),
    BondSpec(18, 19, 1),
    BondSpec(19, 20, 2),
    BondSpec(19, 35, 1),
    BondSpec(20, 21, 1)
  };
  for (int i=0; i<38; ++i) { c->AddBond(bonds[i]); }
  return c;
}
