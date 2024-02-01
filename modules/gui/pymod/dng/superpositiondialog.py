#------------------------------------------------------------------------------
# This file is part of the OpenStructure project <www.openstructure.org>
#
# Copyright (C) 2008-2020 by the OpenStructure authors
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------
#
# Authors: Stefan Bienert 
#
from PyQt5 import QtCore, QtGui, QtWidgets
from ost.mol.alg import Superpose
from ost import mol

class ChainComboBox(QtWidgets.QComboBox):
  def __init__(self, ent, gfx, parent=None):
    # class variables
    self.all_chains = 'All'
    QtWidgets.QComboBox.__init__(self, parent)
    self.entity = ent
    self.addItem(self.all_chains)
    for chain in self.entity.chains:
      self.addItem(chain.name)
    if self.count()>0:
      self.setCurrentIndex(0)
    if gfx:
      self.gfx = gfx
      self.highlighted.connect(self._HighlightChain)
    else:
      self.gfx = None

  def focusOutEvent (self, event):
    if self.gfx:
      self.gfx.selection = None

  def SetItems(self, ent, gfx):
    self.clear()
    self.entity = ent
    self.addItem(self.all_chains)
    for chain in self.entity.chains:
      self.addItem(chain.name)
    if self.count()>0:
      self.setCurrentIndex(0)
    if gfx:
      self.gfx = gfx

  def _HighlightChain(self, chain_idx):
    chain = self.itemText(chain_idx)
    if chain != 'All':
      self.gfx.SetSelection(self.entity.Select('cname="%s"' % str(chain)))
    else:
      self.gfx.SetSelection(self.entity.Select(''))

  def _GetSelectedChain(self):
    if self.currentIndex() == -1:
      return mol.EntityHandle()
    elif self.currentText() == self.all_chains:
      return self.entity
    return self.entity.Select('cname="%s"' % str(self.currentText()))

  def _SetSelectedChain(self, chain):
    if hasattr(chain, 'name'):
      name = chain.name
    else:
      name = str(chain)
    for i in range(self.count()):
      if self.itemText(i) == name:
        self.setCurrentIndex(i)
        break
  selected_chain = property(_GetSelectedChain, _SetSelectedChain)

class SuperpositionDialog(QtWidgets.QDialog):
  """
  Provides a graphical user interface to structurally superpose two entities.
  Uses function :func:`~ost.mol.alg.Superpose`. The RMSD of two superposed
  molecules will be stored in attribute ``rmsd``. An index for the selected
  reference molecule will be stored in attribute ``reference``.

  :param ent_one: The first entity
  :type ent_one: :class:`~ost.mol.EntityView`, :class:`~ost.mol.EntityHandle`
                 or :class:`~ost.gfx.Entity`
  :param ent_two: The second entity
  :type ent_two: :class:`~ost.mol.EntityView`, :class:`~ost.mol.EntityHandle`
                 or :class:`~ost.gfx.Entity`

  **Example Usage:**

  .. code-block:: python

    e1=io.LoadPDB('examples/code_fragments/entity/pdb1ake.ent')
    e2=io.LoadPDB('examples/code_fragments/entity/pdb4ake.ent')

    sd = ost.gui.dng.superpositiondialog.SuperpositionDialog(e1, e2)

    g1=gfx.Entity('G1', e1)
    g2=gfx.Entity('G2', e2)
    scene.Add(g1)
    scene.Add(g2)

    if sd.reference == 0:
      scene.CenterOn(g1)
    else:
      scene.CenterOn(g2)

    if sd.rmsd != None:
      LogScript('RMSD: %.3f'%sd.rmsd)
  """

  def __init__(self, ent_one, ent_two, parent=None):
    # class variables
    self.rmsd_superposed_atoms = None
    self.rmsd = None
    self.fraction_superposed = None
    self.superposition_error = None
    self._mmethod_dict = {'number': 'number',
                          'index': 'index',
                          'local alignment': 'local-aln',
                          'global alignment': 'global-aln'}
    self.gfx_one = None
    self.gfx_two = None
    self.gfx_select_one = None
    self.gfx_select_two = None
    QtWidgets.QDialog.__init__(self, parent)
    self.setWindowTitle('Superpose Structures')
    if not isinstance(ent_one, mol.EntityHandle) and \
       not isinstance(ent_one, mol.EntityView):
      n_one = ent_one.GetName()
      self.gfx_one = ent_one
      self.gfx_select_one = self.gfx_one.GetSelection()
      self.ent_one = ent_one.GetView()
    else:
      if isinstance(ent_one, mol.EntityHandle):
        n_one = ent_one.GetName()
      elif isinstance(ent_one, mol.EntityView):
        n_one = ent_one.GetHandle().GetName()
      self.ent_one = ent_one
    if len(n_one) == 0:
      n_one = '1'
    if not isinstance(ent_two, mol.EntityHandle) and \
       not isinstance(ent_two, mol.EntityView):
      n_two = ent_two.GetName()
      self.gfx_two = ent_two
      self.gfx_select_two = self.gfx_two.GetSelection()
      self.ent_two = ent_two.GetView()
    else:
      if isinstance(ent_two, mol.EntityHandle):
        n_two = ent_two.GetName()
      elif isinstance(ent_two, mol.EntityView):
        n_two = ent_two.GetHandle().GetName()
      self.ent_two = ent_two
    if len(n_two) == 0:
      n_two = '2'
    if n_one == n_two:
      n_one = n_one + ' 1'
      n_two = n_two + ' 2' 
    layout = QtWidgets.QGridLayout(self)
    # select reference molecule
    self.reference = 0;
    self._reference = self._ReferenceSelection(n_one, n_two)
    grow = 0
    layout.addWidget(QtWidgets.QLabel("reference"), grow, 0)
    layout.addWidget(self._reference, grow, 1)
    grow += 1
    # chains
    self._chain_one = ChainComboBox(self.ent_one, self.gfx_one, self)
    self._chain_two = ChainComboBox(self.ent_two, self.gfx_two, self)
    layout.addWidget(QtWidgets.QLabel("reference chain"), grow, 0)
    layout.addWidget(self._chain_one, grow, 1)
    grow += 1
    layout.addWidget(QtWidgets.QLabel("chain"), grow, 0)
    layout.addWidget(self._chain_two, grow, 1)
    grow += 1
    # link chain and reference selection
    self._reference.currentIndexChanged.connect(self._ChangeChainSelection)
    # match methods
    self._methods = self._MatchMethods()
    layout.addWidget(QtWidgets.QLabel('match residues by'), grow, 0)
    grow += 1
    layout.addWidget(self._methods)

    # iterative
    self._iterative=None
    self._it_box, self._it_in, self._dist_in = self._ItBox()
    layout.addWidget(self._it_box, grow, 0)
    # atoms
    self._atoms = self._FetchAtoms(self._methods.size(),
                                   self.ent_one,
                                   self.ent_two)
    self._atmselectbx, self._atmselectgrp = self._AtomSelectionBox()
    layout.addWidget(self._atmselectbx, grow, 1)
    grow += 1
    # buttons
    ok_button = QtWidgets.QPushButton("Superpose")
    ok_button.clicked.connect(self.accept)
    cancel_button = QtWidgets.QPushButton("Cancel")
    hbox_layout = QtWidgets.QHBoxLayout()
    hbox_layout.addStretch(1)
    layout.addLayout(hbox_layout, grow, 0, 1, 2)
    grow += 1
    cancel_button.clicked.connect(self.reject)
    self.accepted.connect(self._Superpose)
    hbox_layout.addWidget(cancel_button, 0)
    hbox_layout.addWidget(ok_button, 0)
    ok_button.setDefault(True)
    self.exec_()
    # restore old selections
    if self.gfx_one:
      self.gfx_one.SetSelection(self.gfx_select_one)
    if self.gfx_two:
      self.gfx_two.SetSelection(self.gfx_select_two)

  def _Superpose(self):
    view_one = self._chain_one.selected_chain
    view_two = self._chain_two.selected_chain
    atoms = self._GetAtomSelection()
    try:
      sp = Superpose(view_two, view_one,
                     self._mmethod_dict[str(self._methods.currentText())],
                     atoms, iterative=self._iterative, 
                     max_iterations=self._it_in.value(), 
                     distance_threshold=self._dist_in.value())
    except Exception as e:
      # mark as failed by setting superposition_error and let caller handle it
      self.superposition_error = str(e)
      return
    self.rmsd = sp.rmsd
    if self._iterative:
      self.rmsd_superposed_atoms = sp.rmsd_superposed_atoms
      self.fraction_superposed = sp.fraction_superposed

  def _toggle_atoms(self, checked):
    if checked:
      self._atoms.setEnabled(True)
    else:
      self._atoms.setEnabled(False)

  def _AtomSelectionBox(self):
    bt1 = QtWidgets.QRadioButton('All')
    bt2 = QtWidgets.QRadioButton('Backbone')
    bt3 = QtWidgets.QRadioButton('CA')
    self.cstmbtntxt = 'Custom'
    custom_rbutton = QtWidgets.QRadioButton(self.cstmbtntxt)
    group = QtWidgets.QButtonGroup()
    group.addButton(bt1)
    group.addButton(bt2)
    group.addButton(bt3)
    group.addButton(custom_rbutton)
    bt1.setChecked(True)
    vbox_layout = QtWidgets.QVBoxLayout()
    vbox_layout.addWidget(bt1)
    vbox_layout.addWidget(bt2)
    vbox_layout.addWidget(bt3)
    vbox_layout.addWidget(custom_rbutton)
    vbox_layout.addWidget(self._atoms)
    custom_rbutton.toggled.connect(self._toggle_atoms)
    box = QtWidgets.QGroupBox("atom selection")
    box.setLayout(vbox_layout)
    return box, group

  def _GetAtomSelection(self):
    checkedbtn = self._atmselectgrp.checkedButton()
    if str(checkedbtn.text()) != self.cstmbtntxt:
      return str(checkedbtn.text())
    slctn_model = self._atoms.selectionModel()
    dt_model = slctn_model.model()
    atms = list()
    for idx in slctn_model.selectedRows():
      slctn = dt_model.data(idx, Qt.DisplayRole).toString()
      atms.append(str(slctn))
    return atms

  def _FetchAtoms(self, dim, ent_a, ent_b):
    # fetch list of atoms: only those which are in both entities are considered
    atm_dict = {}
    for atm in ent_a.GetAtomList():
      atm_dict[atm.name] = 0
    for atm in ent_b.GetAtomList():
      if atm.name in atm_dict:
        atm_dict[atm.name] = 1
    atmlst = list()
    for atm in sorted(atm_dict.keys()):
      if atm_dict[atm]:
        atmlst.append(atm)
    elems = QtCore.QStringListModel(atmlst)
    atoms = QtWidgets.QListView(self)
    dim.setHeight(3*dim.height())
    atoms.setFixedSize(dim)
    atoms.setModel(elems)
    atoms.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
    atoms.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    atoms.setEnabled(False)
    return atoms

  def _ReferenceSelection(self, name_a, name_b):
    cbox = QtWidgets.QComboBox()
    cbox.addItem(name_a)
    cbox.addItem(name_b)
    if cbox.count() > 0:
      cbox.setCurrentIndex(0)
    return cbox

  def _toggle_iterative(self, checked):
    if checked:
      self._it_in.setEnabled(True)
      self._dist_in.setEnabled(True)
      self._iterative=True
    else:
      self._it_in.setEnabled(False)
      self._dist_in.setEnabled(False)
      self._iterative=False

  def _ItBox(self):
    bt1 = QtWidgets.QRadioButton("On")
    iteration_label=QtWidgets.QLabel("Max Iterations: ")
    distance_label=QtWidgets.QLabel("Dist Thresh: ")
    iteration_in=QtWidgets.QSpinBox()
    iteration_in.setRange(1,30)
    iteration_in.setValue(8)
    distance_in=QtWidgets.QDoubleSpinBox()
    distance_in.setRange(1.0,10.0)
    distance_in.setValue(3.0)
    distance_in.setDecimals(1)
    distance_in.setSingleStep(0.5)
    iteration_in.setEnabled(False)
    distance_in.setEnabled(False)
    bt1.setChecked(False)
    self._iterative=False
    vbox_layout = QtWidgets.QVBoxLayout()
    vbox_layout.addWidget(bt1)
    vbox_layout.addWidget(iteration_label)
    vbox_layout.addWidget(iteration_in)
    vbox_layout.addWidget(distance_label)
    vbox_layout.addWidget(distance_in)
    vbox_layout.addSpacing(50)
    bt1.toggled.connect(self._toggle_iterative)
    box = QtWidgets.QGroupBox("Iterative")
    box.setLayout(vbox_layout)
    return box,iteration_in, distance_in

  def _ChangeChainSelection(self, index):
    if index == 0:
      self._chain_one.SetItems(self.ent_one, self.gfx_one)
      self._chain_two.SetItems(self.ent_two, self.gfx_two)
      self.reference = 0;
    elif index == 1:
      self._chain_one.SetItems(self.ent_two, self.gfx_two)
      self._chain_two.SetItems(self.ent_one, self.gfx_one)
      self.reference = 1;

  def _MatchMethods(self):
    methods=QtWidgets.QComboBox(self)
    for method in sorted(self._mmethod_dict):
      methods.addItem(method)
    return methods
