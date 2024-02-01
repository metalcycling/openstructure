//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
#include <boost/python.hpp>
using namespace boost::python;

#include <ost/mol/mol.hh>
using namespace ost;
using namespace ost::mol;

struct WrappedVisitor : EntityVisitor
{
    WrappedVisitor(PyObject *p)
        : self(p) {}

    WrappedVisitor(PyObject *p, const EntityVisitor& e)
        : EntityVisitor(e), self(p) {}

    virtual bool VisitAtom(const AtomHandle& atom) {
      return call_method<bool, AtomHandle>(self, "VisitAtom", atom);
    }
    bool VisitAtomDefault(AtomHandle atom) {
      return true;
    }

    virtual bool VisitResidue(const ResidueHandle& residue) {
      return call_method<bool, ResidueHandle>(self, "VisitResidue", residue);
    }

    bool VisitResidueDefault(ResidueHandle residue) {
      return true;
    }
    virtual bool VisitChain(const ChainHandle& chain) {
      return call_method<bool, ChainHandle>(self, "VisitChain", chain);
    }
    bool VisitChainDefault(ChainHandle chain) {
      return true;
    }
    virtual bool VisitBond(const BondHandle& bond) {
      return call_method<bool, BondHandle>(self, "VisitBond", bond);
    }
    bool VisitBondDefault(BondHandle chain) {
      return true;
    }
 private:
    PyObject* self;
};

void export_Visitor()
{
  class_<EntityVisitor, WrappedVisitor>("EntityVisitor", init<>())
    .def("VisitAtom", &WrappedVisitor::VisitAtomDefault)
    .def("VisitResidue", &WrappedVisitor::VisitResidueDefault)
    .def("VisitChain", &WrappedVisitor::VisitChainDefault)
    .def("VisitBond", &WrappedVisitor::VisitBondDefault)
  ;
}
