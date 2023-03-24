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
#ifndef OST_CONOP_HEURISTIC_HH
#define OST_CONOP_HEURISTIC_HH

#include <ost/mol/entity_handle.hh>
#include "minimal_compound_lib.hh"
#include "diag.hh"
#include "processor.hh"

namespace ost { namespace conop {


class HeuristicProcessor;

typedef boost::shared_ptr<HeuristicProcessor> HeuristicProcessorPtr;

class DLLEXPORT_OST_CONOP HeuristicProcessor  : public Processor {
public:
  HeuristicProcessor():
    Processor(),
    lib_()
  {
  }
  virtual ProcessorPtr Copy() const {
    return ProcessorPtr(new HeuristicProcessor(*this));
  }
  HeuristicProcessor(bool bf, bool at, bool cn, bool aa, bool ch,
                     ConopAction zo): 
    Processor(bf, at, cn, aa, ch, zo),
    lib_()
  {}

  virtual String ToString() const;
protected:
  void ProcessUnkResidue(DiagnosticsPtr diags, mol::ResidueHandle res) const;
  virtual void DoProcess(DiagnosticsPtr diags, 
                         mol::EntityHandle ent) const;
private:
  MinimalCompoundLib lib_;
};


}}

#endif

