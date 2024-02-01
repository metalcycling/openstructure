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

#include <ost/log.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/alg/consistency_checks.hh>

namespace ost { namespace mol { namespace alg {
  
bool ResidueNamesMatch(const EntityView& probe, const EntityView& reference, 
                       bool log_as_error)
{
  bool return_value = true;
  ChainViewList ref_chains = reference.GetChainList();
  ChainViewList probe_chains = probe.GetChainList();
  for (unsigned int rci = 0; rci < ref_chains.size(); ++rci) {
     ChainView ref_chain = ref_chains[rci];
     if (rci < probe_chains.size()) {
      ChainView probe_chain = probe_chains[rci];     
      ResidueViewList ref_residues = ref_chain.GetResidueList();       
      for (ResidueViewList::iterator rri = ref_residues.begin(),
           rre = ref_residues.end(); rri != rre; ++rri) {
        ResidueView probe_residue = probe_chain.FindResidue(rri->GetNumber());
        if (probe_residue.IsValid()) {
          if (probe_residue.GetName()!=rri->GetName()) {
            return_value = false;
            if (log_as_error) {
              LOG_ERROR("Name mismatch for residue "
                        << probe_residue.GetNumber()
                        << ": in the reference structure(s) is "
                        << rri->GetName() << ", in the model "
                        << probe_residue.GetName());
            } else {
              LOG_WARNING("Name mismatch for residue "
                          << probe_residue.GetNumber()
                          << ": in the reference structure(s) is "
                          << rri->GetName() << ", in the model "
                          << probe_residue.GetName());
            }
          } 
        }            
      }
    }
  }          
  return return_value;
}
 
}}}



