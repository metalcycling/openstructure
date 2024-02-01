#include "minimal_compound_lib.hh"
namespace ost { namespace conop {

#include "standard_compounds.cc"

CompoundMap MinimalCompoundLib::compounds_ = MinimalCompoundLib::InitCompounds();

CompoundMap MinimalCompoundLib::InitCompounds() {
  CompoundMap c;
  c["ALA"] = MakeALACompound();
  c["CYS"] = MakeCYSCompound();
  c["ASP"] = MakeASPCompound();
  c["GLU"] = MakeGLUCompound();
  c["PHE"] = MakePHECompound();
  c["GLY"] = MakeGLYCompound();
  c["HIS"] = MakeHISCompound();
  c["ILE"] = MakeILECompound();
  c["LYS"] = MakeLYSCompound();
  c["LEU"] = MakeLEUCompound();
  c["MET"] = MakeMETCompound();
  c["ASN"] = MakeASNCompound();
  c["PRO"] = MakePROCompound();
  c["GLN"] = MakeGLNCompound();
  c["ARG"] = MakeARGCompound();
  c["SER"] = MakeSERCompound();
  c["THR"] = MakeTHRCompound();
  c["VAL"] = MakeVALCompound();
  c["TRP"] = MakeTRPCompound();
  c["TYR"] = MakeTYRCompound();
  c["G"] = MakeGCompound();
  c["T"] = MakeTCompound();
  c["A"] = MakeACompound();
  c["C"] = MakeCCompound();
  c["U"] = MakeUCompound();
  c["DG"] = MakeDGCompound();
  c["DT"] = MakeDTCompound();
  c["DU"] = MakeDUCompound();
  c["DC"] = MakeDCCompound();
  c["DA"] = MakeDACompound();
  return c;
}


CompoundPtr MinimalCompoundLib::FindCompound(const String& id, 
                                             Compound::Dialect dialect) const
{
   CompoundMap::const_iterator i = MinimalCompoundLib::compounds_.find(id); 
   if (i != MinimalCompoundLib::compounds_.end()) {
     return i->second;
   }
   return CompoundPtr();
}

}}
