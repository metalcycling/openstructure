#include "compound_lib_base.hh"
#include <ost/message.hh>

namespace ost { namespace conop {

bool CompoundLibBase::IsResidueComplete(const ost::mol::ResidueHandle& res, bool check_hydrogens, Compound::Dialect dialect) const{

  CompoundPtr compound = this->FindCompound(res.GetName(), dialect);
  
  if(!compound){
    std::stringstream ss;
    ss << "Could not find " << res.GetName() << " in compound library"; 
    throw Error(ss.str());
  }

  const AtomSpecList& a_spec = compound->GetAtomSpecs();

  for(AtomSpecList::const_iterator it=a_spec.begin(); it!=a_spec.end();++it){

    if (it->is_leaving) continue;

    if((it->element=="H" || it->element=="D") && (!check_hydrogens)) continue;

    if(res.FindAtom(it->name).IsValid() || res.FindAtom(it->alt_name).IsValid()) continue;

    return false; 

  }

  return true;
}

}}//namespace
