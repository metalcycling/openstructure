#ifndef OST_CONOP_COMPOUND_LIB_BASE_HH
#define OST_CONOP_COMPOUND_LIB_BASE_HH

#include "compound.hh"
#include <ost/mol/residue_handle.hh>

namespace ost { namespace conop {

class CompoundLibBase;
typedef boost::shared_ptr<CompoundLibBase> CompoundLibBasePtr;

class DLLEXPORT_OST_CONOP CompoundLibBase {
public:
  virtual ~CompoundLibBase() {}
  virtual CompoundPtr FindCompound(const String& id, 
                                   Compound::Dialect dialect) const = 0;

  bool IsResidueComplete(const ost::mol::ResidueHandle& res, 
                         bool check_hydrogens, 
                         Compound::Dialect) const;
};

}}
#endif

