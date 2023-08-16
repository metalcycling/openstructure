
#ifndef OST_CONOP_MINIMAL_COMPOUND_LIB_HH
#define OST_CONOP_MINIMAL_COMPOUND_LIB_HH

#include "compound_lib_base.hh"

namespace ost { namespace conop {

class MinimalCompoundLib;
typedef boost::shared_ptr<MinimalCompoundLib> MinimalCompoundLibPtr;

// a minimal compound lib containing the definitions of the 20 standard 
// amino acids and standard nucleotides
class DLLEXPORT_OST_CONOP MinimalCompoundLib : public CompoundLibBase {
public:
  MinimalCompoundLib():
    CompoundLibBase()
  {}
  virtual CompoundPtr FindCompound(const String& id, 
                                   Compound::Dialect dialect) const;
private:
  static CompoundMap InitCompounds();
  // since this information is never going to change, it is shared 
  // between instances of minimal compound lib.
  static CompoundMap compounds_;
};


}}
#endif
