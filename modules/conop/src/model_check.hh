#ifndef OST_CONOP_MODEL_CHECK_HH
#define OST_CONOP_MODEL_CHECK_HH

#include <ost/mol/entity_handle.hh>
#include <ost/conop/diag.hh>
#include <ost/conop/compound_lib.hh>
namespace ost { namespace conop {

class DLLEXPORT_OST_CONOP Checker {
public:
  Checker(CompoundLibPtr lib, const mol::EntityHandle& ent,
          Diagnostics& diags): lib_(lib), ent_(ent), diags_(diags),
          checked_unk_res_(false), residues_(ent_.GetResidueList())
  {}
  void CheckForUnknownAtoms();
  void CheckForCompleteness(bool require_hydrogens=false);
  void CheckForNonStandard();
  mol::AtomHandleList GetHydrogens();
  mol::AtomHandleList GetZeroOccupancy();
  
private:
  CompoundLibPtr            lib_;
  mol::EntityHandle         ent_;
  Diagnostics&              diags_;
  bool                      checked_unk_res_;
  mol::ResidueHandleList    residues_;
};

}} /* ost::conop */
#endif
