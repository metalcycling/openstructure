#include <ost/mol/alg/distance_test_common.hh>

namespace ost { namespace mol { namespace alg {


String UniqueAtomIdentifier::GetQualifiedAtomName() const
{
  std::stringstream output_string_s;
  output_string_s << chain_ << "." << residue_name_ << residue_ << "." << atom_;
  return (output_string_s.str());
}

bool UniqueAtomIdentifier::operator==(const UniqueAtomIdentifier& rhs) const
{
  return chain_==rhs.chain_ && residue_==rhs.residue_ && atom_==rhs.atom_;
}

bool UniqueAtomIdentifier::operator<(const UniqueAtomIdentifier& rhs) const
{
  int cc=chain_.compare(rhs.chain_);
  if (cc) {
    return cc<0;
  }
  if (residue_<rhs.residue_) {
    return true;
  } else if (residue_>rhs.residue_) {
   return false;
  }
  return atom_.compare(rhs.atom_)<0;
}

// helper function
String SwappedName(const String& name)
{
  if (name=="OE1") return "OE2";
  if (name=="OE2") return "OE1";

  if (name=="OD1") return "OD2";
  if (name=="OD2") return "OD1";

  if (name=="CG1") return "CG2";
  if (name=="CG2") return "CG1";

  if (name=="CE1") return "CE2";
  if (name=="CE2") return "CE1";

  if (name=="CD1") return "CD2";
  if (name=="CD2") return "CD1";

  if (name=="NH1") return "NH2";
  if (name=="NH2") return "NH1";

  return name;
}

// helper function
bool Swappable(const String& rname, const String& aname)
{
  if (rname=="GLU") {
    return (aname=="OE1" || aname=="OE2");
  }
  if (rname=="ASP") {
    return (aname=="OD1" || aname=="OD2");
  }
  if (rname=="VAL") {

    return (aname=="CG1" || aname=="CG2");
  }
  if (rname=="TYR" || rname=="PHE") {
    return (aname=="CD1" || aname=="CD2" || aname=="CE1" || aname=="CE2");
  }
  if (rname=="LEU") {
    return (aname=="CD1" || aname=="CD2");
  }
  if (rname=="ARG") {
    return (aname=="NH1" || aname=="NH2");
  }
  return false;
}

}}}
