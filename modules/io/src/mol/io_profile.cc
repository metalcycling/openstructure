#include "io_profile.hh"

namespace ost { namespace io {

bool IOProfileRegistry::alive = false;
IOProfileRegistry& IOProfileRegistry::Instance()
{
  static IOProfileRegistry reg;
  return reg;
}

IOProfileRegistry::IOProfileRegistry():
  profiles_()
{
  IOProfileRegistry::alive = true;
  profiles_["DEFAULT"]=IOProfile();
}

}}
