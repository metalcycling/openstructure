#include <boost/python.hpp>

#include <ost/mol/entity_property_mapper.hh>

using namespace ost;
using namespace ost::mol;
using namespace boost::python;

namespace {


Real (EntityPropertyMapper::*get_ah_a)(const AtomHandle&) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_av_a)(const AtomView&) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_rh_a)(const ResidueHandle&) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_rv_a)(const ResidueView&) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_ch_a)(const ChainHandle&) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_cv_a)(const ChainView&) const=&EntityPropertyMapper::Get;


Real (EntityPropertyMapper::*get_ah_b)(const AtomHandle&, Real) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_av_b)(const AtomView&, Real) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_rh_b)(const ResidueHandle&, Real) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_rv_b)(const ResidueView&, Real) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_ch_b)(const ChainHandle&, Real) const=&EntityPropertyMapper::Get;
Real (EntityPropertyMapper::*get_cv_b)(const ChainView&, Real) const=&EntityPropertyMapper::Get;

EntityPropertyMapper create_epm(const String& prop_name, char level)
{
  switch(level) {
    case 'A':
    case 'a':
      return EntityPropertyMapper(prop_name, Prop::ATOM);
    case 'R':
    case 'r':
      return EntityPropertyMapper(prop_name, Prop::RESIDUE);
    case 'C':
    case 'c':
      return EntityPropertyMapper(prop_name, Prop::CHAIN);
    case 'U':
    case 'u':
      return EntityPropertyMapper(prop_name, Prop::UNSPECIFIED);
    default:
      throw Error(String("unknown property level '")+level+"'");
  }
}

}

void export_EntityPropertyMapper()
{
  class_<EntityPropertyMapper>("EntityPropertyMapper", 
                               init<const String&, 
                                    Prop::Level>(arg("level")=Prop::UNSPECIFIED))
    .def("__init__", &create_epm)
    .def("Get", get_ah_a)
    .def("Get", get_av_a)
    .def("Get", get_rh_a)
    .def("Get", get_rv_a)
    .def("Get", get_ch_a)
    .def("Get", get_cv_a)
    .def("Get", get_ah_b)
    .def("Get", get_av_b)
    .def("Get", get_rh_b)
    .def("Get", get_rv_b)
    .def("Get", get_ch_b)
    .def("Get", get_cv_b)
  ;  
}
