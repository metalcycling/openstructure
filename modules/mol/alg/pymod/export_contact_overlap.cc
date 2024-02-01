#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <ost/mol/alg/contact_overlap.hh>
#include <ost/mol/alg/domain_find.hh>

using namespace boost::python;
using namespace ost;
using namespace ost::mol::alg;

SimilarityMatrix (*com_a)(const ost::seq::AlignmentHandle&,
                          Real, int, int)=&ContactOverlapMap;
SimilarityMatrix (*com_b)(const ost::seq::AlignmentHandle&,
                          Real, OverlapMode, Real)=&ContactOverlapMap;
SimilarityMatrix (*com_c)(const DistanceMatrix&, const DistanceMatrix&,
                          Real, OverlapMode, Real)=&ContactOverlapMap;

void export_contact_overlap()
{
  class_<SimilarityMatrix>("SimilarityMatrix", init<int, Real>())
    .def("Get", &SimilarityMatrix::Get,
         return_value_policy<copy_const_reference>())
    .def("Set", &SimilarityMatrix::Set)
    .def("GetSize", &SimilarityMatrix::GetSize)
    .add_property("size", &SimilarityMatrix::GetSize)
    .def("AsImage", &SimilarityMatrix::AsImage)
    .def("Save", &SimilarityMatrix::Save)
    .def("Load", &SimilarityMatrix::Load).staticmethod("Load")
  ;

  class_<DistanceMatrix, DistanceMatrixPtr>("DistanceMatrix", no_init)
    .def("Create", &DistanceMatrix::Create).staticmethod("Create")
    .def("FromAln", &DistanceMatrix::FromAln,
        (arg("aln"), arg("index")=1)).staticmethod("FromAln")
    .def("Get", &DistanceMatrix::Get,
         return_value_policy<copy_const_reference>())
    .def("Set", &DistanceMatrix::Set)

    .def("GetSize", &DistanceMatrix::GetSize)
    .add_property("size", &DistanceMatrix::GetSize)
  ;

  enum_<OverlapMode>("OverlapMode")
    .value("OM_RELATIVE", OM_RELATIVE)
    .value("OM_ABSOLUTE", OM_ABSOLUTE)
    .export_values()
  ;

  def("DomainsFromCOM", &DomainsFromCOM, (arg("idx_a")=-1, arg("idx_b")=-1));

  def("ContactOverlap", &ContactOverlap, (arg("dist_mat_a"), arg("dist_mat_b"),
                                          arg("tolerance")=0.1,
                                          arg("mode")=OM_RELATIVE,
                                          arg("max_dist")=1e6));
  class_<ContactList>("ContactList", init<>())
    .def(vector_indexing_suite<ContactList>())
  ;
  
  def("Contacts", &Contacts);
  def("CalcContactNodeDegree", &CalcContactNodeDegree);
  def("LocalDistanceTest", &LocalDistanceTest);
  def("ContactOverlapMap", com_a);
  def("ContactOverlapMap", com_b, (arg("mode")=OM_RELATIVE,
      arg("max_dist")=1e6));
  def("ContactOverlapMap", com_c, (arg("mode")=OM_RELATIVE,
      arg("max_dist")=1e6));
  class_<Domain>("Domain", no_init)
    .def_readwrite("views", &Domain::views)
    .def_readwrite("views", &Domain::views)
    .def_readwrite("inter_s", &Domain::inter_s)
    .def_readwrite("inter_d", &Domain::inter_d)
    .def_readwrite("intra_s", &Domain::intra_s)
    .def_readwrite("intra_d", &Domain::intra_d)
  ;
  class_<std::vector<Domain> >("DomainList", init<>())
    .def(vector_indexing_suite<std::vector<Domain> >())
  ;

  class_<Domains, DomainsPtr>("Domains", no_init)
    .def("FromCOM", Domains::FromCOM, 
         (arg("sim"), arg("aln"), arg("idx_a")=-1, 
          arg("idx_b")=-1)).staticmethod("FromCOM")
    .def("FromAln", Domains::FromAln,
         (arg("aln"), arg("tolerance")=1.5,
          arg("radius")=1e6, arg("threshold")=0.5, arg("idx_a")=0, arg("idx_b")=1)).staticmethod("FromAln")
    .def("FromDistMats", Domains::FromDistMats,
         (arg("dmat_a"), arg("dmat_b"), arg("tolerance")=1.5, 
          arg("radius")=1e6, arg("threshold")=0.5, 
          arg("num_threads")=1,
          arg("defined_only")=true,
          arg("adj_map")=false)).staticmethod("FromDistMats")
    .add_property("domains", make_function(&Domains::GetDomains,
                  return_value_policy<copy_const_reference>()))
    .add_property("components", make_function(&Domains::GetComponents,
                  return_value_policy<copy_const_reference>()))
    .add_property("adj_map", &Domains::GetAdjMap)
    .add_property("num_free", &Domains::GetNumFreeResidues)
    .add_property("free_d", &Domains::GetFreeD)
    .add_property("free_s", &Domains::GetFreeS)
  ;
}

