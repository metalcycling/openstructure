#ifndef OST_EXPORT_HELPER_VECTOR_HH
#define OST_EXPORT_HELPER_VECTOR_HH

#include <boost/python/def_visitor.hpp>
#include <boost/python/suite/indexing/container_utils.hpp>
/*
  Author: Marco Biasini
 */

namespace geom {

namespace bp=boost::python;
template <typename Container>
class VectorAdditions : 
  public bp::def_visitor<VectorAdditions<Container> > {
public:
 typedef typename Container::value_type value_type;
 typedef Container container_type;
 template <class Class>
 void visit(Class& cl) const
 {
   cl
     .def("__str__", &to_string)
     .def("__init__", make_constructor(&from_iter))
   ;
 }
private:
  static boost::shared_ptr<Container> from_iter(const bp::object& obj)
  {
    boost::shared_ptr<Container> c(new Container);
    bp::container_utils::extend_container(*c.get(), obj);
    return c;
  }
  
  static std::string to_string(Container& cl)
  {
    std::stringstream ss;
    ss << "[";
    bool first=true;
    for (typename Container::const_iterator
         i=cl.begin(), e=cl.end(); i!=e; ++i) {
      if (first) {
        first=false;
      } else {
        ss << ", ";
      }
      ss << (*i);
    }
    ss << "]";
    return ss.str();
  }
};

}

#endif
