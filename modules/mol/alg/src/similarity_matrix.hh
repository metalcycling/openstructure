#ifndef OST_MOL_ALG_SIMILARITY_MATRIX_HH
#define OST_MOL_ALG_SIMILARITY_MATRIX_HH

#include <limits>
#include <fstream>
#include <ost/tri_matrix.hh>
#include <ost/img/image.hh>
#include "module_config.hh"

namespace ost { namespace mol { namespace alg {

// triangular distance matrix
class DLLEXPORT_OST_MOL_ALG SimilarityMatrix : public ost::TriMatrix<Real> {
public:
  /// \brief initialise similarity matrix with n datapoints.
  ///
  /// The distances are initialised with the default value
  ///
  /// \param n is the number of datapoints
  /// \param def_val is the default similarity
  SimilarityMatrix(int n, Real def_val=std::numeric_limits<Real>::min()):
    ost::TriMatrix<Real>(n, def_val)
  { }

  ost::img::ImageHandle AsImage(bool full=false) const;

  static SimilarityMatrix Load(const String& filename);

  void Save(const String& filename);

};

}}}

#endif

