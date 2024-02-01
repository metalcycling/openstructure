#ifndef OST_TRI_MATRIX_HH
#define OST_TRI_MATRIX_HH

#include <vector>
#include <cassert>
#include <ost/module_config.hh>


namespace ost {

/// \brief triangular matrix template
template <typename T>
class TEMPLATE_EXPORT TriMatrix {
public:
  TriMatrix(int n, const T& def_val=T()):
    data_((n*(n+1))/2, def_val), n_(n)
  { }
  
  void Set(int i, int j, const T& sim)
  {
    data_[this->GetIndex(i, j)]=sim;
  }
  
  const T& Get(int i, int j) const
  {
    return data_[this->GetIndex(i, j)];
  }
  
  T& operator()(int i, int j)
  {
    return data_[this->GetIndex(i, j)];
  }
  
  const T& operator()(int i, int j) const
  {
    return data_[this->GetIndex(i, j)];
  }  
  
  int GetSize() const 
  {
    return n_;
  }
  std::vector<T>& Data()
  {
    return data_;
  }
private:
  int GetIndex(int i, int j) const {
    assert(i<n_);
    assert(j<n_);
    if (j>i) 
      std::swap(j, i);
    return ((2*n_-j+1)*j)/2+i-j;
  }
  std::vector<T> data_;
  int n_;
};

}

#endif
