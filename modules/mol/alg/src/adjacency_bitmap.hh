#ifndef OST_MOL_ALG_ADJACENCY_BITMAP_HH
#define OST_MOL_ALG_ADJACENCY_BITMAP_HH

#include <cstring>
#include <boost/shared_ptr.hpp>
#include <ost/img/image_handle.hh>
#include <ost/tri_matrix.hh>
#include <ost/seq/sequence_handle.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/stdint.hh>
#include "similarity_matrix.hh"
#include "contact_overlap.hh"
#include "module_config.hh"

namespace ost { namespace mol { namespace alg {


struct OverlapResult {
  OverlapResult(uint16_t nom, uint16_t denom, uint16_t def): 
    nominator(nom), denominator(denom), defined(def) {}
  OverlapResult(): nominator(0), denominator(0), defined(0) {}
  uint16_t nominator;
  uint16_t denominator;
  uint16_t defined;
};

// A highly optimized representation of the adjacency matrix
class AdjacencyBitmap {
public:
  explicit AdjacencyBitmap(size_t n_vertices): 
    num_vertices_(n_vertices),
    storage_size_(StorageSize(n_vertices))
  { 
    data_ = static_cast<uint64_t*>(malloc(this->num_bytes()));
    end_ = data_+storage_size_*num_vertices_;
    memset(data_, 0, this->num_bytes());
  }

  AdjacencyBitmap(const DistanceMatrix& dmat_a, 
                  const DistanceMatrix& dmat_b,
                  Real threshold, Real radius,
                  int start, int end):
    num_vertices_(end-start),
    storage_size_(StorageSize(end-start)),
    data_(static_cast<uint64_t*>(malloc(this->num_bytes()))),
    end_(data_+storage_size_*num_vertices_) {

    for (int i = start; i < end; ++i) {
      this->set(i-start, i-start, false, false);
      for (int j = i+1; j < end; ++j) {
        Real da = dmat_a.Get(i, j);
        Real db = dmat_b.Get(i, j);
        bool defined = da>=0 && db>=0 && (db < radius || da < radius);
        bool agree = false;
        if (defined) {
          agree = std::abs(da-db)<threshold;
        }
        this->set(i-start, j-start, agree, defined);
        this->set(j-start, i-start, agree, defined);
      }
    }
  }

  AdjacencyBitmap(const SimilarityMatrix& sim, Real threshold):
    num_vertices_(sim.GetSize()),
    storage_size_(StorageSize(sim.GetSize())),
    data_(static_cast<uint64_t*>(malloc(this->num_bytes()))),
    end_(data_+storage_size_*num_vertices_) {
    memset(data_, 0, this->num_bytes());
    for (int i = 0; i < sim.GetSize(); ++i) {
      this->set(i, i, false, false);
      for (int j = i+1; j < sim.GetSize(); ++j) {
        this->set(i, j, sim.Get(i, j)>threshold, sim.Get(i, j)>=0.0);
        this->set(j, i, sim.Get(i, j)>threshold, sim.Get(i, j)>=0.0);
      }
    }
  }

  AdjacencyBitmap(const AdjacencyBitmap& rhs):
    num_vertices_(rhs.num_vertices_),
    storage_size_(rhs.storage_size_),
    data_(static_cast<uint64_t*>(malloc(this->num_bytes()))),
    end_(data_+storage_size_*num_vertices_) {
     memcpy(data_, rhs.data_, this->num_bytes());
  }

  ~AdjacencyBitmap() { 
    if (data_)
      free(data_);
  }

  // calculates the overlap between the direct neighbours of vertex i 
  // and j using pure awesomeness and some bit trickery.
  OverlapResult Overlap(size_t vert_i, size_t vert_j) const; 

  void set(size_t i, size_t j, bool val, bool defined=true) {
    size_t byte_index = i*storage_size_ + j/(4*sizeof(uint64_t));
    assert(byte_index < this->num_bytes());
    size_t bit_index = (2*j) % BITS;
    // the following two variables are required to force use of 
    // 64 bit integer for bit operations.
    uint64_t one = 1;
    uint64_t two = 2;
    assert(bit_index < sizeof(uint64_t)*8);

    if (val) {
      data_[byte_index] |=  (one << bit_index);
    } else {
      data_[byte_index] &= ~(one << bit_index);
    }

    if (defined) {
      data_[byte_index] |=  (two << bit_index);
     } else {
      data_[byte_index] &= ~(two << bit_index);
     }
  }

  bool get(size_t i, size_t j) const {
    uint64_t one = 1;
    size_t byte_index = i*storage_size_ + j/(4*sizeof(uint64_t));
    assert(byte_index < storage_size_*num_vertices_);
    size_t bit_index = (2*j) % BITS;
    assert(bit_index < sizeof(uint64_t)*8);
    return data_[byte_index] & (one << bit_index);
  }

  bool defined(size_t i, size_t j) const {
    uint64_t two = 2;
    size_t byte_index = i*storage_size_ + j/(4*sizeof(uint64_t));
    assert(byte_index < storage_size_*num_vertices_);
    size_t bit_index = (2*j) % BITS;
    assert(bit_index < sizeof(uint64_t)*8);
    return data_[byte_index] & (two << bit_index);
  }

  size_t size() const {
    return num_vertices_;
  }

  size_t num_bytes() const {
    return storage_size_*num_vertices_*sizeof(uint64_t);
  }
  void dump() const;
private:

  AdjacencyBitmap();
  AdjacencyBitmap& operator=(const AdjacencyBitmap&);

  const static size_t BITS=sizeof(uint64_t)*8;
  size_t   num_vertices_;
  size_t   storage_size_;
  uint64_t* data_;
  uint64_t* end_;

  static const uint8_t NUMBER_OF_BITS_SET[];
  static const uint8_t NUMBER_OF_EDGES_SET[];
  static size_t StorageSize(size_t verts) { 
    return verts/(4*sizeof(uint64_t)) + (verts % (sizeof(uint64_t)*4) ? 1 : 0);
  }
};


}}}
#endif

