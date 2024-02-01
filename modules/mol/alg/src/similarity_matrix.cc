#include "similarity_matrix.hh"

namespace ost { namespace mol { namespace alg {


ost::img::ImageHandle SimilarityMatrix::AsImage(bool full) const
{
  ost::img::Size img_size(this->GetSize(), this->GetSize());
  ost::img::ImageHandle im=ost::img::CreateImage(img_size);
  for (int i=0; i<this->GetSize(); ++i) {

    for (int j=full ? 0 : i+1; j<this->GetSize(); ++j) {
      im.SetReal(ost::img::Point(i, j), this->Get(i, j));
      im.SetReal(ost::img::Point(j, i), this->Get(i, j));
    }
  }
  return im;
}

SimilarityMatrix SimilarityMatrix::Load(const String& filename)
{
  std::ifstream istream(filename.c_str(), std::ios::binary);
  size_t the_size;
  istream.read(reinterpret_cast<char*>(&the_size), sizeof(size_t));
  SimilarityMatrix sim(the_size);
  istream.read(reinterpret_cast<char*>(&sim.Data().front()),
    sizeof(Real)*sim.Data().size());
  return sim;
}

void SimilarityMatrix::Save(const String& filename)
{
  std::ofstream of(filename.c_str(), std::ios::binary);
  size_t the_size=this->GetSize();
  of.write(reinterpret_cast<const char*>(&the_size), sizeof(size_t));
  of.write(reinterpret_cast<const char*>(&this->Data().front()),
    sizeof(Real)*this->Data().size());
}

}}}

