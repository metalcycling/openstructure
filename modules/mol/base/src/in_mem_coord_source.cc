#include "in_mem_coord_source.hh"

namespace ost { namespace mol {


InMemCoordSource::InMemCoordSource(const AtomHandleList& atoms):
  CoordSource(atoms)
{
  this->SetMutable(true);
}
  
uint InMemCoordSource::GetFrameCount() const
{
  return frames_.size();
}

CoordFramePtr InMemCoordSource::GetFrame(uint frame_id) const
{
  return frame_id>=frames_.size() ? CoordFramePtr() : frames_[frame_id];
}

void InMemCoordSource::AddFrame(const CoordFramePtr& frame)
{
  frames_.push_back(frame);
}

void InMemCoordSource::AddFrame(const std::vector<geom::Vec3>& coords)
{
  CoordFramePtr fp(new CoordFrame(coords));
  frames_.push_back(fp);
}

void InMemCoordSource::AddFrame(const std::vector<geom::Vec3>& coords,const geom::Vec3& cell_size,const geom::Vec3& cell_angles)
{
  CoordFramePtr fp(new CoordFrame(coords,cell_size,cell_angles));
  frames_.push_back(fp);
}
  
void InMemCoordSource::InsertFrame(int pos, const std::vector<geom::Vec3>& coords)
{
  CoordFrameList::iterator it = frames_.begin();
  for(int i=0;i<pos;++i) ++it;
  CoordFramePtr fp(new CoordFrame(coords));
  frames_.insert(it,fp);
}

}}
