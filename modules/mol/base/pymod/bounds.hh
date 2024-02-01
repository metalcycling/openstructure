#ifndef OST_MOL_EXPORT_BOUNDS_HH
#define OST_MOL_EXPORT_BOUNDS_HH

#include <ost/log.hh>
  
template <typename H>
geom::Vec3 geom_center(const H& h)
{
  return h.GetBounds().GetCenter();
}

template <typename H>
geom::Vec3 geom_size(const H& h)
{
  WARN_DEPRECATED("GetBoundarySize()/boundary_size is deprecated. Use "
                  "bounds.size instead")
  return h.GetBounds().GetSize();
}

template <typename H>
geom::Vec3 geom_start(const H& h)
{
  WARN_DEPRECATED("GetGeometricStart()/geometric_start is deprecated. Use "
                  "bounds.min instead")
  return h.GetBounds().GetMin();
}

template <typename H>
geom::Vec3 geom_end(const H& h)
{
  WARN_DEPRECATED("GetGeometricEnd()/geometric_end is deprecated. Use "
                  "bounds.max instead")
  return h.GetBounds().GetMax();
}

#endif
