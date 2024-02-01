//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
#ifndef GEOM_POINT_CLOUD_H
#define GEOM_POINT_CLOUD_H

#include <vector>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>

#include "vec3.hh"

namespace geom {

/*
  the pointcloud and its manager are meant to provide
  name-based lookup for the selection mechanism. They are
  defined here in abstract terms, and are then used in e.g.
  the gfx modules - a graphical object ISA point cloud, and 
  the gfx object manager ISA point cloud manager. This way
  named graphical objects can be used in a selection statement,
  when the query is receiving a pointer to a point cloud manager.
*/

class PointCloud
{
public:
  virtual ~PointCloud();
  virtual std::vector<Vec3> GetPoints() const;
  virtual Vec3 GetCenter() const;
  virtual Real GetRadius() const;
  
protected:

  void AddPoint(const Vec3& v);
 

private:

  std::vector<Vec3> point_list_;
};

typedef boost::shared_ptr<PointCloud> PointCloudPtr;


class PointCloudManager
{
public:
  virtual ~PointCloudManager();

  PointCloudPtr GetPointCloud(const std::string& name);

protected:
  void AddPointCloud(const std::string& name, const PointCloudPtr& p);

private:
  std::map<std::string, PointCloudPtr> pmap_;
};


} // ns

#endif
