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
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>

#include <boost/filesystem/fstream.hpp>

#include "composite3.hh"
#include "vecmat3_op.hh"


namespace geom {

/************************************************************
*  Line3                                                    *
*************************************************************/
Line3::Line3(const Vec3& from, const Vec3& to):ori_(from)
{
  dir_=(to-from)/Length(to-from);
}

Line3::Line3():
  ori_(Vec3(0,0,0)),dir_(Vec3(1,0,0)) 
{}

Vec3 Line3::At(Real r) const   {return ori_+r*dir_;};
Vec3 Line3::GetOrigin() const    {return ori_;};
Vec3 Line3::GetDirection() const {return dir_;};

std::ostream& operator<<(std::ostream& s, const Line3& l)
{
  s << "[" << l.GetOrigin() << " -> " << l.GetDirection() << "]";
  return s;
}



/************************************************************
*  Plane                                                    *
*************************************************************/

Plane::Plane(const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
  Vec3 ip1=p2-p1;
  Vec3 ip2=p3-p1;
  n_=Normalize(Cross(ip1,ip2));
  p_=-Dot(n_,p1);
}

Plane::Plane(const Vec3& p, const Vec3& n)
{
  n_=Normalize(n);
  p_=-Dot(n_,p);
}

Plane::Plane(const Line3& l,const Vec3& p)
{
  Vec3 ip=l.GetOrigin()-p;
  n_=Normalize(Cross(ip,l.GetDirection()));
  p_=-Dot(n_,p);
}

Plane::Plane(Real a, Real b, Real c, Real d)
{
  n_=Vec3(a,b,c);
  p_=d/Length(n_);
  n_=Normalize(n_);
}

Plane::Plane( Real x, Real y, Real z)
{
  n_=Vec3(1.0/x,1.0/y,1.0/z);
  p_=-1.0/Length(n_);
  n_=Normalize(n_);
}

Vec3 Plane::GetNormal() const
{
  return n_;
};

void Plane::SetNormal(const Vec3& n)
{
  n_=Normalize(n);
}

Real Plane::GetP() const
{
  return p_;
}

Vec3 Plane::At(Real x,Real y) const
{
  return n_[2]!=0.0 ? Vec3(x,y,-(n_[0]*x+n_[1]*y+p_)/n_[2]) : Vec3();
}


/************************************************************
*  Sphere                                                    *
*************************************************************/
Sphere::Sphere():
  origin_(Vec3(0.0,0.0,0.0)),
  radius_(1.0)
{
}
Sphere::Sphere(const Vec3& origin, Real r):
  origin_(origin),
  radius_(std::fabs(r))
{
}
Vec3 Sphere::GetOrigin() const
{
  return origin_;
}
Real Sphere::GetRadius() const
{
  return radius_;
}
void Sphere::SetOrigin(const Vec3& v)
{
  origin_=v;
}
void Sphere::SetRadius(Real r)
{
  radius_=r;
}


/************************************************************
*  Rotation3                                                    *
*************************************************************/


Rotation3::Rotation3():
  q_(1.0,0.0,0.0,0.0),
  origin_()
{
}
Rotation3::Rotation3(Real phi, Real theta, Real psi, const Vec3& origin):
  q_(generate_from_eulers(phi,theta,psi)),
  origin_(origin)
{
}
Rotation3::Rotation3(const Vec3& axis, Real angle, const Vec3& origin):
  q_(generate_from_axis_angle(axis,angle)),
  origin_(origin)
{
}
Rotation3::Rotation3(const Line3& line, Real angle):
  q_(generate_from_axis_angle(line.GetDirection(),angle)),
  origin_(line.GetOrigin())
{
}

Rotation3::Rotation3(const Mat3& rot, const Vec3& origin):
  q_(generate_from_matrix(rot)),
  origin_(origin)
{
}
Rotation3::Rotation3(const Quat& q, const Vec3& origin):
  q_(Normalize(q)),
  origin_(origin)
{
}
Vec3 Rotation3::GetOrigin() const
{
  return origin_;
}
Real Rotation3::GetPhi() const
{
  if(GetTheta()<EPSILON){
    return GetRotationAngle()*
        (q_.GetAxis().GetZ()/std::fabs(q_.GetAxis().GetZ()));
  }else{
    return reduce_angle(M_PI-atan2(q_.x*q_.z+
                                   q_.y*q_.w,
                                   q_.y*q_.z-
                                   q_.x*q_.w));
  }
}
Real Rotation3::GetTheta() const
{
  return reduce_angle(acos(-q_.x*q_.x-
                           q_.y*q_.y+
                           q_.z*q_.z+
                           q_.w*q_.w));
}
Real Rotation3::GetPsi() const
{
  if(GetTheta()<EPSILON){
    return 0.0;
  }else{
    return reduce_angle(atan2(q_.x*q_.z-
                              q_.y*q_.w,
                              q_.y*q_.z+
                              q_.x*q_.w));
  }

}
Quat Rotation3::GetQuat() const
{
  return q_;
}
Vec3 Rotation3::GetRotationAxis() const
{
  return q_.GetAxis();
}
Real Rotation3::GetRotationAngle() const
{
  //return 2.0*acos(q_.GetAngle());
  return q_.GetAngle();
}
Mat3 Rotation3::GetRotationMatrix() const
{
  return q_.ToRotationMatrix();
}

void Rotation3::SetOrigin(const Vec3& o)
{
  origin_=o;
}
void Rotation3::SetPsi(Real psi)
{
  q_=generate_from_eulers(GetPhi(),GetTheta(),psi);
}
void Rotation3::SetTheta(Real theta)
{
  q_=generate_from_eulers(GetPhi(),theta,GetPsi());
}
void Rotation3::SetPhi(Real phi)
{
  q_=generate_from_eulers(phi,GetTheta(),GetPsi());
}
void Rotation3::SetQuat(const Quat& q)
{
  q_=Normalize(q);
}
void Rotation3::SetRotationAxis(const Vec3& v)
{
  q_=generate_from_axis_angle(v,GetRotationAngle());
}
void Rotation3::SetRotationAngle(Real angle)
{
  q_=generate_from_axis_angle(GetRotationAxis(),angle);
}
void Rotation3::SetRotationMatrix(const Mat3& rot)
{
  q_=generate_from_matrix(rot);
}
Vec3 Rotation3::Apply(const Vec3& v) const
{
  return origin_+q_.Rotate(v-origin_);
}

bool Rotation3::operator==(const Rotation3& rhs) const
{
  return q_==rhs.q_;
}

Quat Rotation3::generate_from_eulers(Real phi, Real theta, Real psi)
{
  return Quat(cos(0.5*(phi+psi))*cos(0.5*theta),
                    cos(0.5*(phi-psi))*sin(0.5*theta),
                    sin(0.5*(phi-psi))*sin(0.5*theta),
                    sin(0.5*(phi+psi))*cos(0.5*theta));
}
Quat Rotation3::generate_from_axis_angle(const Vec3& axis, Real angle)
{
  //return Quat(cos(angle/2.0),sin(angle/2.0)*Normalize(axis));
  return Quat(angle, axis);
}

/*
  helper functions, adapted from the boost Quat example
  (C) Copyright Hubert Holin 2001.
*/

Vec3 Rotation3::find_invariant_vector(Mat3 rot)
{
  rot(0,0) -= 1.0;
  rot(1,1) -= 1.0;
  rot(2,2) -= 1.0;

  Real minors[9];
  for(int i=2;i>=0;--i){
    for(int j=2;j>=0;--j){
      minors[3*(2-i)+(2-j)]=Minor(rot,i,j);
    }
  }
  Real* where = ::std::max_element(minors, minors+9);

  Real det = *where;

  // assume this never happens
  // if (det <= numeric_limits<Real>::epsilon())

  Vec3 result;

  switch (where-minors) {
  case 0:
    result.SetZ(1.0);
    result.SetX(minors[2]/det);
    result.SetY(-minors[1]/det);
    break;
  case 1:
    result.SetY(1.0);
    result.SetX(-minors[2]/det);
    result.SetZ(-minors[0]/det);
    break;
  case 2:
    result.SetX(1.0);
    result.SetY(-minors[1]/det);
    result.SetZ(minors[0]/det);
    break;
  case 3:
    result.SetZ(1.0);
    result.SetX(minors[5]/det);
    result.SetY(-minors[4]/det);
    break;
  case 4:
    result.SetY(1.0);
    result.SetX(-minors[5]/det);
    result.SetZ(-minors[3]/det);
    break;
  case 5:
    result.SetX(1.0);
    result.SetY(-minors[4]/det);
    result.SetZ(minors[3]/det);
    break;
  case 6:
    result.SetZ(1.0);
    result.SetX(minors[8]/det);
    result.SetY(-minors[7]/det);
    break;
  case 7:
    result.SetY(1.0);
    result.SetX(-minors[8]/det);
    result.SetZ(-minors[6]/det);
    break;
  case 8:
    result.SetX(1.0);
    result.SetY(-minors[7]/det);
    result.SetZ(minors[6]/det);
    break;
  default:
    assert(false); // should never get here
  }

  result = Normalize(result);
  return result;
}


Vec3 Rotation3::find_orthogonal_vector(const Vec3& xyz)
{
  Real vecnormsqr = Length2(xyz);

  // assume this never happens
  //if (vecnormsqr <= numeric_limits<Real>::epsilon())

  Real components[3] = {
    std::fabs(xyz.GetX()),
    std::fabs(xyz.GetY()),
    std::fabs(xyz.GetZ())
  };

  Real*  where = ::std::min_element(components, components+3);

  Vec3 result;

  switch (where-components) {
  case 0:
    if (*where <= std::numeric_limits<Real>::epsilon()) {
      result=Vec3(1.0,0.0,0.0);
    } else {
      Real lambda = -xyz.GetX()/vecnormsqr;
      result = Vec3(1.0 + lambda*xyz.GetX(),lambda*xyz.GetY(),lambda*xyz.GetZ());
    }
    break;
  case 1:
    if (*where <= std::numeric_limits<Real>::epsilon()) {
      result=Vec3(0.0,1.0,0.0);
    } else {
      Real lambda = -xyz.GetY()/vecnormsqr;
      result = Vec3(lambda*xyz.GetX(),1.0+lambda*xyz.GetY(),lambda*xyz.GetZ());
    }
    break;
  case 2:
    if (*where <= std::numeric_limits<Real>::epsilon()) {
      result=Vec3(0.0,0.0,1.0);
    } else {
      Real lambda = -xyz.GetZ()/vecnormsqr;
      result = Vec3(lambda*xyz.GetX(),lambda*xyz.GetY(),1.0+lambda*xyz.GetZ());
    }
    break;
  default:
    assert(false); // should never get here
  }
  result = Normalize(result);
  return result;
}

Vec3 Rotation3::find_vector_for_BOD(const Vec3& xyz, const Vec3& uvw)
{
  Vec3 nrvo(+xyz.GetY()*uvw.GetZ()-xyz.GetZ()*uvw.GetY(),
      -xyz.GetX()*uvw.GetZ()+xyz.GetZ()*uvw.GetX(),
      +xyz.GetX()*uvw.GetY()-xyz.GetY()*uvw.GetX());
  return nrvo;
}

Quat Rotation3::generate_from_matrix(const Mat3& rot)
{
  if ((std::fabs(rot(0,0) - 1.0) <= std::numeric_limits<Real>::epsilon()) &&
      (std::fabs(rot(1,1) - 1.0) <= std::numeric_limits<Real>::epsilon()) &&
      (std::fabs(rot(2,2) - 1.0) <= std::numeric_limits<Real>::epsilon())) {
    return Quat(1.0,0.0,0.0,0.0);
  }

  Real cos_theta = (rot(0,0)+rot(1,1)+rot(2,2)-1.0)*0.5;
  Real stuff = (cos_theta+1.0)*0.5;
  Real cos_theta_sur_2 = sqrt(stuff);
  Real sin_theta_sur_2 = sqrt(1.0-stuff);

  Vec3 xyz = find_invariant_vector(rot);
  Vec3 uvw = find_orthogonal_vector(xyz);
  Vec3 rst = find_vector_for_BOD(xyz,uvw);

  Real ru = rot(0,0)*uvw.GetX()+rot(0,1)*uvw.GetY()+rot(0,2)*uvw.GetZ();
  Real rv = rot(1,0)*uvw.GetX()+rot(1,1)*uvw.GetY()+rot(1,2)*uvw.GetZ();
  Real rw = rot(2,0)*uvw.GetX()+rot(2,1)*uvw.GetY()+rot(2,2)*uvw.GetZ();

  Real angle_sign_determinator = rst.GetX()*ru+rst.GetY()*rv+rst.GetZ()*rw;

  Quat q;

  if (angle_sign_determinator > +std::numeric_limits<Real>::epsilon())
  {
    q = Quat(cos_theta_sur_2, +xyz.GetX()*sin_theta_sur_2,
             +xyz.GetY()*sin_theta_sur_2, +xyz.GetZ()*sin_theta_sur_2);
  } else if (angle_sign_determinator < -std::numeric_limits<Real>::epsilon())
  {
    q = Quat(cos_theta_sur_2, -xyz.GetX()*sin_theta_sur_2,
             -xyz.GetY()*sin_theta_sur_2, -xyz.GetZ()*sin_theta_sur_2);
  } else {
    Real desambiguator = uvw.GetX()*ru+uvw.GetY()*rv+uvw.GetZ()*rw;
    if(desambiguator >= 1.0) {
      q = Quat(0.0, +xyz.GetX(), +xyz.GetY(), +xyz.GetZ());
    } else {
      q = Quat(0, -xyz.GetX(), -xyz.GetY(), -xyz.GetZ());
    }
  }

  return q;
}

Real Rotation3::reduce_angle(Real ang) const
{
  while(ang>=2.0*M_PI){
    ang-=2.0*M_PI;
  }
  while(ang<0){
    ang+=2.0*M_PI;
  }
  return ang;
}

Rotation3List ImportEulerAngles (const boost::filesystem::path& loc)
{
  boost::filesystem::ifstream infile(loc);
  if(!infile) throw std::string("could not open ")+loc.string();
  std::string line;
  Rotation3List rot_list;
  while (std::getline(infile,line)){
    std::istringstream linestream(line);
    Real psi, theta, phi;
    linestream >> psi >> theta >> phi;
    if (linestream.fail()) throw
        std::string("invalid data format or corrupted file");
    Real psi_rad = psi*M_PI/180.0;
    Real theta_rad = theta*M_PI/180.0;
    Real phi_rad = phi*M_PI/180.0;

    Rotation3 rot(psi_rad,theta_rad,phi_rad);
    rot_list.push_back(rot);
  }
  return rot_list;
}

void ExportEulerAngles (const Rotation3List& rot_list,
                        const boost::filesystem::path& loc)
{
  boost::filesystem::ofstream outfile(loc);
  if(!outfile) throw std::string("could not open ")+loc.string();
  for(Rotation3List::const_iterator iter_list = rot_list.begin();
      iter_list != rot_list.end(); ++iter_list)
  {
    Real psi = iter_list->GetPsi();
    Real theta = iter_list->GetTheta();
    Real phi = iter_list->GetPhi();
    outfile << std::left << std::showpoint << std::fixed << std::setprecision(4)
            << std::setw(10) << psi*180.0/M_PI
            << std::setw(10) << theta*180.0/M_PI
            << std::setw(10) << phi*180.0/M_PI << std::endl;
  }
}

Cuboid::Cuboid(const Vec3& center, const CuboidAxis& a, 
               const CuboidAxis& b, const CuboidAxis& c):
  center_(center)
{
  axes_[0]=a;
  axes_[1]=b;
  axes_[2]=c;
}

Cuboid::Cuboid()
{
}

CuboidAxis::CuboidAxis(const Vec3& axis, Real half_extent):
  axis_(Normalize(axis)), half_extent_(half_extent)
{ }
} // ns
