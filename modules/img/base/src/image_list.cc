//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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

/*
  Author:Andreas Schenk
*/
#include <sstream>

#include "image_list.hh"

namespace ost { namespace img {
ImageList::ImageList():std::vector<ImageHandle>()
{
}
ImageList::ImageList(const ImageList& il):std::vector<ImageHandle>(il)
{
}
ImageList::ImageList(ImageList::iterator start,ImageList::iterator end):std::vector<ImageHandle>(start,end)
{
}
    
ImageHandle ImageList::GetGallery(unsigned int columns,bool border)
{ 
    int bsize=0;
    if(border){
      bsize=1;
    }
    if(columns==0){
        columns=static_cast<int>(round(sqrt(static_cast<Real>(size()))));
    }
    unsigned int rows=static_cast<int>(round(static_cast<Real>(size())/static_cast<Real>(columns)));
    Size s=get_max_size_();
    DataType datatype=get_data_type_();
    if(s[2]>1){
      throw("ImageList::GetGallery() not implemented for 3D datasets");
    }
    ImageHandle result=CreateImage(Extent(Size(columns*(bsize+s[0])+bsize,rows*(bsize+s[1])+bsize)),datatype,SPATIAL);
    for(unsigned int i=0;i<size();++i){
        unsigned int icol=i%columns;
        unsigned int irow=i/columns; //integer division intended
        ImageHandle tmp=CreateImage(operator[](i).GetExtent(),datatype,SPATIAL);
        tmp.Paste(operator[](i));
        tmp.SetSpatialOrigin(Point(icol*(bsize+s[0])+bsize,irow*(bsize+s[1])+bsize));
        result.Paste(tmp);
    }
    return result;
}
ImageHandle ImageList::GetImageStack()
{
    Size s=get_max_size_();
    DataType datatype=get_data_type_();
    if(s[2]>1){
      throw("not implemented in 3D");
    }
    ImageHandle result=CreateImage(Extent(Size(s[0],s[1],size())),datatype,SPATIAL);
    for(unsigned int i=0;i<size();++i){
        ImageHandle tmp=CreateImage(operator[](i).GetExtent(),datatype,SPATIAL);
        tmp.Paste(operator[](i));
        tmp.SetSpatialOrigin(Point(0,0,i));
        result.Paste(tmp);
    }
    return result;
}

Size ImageList::get_max_size_()
{
  Size result;
  for(const_iterator it=begin();it!=end();++it){
    Size s=it->GetSize();      
    result[0]=std::max<int>(s[0],result[0]);
    result[1]=std::max<int>(s[1],result[1]);
    result[2]=std::max<int>(s[2],result[2]);
  }
  return result;
}
DataType ImageList::get_data_type_()
{
  DataType datatype=REAL;
  for(const_iterator it=begin();it!=end();++it){
    if(it->GetType()==COMPLEX){
      datatype=COMPLEX;
    }
  }
  return datatype;
}

void ImageList::Apply(NonModAlgorithm& a) const
{
    for(const_iterator it=begin();it!=end();++it){
        (*it).Apply(a);    
    }
}
void ImageList::ApplyIP(NonModAlgorithm& a) const
{
    for(const_iterator it=begin();it!=end();++it){
        (*it).ApplyIP(a);            
    }
}
void ImageList::ApplyIP(ModIPAlgorithm& a)
{
    for(iterator it=begin();it!=end();++it){
        (*it).ApplyIP(a);    
    }
}
ImageList ImageList::Apply(ModIPAlgorithm& a) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).Apply(a));    
    }
    return result;
}
void ImageList::ApplyIP(const ConstModIPAlgorithm& a)
{
    for(iterator it=begin();it!=end();++it){
        (*it).ApplyIP(a);    
    }
}
ImageList ImageList::Apply(const ConstModIPAlgorithm& a) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).Apply(a));    
    }
    return result;
}
void ImageList::ApplyIP(ModOPAlgorithm& a)
{
    for(iterator it=begin();it!=end();++it){
        (*it).ApplyIP(a);    
    }
}
ImageList ImageList::Apply(ModOPAlgorithm& a) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).Apply(a));    
    }
    return result;
}
void ImageList::ApplyIP(const ConstModOPAlgorithm& a)
{
    for(iterator it=begin();it!=end();++it){
        (*it).ApplyIP(a);    
    }
}
ImageList ImageList::Apply(const ConstModOPAlgorithm& a) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).Apply(a));    
    }
    return result;
}
void ImageList::StateApply(ImageStateNonModVisitorBase& v) const
{
    for(const_iterator it=begin();it!=end();++it){
        (*it).StateApply(v);    
    }
}
void ImageList::StateApplyIP(ImageStateModIPVisitorBase& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it).StateApplyIP(v);    
    }
}
ImageList ImageList::StateApply(ImageStateModIPVisitorBase& v) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).StateApply(v));    
    }
    return result;
}
void ImageList::StateApplyIP(const ImageStateConstModIPVisitorBase& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it).StateApplyIP(v);    
    }
}
ImageList ImageList::StateApply(const ImageStateConstModIPVisitorBase& v) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).StateApply(v));    
    }
    return result;
}
void ImageList::StateApplyIP(ImageStateModOPVisitorBase& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it).StateApplyIP(v);    
    }
}
ImageList ImageList::StateApply(ImageStateModOPVisitorBase& v) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).StateApply(v));    
    }
    return result;
}
void ImageList::StateApplyIP(const ImageStateConstModOPVisitorBase& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it).StateApplyIP(v);    
    }
}
ImageList ImageList::StateApply(const ImageStateConstModOPVisitorBase& v) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).StateApply(v));    
    }
    return result;
}
void ImageList::StateApplyIP(ImageStateMorphVisitorBase& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it).StateApplyIP(v);    
    }
}
ImageList ImageList::StateApply(ImageStateMorphVisitorBase& v) const
{
    ImageList result;
    for(const_iterator it=begin();it!=end();++it){
        result.push_back((*it).StateApply(v));    
    }
    return result;
}
ImageList& ImageList::operator+=(Real v)
{
    for(iterator it=begin();it!=end();++it){
        (*it)+=v;    
    }
    return *this;
}
ImageList& ImageList::operator+=(const Complex& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it)+=v;    
    }
    return *this;
}
ImageList& ImageList::operator-=(Real v)
{
    for(iterator it=begin();it!=end();++it){
        (*it)-=v;    
    }
    return *this;
}
ImageList& ImageList::operator-=(const Complex& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it)-=v;    
    }
    return *this;
}
ImageList& ImageList::operator*=(Real v)
{
    for(iterator it=begin();it!=end();++it){
        (*it)*=v;    
    }
    return *this;
}
ImageList& ImageList::operator*=(const Complex& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it)*=v;    
    }
    return *this;
}
ImageList& ImageList::operator/=(Real v)
{
    for(iterator it=begin();it!=end();++it){
        (*it)/=v;    
    }
    return *this;
}
ImageList& ImageList::operator/=(const Complex& v)
{
    for(iterator it=begin();it!=end();++it){
        (*it)/=v;    
    }
    return *this;
}
ImageList& ImageList::operator+=(const ConstImageHandle& h)
{
    for(iterator it=begin();it!=end();++it){
        (*it)+=h;    
    }
    return *this;
}
ImageList& ImageList::operator-=(const ConstImageHandle& h)
{
    for(iterator it=begin();it!=end();++it){
        (*it)-=h;    
    }
    return *this;
}
ImageList& ImageList::operator*=(const ConstImageHandle& h)
{
    for(iterator it=begin();it!=end();++it){
        (*it)*=h;    
    }
    return *this;
}
ImageList& ImageList::operator/=(const ConstImageHandle& h)
{
    for(iterator it=begin();it!=end();++it){
        (*it)/=h;    
    }
    return *this;
}
}}//ns
