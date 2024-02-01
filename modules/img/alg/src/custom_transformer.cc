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

#include "custom_transformer.hh"

namespace ost { namespace img { namespace tf {

namespace {

Mat4 calc_inverse(const Mat4& in)
{
  Mat4 m(in);
  int indxc[]={0,0,0,0};
  int indxr[]={0,0,0,0};
  int ipiv[]={0,0,0,0};
  int icol=0, irow=0;
  Real dum=0.0, pivinv=0.0;
  Real b[][1]={{0.0},{0.0},{0.0},{0.0}};
  
  for(int i=0;i<4;++i) {
    Real big=0.0;
    for(int j=0;j<4;++j) {
      if(ipiv[j]!=1) {
	for(int k=0;k<4;++k) {
	  if(ipiv[k]==0) {
	    if(std::fabs(m(j,k)) >= big) {
	      big=std::fabs(m(j,k));
	      irow=j;
	      icol=k;
	    }
	  } else if(ipiv[k] >1) {
	    throw CustomTransformerException("error in matrix inversion");
	  }
	}
      }
    }
    ++(ipiv[icol]);
    if(irow!=icol) {
      for(int l=0;l<4;++l) std::swap(m(irow,l),m(icol,l));
      for(int l=0;l<1;++l) std::swap(b[irow][l], b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if(m(icol,icol)==0.0) {
      throw CustomTransformerException("error in matrix inversion");
    }
    pivinv=1.0/m(icol,icol);
    m(icol,icol)=1.0;
    for(int l=0;l<4;++l) m(icol,l)*=pivinv;
    for(int l=0;l<1;++l) b[icol][l]*=pivinv;

    for(int ll=0;ll<4;++ll) {
      if(ll!=icol) {
	dum=m(ll,icol);
	m(ll,icol)=0.0;
	for(int l=0;l<4;++l) m(ll,l) -=m(icol,l)*dum;
	for(int l=0;l<1;++l) b[ll][l] -=b[icol][l]*dum;
      }
    }
  }
  for(int l=3;l>=0;--l) {
    if(indxr[l]!=indxc[l]) {
      for(int k=0;k<4;++k) {
	std::swap(m(k,indxr[l]),m(k,indxc[l]));
      }
    }
  }

  return m;
}

} // anon ns

CustomTransformer::CustomTransformer(const Mat4& tmat):
  tmat_(tmat),imat_(calc_inverse(tmat)) 
{}

CustomTransformer::CustomTransformer(const Mat4& tmat, const Mat4& imat):
  tmat_(tmat),imat_(imat) 
{}

}}}

