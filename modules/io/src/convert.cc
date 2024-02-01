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

#include <ost/stdint.hh>
#include <ost/log.hh>
#include <ost/io/io_exception.hh>

#include "convert.hh"

namespace ost { namespace io {

namespace {

template<int CONVERSIONSOURCE,int CONVERSIONTARGET, class VALUETYPE>
class ConvertHelper{
public:
  static VALUETYPE Convert(const VALUETYPE& value)
  {
    throw IOException("invalid conversion");
  }
  static void ConvertIP(VALUETYPE* value)
  {
    throw IOException("invalid conversion");
  }
};

template<int CONVERSIONSOURCETARGET, class VALUETYPE>
class ConvertHelper<CONVERSIONSOURCETARGET,CONVERSIONSOURCETARGET,VALUETYPE>{
public:
  static VALUETYPE Convert(const VALUETYPE& value)
  {
    return value;
  }
  static void ConvertIP(VALUETYPE* value)
  {
    //no op
  }
};



// specialization for source!=target big/little endian conversion
template<class VALUETYPE>
VALUETYPE endian_convert_helper(const VALUETYPE& value)
{
  VALUETYPE result;
  const unsigned char* vp=reinterpret_cast<const unsigned char *>(&value);
  unsigned char* rp=reinterpret_cast<unsigned char *>(&result);
  
  if(sizeof(VALUETYPE)==1) {
    return value;
  }else if(sizeof(VALUETYPE)==2) {
    rp[0]=vp[1];
    rp[1]=vp[0];
  } else if(sizeof(VALUETYPE)==4) {
    rp[0]=vp[3];
    rp[1]=vp[2];
    rp[2]=vp[1];
    rp[3]=vp[0];
  } else if(sizeof(VALUETYPE)==8) {
    rp[0]=vp[7];
    rp[1]=vp[6];
    rp[2]=vp[5];
    rp[3]=vp[4];
    rp[4]=vp[3];
    rp[5]=vp[2];
    rp[6]=vp[1];
    rp[7]=vp[0];
  }
  return result;
}

template<class VALUETYPE>
void endian_ip_convert_helper(VALUETYPE* value)
{
  unsigned char tmp;
  unsigned char* a=reinterpret_cast<unsigned char *>(value);
  
  if(sizeof(VALUETYPE)==1) {
    return ;
  }else if(sizeof(VALUETYPE)==2) {
    tmp=a[0];
    a[0]=a[1];
    a[1]=tmp;
  } else if(sizeof(VALUETYPE)==4) {
    tmp=a[0];
    a[0]=a[3];
    a[3]=tmp;
    tmp=a[1];
    a[1]=a[2];
    a[2]=tmp;
  } else if(sizeof(VALUETYPE)==8) {
    tmp = a[0];
    a[0] = a[7];
    a[7] = tmp;
    tmp = a[1];
    a[1] = a[6];
    a[6] = tmp;
    tmp = a[2];
    a[2] = a[5];
    a[5] =tmp;
    tmp = a[3];
    a[3] = a[4];
    a[4] = tmp;
  }
}

template<class VALUETYPE>
class ConvertHelper<OST_BIG_ENDIAN,OST_LITTLE_ENDIAN,VALUETYPE>{
public:
  static VALUETYPE Convert(const VALUETYPE& value)
  {
    return endian_convert_helper<VALUETYPE>(value);
  }
  static void ConvertIP(VALUETYPE* value)
  {
    endian_ip_convert_helper<VALUETYPE>(value);
  }
};

template<class VALUETYPE>
class ConvertHelper<OST_LITTLE_ENDIAN,OST_BIG_ENDIAN,VALUETYPE>{
public:
  static VALUETYPE Convert(const VALUETYPE& value)
  {
    return endian_convert_helper<VALUETYPE>(value);
  }
  static void ConvertIP(VALUETYPE* value)
  {
    endian_ip_convert_helper<VALUETYPE>(value);
  }
};


// specialization for source!=target vax conversion integer types

// NOTE: the VAX conversion is based on the libvaxdata code by L. M. Baker
// -> a copy of Baker's code was included until OST 1.11 and removed as unused

// byte swapping for int types if local is big endian
template<class VALUETYPE>
class ConvertHelper<OST_BIG_ENDIAN,OST_VAX_DATA,VALUETYPE>{
public:
  static VALUETYPE Convert(const VALUETYPE& value)
  {
    return endian_convert_helper<VALUETYPE>(value);
  }
  static void ConvertIP(VALUETYPE* value)
  {
    endian_ip_convert_helper<VALUETYPE>(value);
  }
};

template<class VALUETYPE>
class ConvertHelper<OST_VAX_DATA,OST_BIG_ENDIAN,VALUETYPE>{
public:
  static VALUETYPE Convert(const VALUETYPE& value)
  {
    return endian_convert_helper<VALUETYPE>(value);
  }
  static void ConvertIP(VALUETYPE* value)
  {
    endian_ip_convert_helper<VALUETYPE>(value);
  }
};


// no op int types if local is little endian
template<class VALUETYPE>
class ConvertHelper<OST_LITTLE_ENDIAN,OST_VAX_DATA,VALUETYPE>{
public:
  static VALUETYPE Convert(const VALUETYPE& value)
  {
    return value;
  }
  static void ConvertIP(VALUETYPE* value)
  {
    //no op
  }
};
template<class VALUETYPE>
class ConvertHelper<OST_VAX_DATA,OST_LITTLE_ENDIAN,VALUETYPE>{
public:
  static VALUETYPE Convert(const VALUETYPE& value)
  {
    return value;
  }
  static void ConvertIP(VALUETYPE* value)
  {
    //no op
  }
};

// specialization for source!=target vax conversion float types
const unsigned int SIGN_BIT=0x80000000;
const unsigned int EXPONENT_MASK=0x7F800000;
const unsigned int EXPONENT_SIZE =8;
const unsigned int VAX_F_EXPONENT_BIAS=(1 << ( EXPONENT_SIZE - 1 ));
const unsigned int IEEE_S_EXPONENT_BIAS=(( 1 << ( EXPONENT_SIZE - 1 ) ) - 1 );
const unsigned int MANTISSA_MASK =0x007FFFFF;
const unsigned int MANTISSA_SIZE=23;
const unsigned int HIDDEN_BIT=( 1 << MANTISSA_SIZE );
const unsigned int EXPONENT_ADJUSTMENT=( 1 + VAX_F_EXPONENT_BIAS - IEEE_S_EXPONENT_BIAS );
const unsigned int IN_PLACE_EXPONENT_ADJUSTMENT =( EXPONENT_ADJUSTMENT << MANTISSA_SIZE );

void to_vax_converter_float_helper(const unsigned int* in, unsigned int * out)
{
  if( (*in & ~SIGN_BIT) == 0){
    *out=0;  // Set IEEE +-zero [e=m=0] to VAX zero [s=e=m=0] 
    return;
  } 
  unsigned int e= *in & EXPONENT_MASK;
  if ( e == EXPONENT_MASK ) {
    // VAX's have no equivalents for IEEE +-Infinity and +-NaN [e=all-1's] 
    // Fixup to VAX +-extrema [e=all-1's] with zero mantissa [m=0] 
    LOG_ERROR("IEEE Infinity or NaN encountered. Fixup to VAX extrema");
    *out = (*in & SIGN_BIT) | EXPONENT_MASK;
  }else{
    e >>= MANTISSA_SIZE;   // Obtain the biased IEEE exponent 
    unsigned int m = *in & MANTISSA_MASK; // Obtain the IEEE mantissa 
    if(e == 0) {                          // Denormalized [e=0, m<>0] 
      m <<= 1; // Adjust representation from 2**(1-bias) to 2**(e-bias) 
      while ( ( m & HIDDEN_BIT ) == 0 ) {
        m <<= 1;
        e -= 1;  // Adjust exponent 
      }
      m &= MANTISSA_MASK; // Adjust mantissa to hidden-bit form 
    }
    if((e += EXPONENT_ADJUSTMENT) <= 0){
      // Silent underflow 
      *out=0;;                                  
    }else if( e > ( 2 * VAX_F_EXPONENT_BIAS - 1 ) ){
      LOG_ERROR("Overflow. Fixup to VAX extrema");
       *out = ( *in & SIGN_BIT ) | ~SIGN_BIT;
    }else{
      // VAX normalized form [e>0] (both mantissas are 23 bits) 
      *out = ( *in & SIGN_BIT ) | ( e << MANTISSA_SIZE ) | m;
    }
  }
}

void from_vax_converter_float_helper(const unsigned int* in, unsigned int * out)
{
  int e = *in & EXPONENT_MASK ;
  if( e == 0 ){
    // If the biased VAX exponent is zero [e=0] 
    if ( ( *in & SIGN_BIT ) == SIGN_BIT ){
      // If negative [s=1] 
      LOG_ERROR("VAX reserved operand fault; fixup to IEEE zero");
    }
    // Set VAX dirty [m<>0] or true [m=0] zero to IEEE +zero [s=e=m=0] 
    *out = 0;
  }else{ 
    // The biased VAX exponent is non-zero [e<>0] 
    e >>= MANTISSA_SIZE;  // Obtain the biased VAX exponent 
    // The  biased  VAX  exponent  has to be adjusted to account for the 
    // right shift of the IEEE mantissa binary point and the  difference 
    // between  the biases in their "excess n" exponent representations. 
    // If the resulting biased IEEE exponent is less than  or  equal  to 
    // zero, the converted IEEE S_float must use subnormal form.         
    if(( e - EXPONENT_ADJUSTMENT ) > 0){
      //Use IEEE normalized form [e>0] 
      //Both mantissas are 23 bits; adjust the exponent field in place
      *out = *in - IN_PLACE_EXPONENT_ADJUSTMENT;
    }else{
      // Use IEEE subnormal form [e=0, m>0]
      *out = ( *in & SIGN_BIT ) | ((HIDDEN_BIT | ( *in & MANTISSA_MASK )) >> ( 1 - e ) );
     }
  }
} 

template<>
class ConvertHelper<OST_BIG_ENDIAN,OST_VAX_DATA,float>{
public:
  static float Convert(const float& value)
  {
    float result=0.0;
    unsigned int tmp=0;
    to_vax_converter_float_helper(reinterpret_cast<const unsigned int *>(&value),&tmp);
    reinterpret_cast<unsigned char *>(&result)[0]=reinterpret_cast<unsigned char *>(&tmp)[1];
    reinterpret_cast<unsigned char *>(&result)[1]=reinterpret_cast<unsigned char *>(&tmp)[0];
    reinterpret_cast<unsigned char *>(&result)[2]=reinterpret_cast<unsigned char *>(&tmp)[3];
    reinterpret_cast<unsigned char *>(&result)[3]=reinterpret_cast<unsigned char *>(&tmp)[2];
    return result;
  }
  static void ConvertIP(float* value)
  {
    unsigned int tmp=0;
    to_vax_converter_float_helper(reinterpret_cast<const unsigned int *>(value),&tmp);
    reinterpret_cast<unsigned char *>(value)[0]=reinterpret_cast<unsigned char *>(&tmp)[1];
    reinterpret_cast<unsigned char *>(value)[1]=reinterpret_cast<unsigned char *>(&tmp)[0];
    reinterpret_cast<unsigned char *>(value)[2]=reinterpret_cast<unsigned char *>(&tmp)[3];
    reinterpret_cast<unsigned char *>(value)[3]=reinterpret_cast<unsigned char *>(&tmp)[2];
  }
};

template<>
class ConvertHelper<OST_VAX_DATA,OST_BIG_ENDIAN,float>{
public:
  static float Convert(const float& value)
  {
    float result=0.0;
    unsigned int tmp=0;
    reinterpret_cast<unsigned char *>(&tmp)[1]=reinterpret_cast<const unsigned char *>(&value)[0];
    reinterpret_cast<unsigned char *>(&tmp)[0]=reinterpret_cast<const unsigned char *>(&value)[1];
    reinterpret_cast<unsigned char *>(&tmp)[3]=reinterpret_cast<const unsigned char *>(&value)[2];
    reinterpret_cast<unsigned char *>(&tmp)[2]=reinterpret_cast<const unsigned char *>(&value)[3];
    from_vax_converter_float_helper(&tmp,reinterpret_cast<unsigned int *>(&result));
    return result;
  }
  static void ConvertIP(float* value)
  {
    unsigned int tmp=0;
    reinterpret_cast<unsigned char *>(&tmp)[1]=reinterpret_cast<unsigned char *>(value)[0];
    reinterpret_cast<unsigned char *>(&tmp)[0]=reinterpret_cast<unsigned char *>(value)[1];
    reinterpret_cast<unsigned char *>(&tmp)[3]=reinterpret_cast<unsigned char *>(value)[2];
    reinterpret_cast<unsigned char *>(&tmp)[2]=reinterpret_cast<unsigned char *>(value)[3];
    from_vax_converter_float_helper(&tmp,reinterpret_cast<unsigned int *>(value));
  }
};

template<>
class ConvertHelper<OST_LITTLE_ENDIAN,OST_VAX_DATA,float>{
public:
  static float Convert(const float& value)
  {
    float result=0.0;
    unsigned int tmp=0;
    to_vax_converter_float_helper(reinterpret_cast<const unsigned int *>(&value),&tmp);
    reinterpret_cast<unsigned short *>(&result)[0]=reinterpret_cast<unsigned short *>(&tmp)[1];
    reinterpret_cast<unsigned short *>(&result)[1]=reinterpret_cast<unsigned short *>(&tmp)[0];
    return result;
  }
  static void ConvertIP(float* value)
  {
    unsigned int tmp=0;
    to_vax_converter_float_helper(reinterpret_cast<const unsigned int *>(value),&tmp);
    reinterpret_cast<unsigned short *>(value)[0]=reinterpret_cast<unsigned short *>(&tmp)[1];
    reinterpret_cast<unsigned short *>(value)[1]=reinterpret_cast<unsigned short *>(&tmp)[0];
  }
};

template<>
class ConvertHelper<OST_VAX_DATA,OST_LITTLE_ENDIAN,float>{
public:
  static float Convert(const float& value)
  {
    float result=0.0;
    unsigned int tmp=0;
    reinterpret_cast<unsigned short *>(&tmp)[1]=reinterpret_cast<const unsigned short *>(&value)[0];
    reinterpret_cast<unsigned short *>(&tmp)[0]=reinterpret_cast<const unsigned short *>(&value)[1];
    from_vax_converter_float_helper(&tmp,reinterpret_cast<unsigned int *>(&result));
    return result;
  }
  static void ConvertIP(float* value)
  {
    unsigned int tmp=0;
    reinterpret_cast<unsigned short *>(&tmp)[1]=reinterpret_cast<unsigned short *>(value)[0];
    reinterpret_cast<unsigned short *>(&tmp)[0]=reinterpret_cast<unsigned short *>(value)[1];
    from_vax_converter_float_helper(&tmp,reinterpret_cast<unsigned int *>(value));
  }
};

template<>
class ConvertHelper<OST_BIG_ENDIAN,OST_VAX_DATA,double>{
public:
  static double Convert(const double& value)
  {
    throw IOException("not implemented");
  }
  static void ConvertIP(double* value)
  {
    throw IOException("not implemented");
  }
};

template<>
class ConvertHelper<OST_VAX_DATA,OST_BIG_ENDIAN,double>{
public:
  static double Convert(const double& value)
  {
    throw IOException("not implemented");
  }
  static void ConvertIP(double* value)
  {
    throw IOException("not implemented");
  }
};

template<>
class ConvertHelper<OST_LITTLE_ENDIAN,OST_VAX_DATA,double>{
public:
  static double Convert(const double& value)
  {
    throw IOException("not implemented");
  }
  static void ConvertIP(double* value)
  {
    throw IOException("not implemented");
  }
};

template<>
class ConvertHelper<OST_VAX_DATA,OST_LITTLE_ENDIAN,double>{
public:
  static double Convert(const double& value)
  {
    throw IOException("not implemented");
  }
  static void ConvertIP(double* value)
  {
    throw IOException("not implemented");
  }
};



} //anon ns

template<int CONVERSIONTYPE, class VALUETYPE>
VALUETYPE Convert<CONVERSIONTYPE,VALUETYPE>::From(const VALUETYPE& value)
{
  throw IOException("invalid conversion");
}

template<int CONVERSIONTYPE, class VALUETYPE>
VALUETYPE Convert<CONVERSIONTYPE,VALUETYPE>::To(const VALUETYPE& value)
{
  throw IOException("invalid conversion");
}

template<int CONVERSIONTYPE, class VALUETYPE>
void Convert<CONVERSIONTYPE,VALUETYPE>::FromIP(VALUETYPE* value)
{
  throw IOException("invalid conversion");
}

template<int CONVERSIONTYPE, class VALUETYPE>
void Convert<CONVERSIONTYPE,VALUETYPE>::ToIP(VALUETYPE* value)
{
  throw IOException("invalid conversion");
}


template<class VALUETYPE>
class Convert<OST_BIG_ENDIAN,VALUETYPE>{
public:
  static VALUETYPE From(const VALUETYPE& value)
  {
    return ConvertHelper<OST_BIG_ENDIAN,OST_LOCAL_ENDIAN,VALUETYPE>::Convert(value);
  }
  static VALUETYPE To(const VALUETYPE& value)
  {
    return ConvertHelper<OST_LOCAL_ENDIAN,OST_BIG_ENDIAN,VALUETYPE>::Convert(value);
  }
  static void FromIP(VALUETYPE* value)
  {
    ConvertHelper<OST_BIG_ENDIAN,OST_LOCAL_ENDIAN,VALUETYPE>::ConvertIP(value);
  }
  static void ToIP(VALUETYPE* value)
  {
    ConvertHelper<OST_LOCAL_ENDIAN,OST_BIG_ENDIAN,VALUETYPE>::ConvertIP(value);
  }
};

template<class VALUETYPE>
class Convert<OST_LITTLE_ENDIAN,VALUETYPE>{
public:
  static VALUETYPE From(const VALUETYPE& value)
  {
    return ConvertHelper<OST_LITTLE_ENDIAN,OST_LOCAL_ENDIAN,VALUETYPE>::Convert(value);
  }
  static VALUETYPE To(const VALUETYPE& value)
  {
    return ConvertHelper<OST_LOCAL_ENDIAN,OST_LITTLE_ENDIAN,VALUETYPE>::Convert(value);
  }
  static void FromIP(VALUETYPE* value)
  {
    ConvertHelper<OST_LITTLE_ENDIAN,OST_LOCAL_ENDIAN,VALUETYPE>::ConvertIP(value);
  }
  static void ToIP(VALUETYPE* value)
  {
    ConvertHelper<OST_LOCAL_ENDIAN,OST_LITTLE_ENDIAN,VALUETYPE>::ConvertIP(value);
  }
};


template<class VALUETYPE>
class Convert<OST_VAX_DATA,VALUETYPE>{
public:
  static VALUETYPE From(const VALUETYPE& value)
  {
    return ConvertHelper<OST_VAX_DATA,OST_LOCAL_ENDIAN,VALUETYPE>::Convert(value);
  }
  static VALUETYPE To(const VALUETYPE& value)
  {
    return ConvertHelper<OST_LOCAL_ENDIAN,OST_VAX_DATA,VALUETYPE>::Convert(value);
  }
  static void FromIP(VALUETYPE* value)
  {
    ConvertHelper<OST_VAX_DATA,OST_LOCAL_ENDIAN,VALUETYPE>::ConvertIP(value);
  }
  static void ToIP(VALUETYPE* value)
  {
    ConvertHelper<OST_LOCAL_ENDIAN,OST_VAX_DATA,VALUETYPE>::ConvertIP(value);
  }
};

template class Convert<OST_BIG_ENDIAN,int8_t>;
template class Convert<OST_BIG_ENDIAN,uint8_t>;
template class Convert<OST_BIG_ENDIAN,int16_t>;
template class Convert<OST_BIG_ENDIAN,uint16_t>;
template class Convert<OST_BIG_ENDIAN,int32_t>;
template class Convert<OST_BIG_ENDIAN,uint32_t>;
template class Convert<OST_BIG_ENDIAN,int64_t>;
template class Convert<OST_BIG_ENDIAN,uint64_t>;
template class Convert<OST_BIG_ENDIAN,float>;
template class Convert<OST_BIG_ENDIAN,double>;

template class Convert<OST_LITTLE_ENDIAN,int8_t>;
template class Convert<OST_LITTLE_ENDIAN,uint8_t>;
template class Convert<OST_LITTLE_ENDIAN,int16_t>;
template class Convert<OST_LITTLE_ENDIAN,uint16_t>;
template class Convert<OST_LITTLE_ENDIAN,int32_t>;
template class Convert<OST_LITTLE_ENDIAN,uint32_t>;
template class Convert<OST_LITTLE_ENDIAN,int64_t>;
template class Convert<OST_LITTLE_ENDIAN,uint64_t>;
template class Convert<OST_LITTLE_ENDIAN,float>;
template class Convert<OST_LITTLE_ENDIAN,double>;

template class Convert<OST_VAX_DATA,int8_t>;
template class Convert<OST_VAX_DATA,uint8_t>;
template class Convert<OST_VAX_DATA,int16_t>;
template class Convert<OST_VAX_DATA,uint16_t>;
template class Convert<OST_VAX_DATA,int32_t>;
template class Convert<OST_VAX_DATA,uint32_t>;
template class Convert<OST_VAX_DATA,int64_t>;
template class Convert<OST_VAX_DATA,uint64_t>;
template class Convert<OST_VAX_DATA,float>;
template class Convert<OST_VAX_DATA,double>;

}} // ns




