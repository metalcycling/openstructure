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
#include <algorithm>
#include "atom_prop.hh"

namespace ost { namespace mol { namespace impl {

std::map<String, AtomProp> AtomProp::defaults_;

void AtomProp::Init()
{
  defaults_["AC"]=AtomProp(2.0000, 227.0000);
  defaults_["AG"]=AtomProp(1.7200, 107.8680);
  defaults_["AL"]=AtomProp(2.0000, 26.9820);
  defaults_["AM"]=AtomProp(2.0000, 243.0000);
  defaults_["AR"]=AtomProp(1.8800, 39.9480);
  defaults_["AS"]=AtomProp(1.8500, 74.9220);
  defaults_["AT"]=AtomProp(2.0000, 210.0000);
  defaults_["AU"]=AtomProp(1.6600, 196.9670);
  defaults_["B"]=AtomProp(2.0000, 10.8110);
  defaults_["BA"]=AtomProp(2.0000, 137.3270);
  defaults_["BE"]=AtomProp(2.0000, 9.0120);
  defaults_["BH"]=AtomProp(2.0000, 264.0000);
  defaults_["BI"]=AtomProp(2.0000, 208.9800);
  defaults_["BK"]=AtomProp(2.0000, 247.0000);
  defaults_["BR"]=AtomProp(1.8500, 79.9040);
  defaults_["C"]=AtomProp(1.7000, 12.0110);
  defaults_["CA"]=AtomProp(2.0000, 40.0780);
  defaults_["CD"]=AtomProp(1.5800, 112.4110);
  defaults_["CE"]=AtomProp(2.0000, 140.1160);
  defaults_["CF"]=AtomProp(2.0000, 251.0000);
  defaults_["CL"]=AtomProp(1.7500, 35.4530);
  defaults_["CM"]=AtomProp(2.0000, 247.0000);
  defaults_["CO"]=AtomProp(2.0000, 58.9330);
  defaults_["CR"]=AtomProp(2.0000, 51.9960);
  defaults_["CS"]=AtomProp(2.0000, 132.9050);
  defaults_["D"]=AtomProp(1.0900, 2.0141);
  defaults_["CU"]=AtomProp(1.4000, 63.5460);
  defaults_["DB"]=AtomProp(2.0000, 262.0000);
  defaults_["DS"]=AtomProp(2.0000, 271.0000);
  defaults_["DY"]=AtomProp(2.0000, 162.5000);
  defaults_["ER"]=AtomProp(2.0000, 167.2600);
  defaults_["ES"]=AtomProp(2.0000, 252.0000);
  defaults_["EU"]=AtomProp(2.0000, 151.9640);
  defaults_["F"]=AtomProp(1.4700, 18.9980);
  defaults_["FE"]=AtomProp(2.0000, 55.8450);
  defaults_["FM"]=AtomProp(2.0000, 257.0000);
  defaults_["FR"]=AtomProp(2.0000, 223.0000);
  defaults_["GA"]=AtomProp(1.8700, 69.7230);
  defaults_["GD"]=AtomProp(2.0000, 157.2500);
  defaults_["GE"]=AtomProp(2.0000, 72.6100);
  defaults_["H"]=AtomProp(1.0900, 1.0080);
  defaults_["HE"]=AtomProp(1.4000, 4.0030);
  defaults_["HF"]=AtomProp(2.0000, 178.4900);
  defaults_["HG"]=AtomProp(1.5500, 200.5900);
  defaults_["HO"]=AtomProp(2.0000, 164.9300);
  defaults_["HS"]=AtomProp(2.0000, 269.0000);
  defaults_["I"]=AtomProp(1.9800, 126.9040);
  defaults_["IN"]=AtomProp(1.9300, 114.8180);
  defaults_["IR"]=AtomProp(2.0000, 192.2170);
  defaults_["K"]=AtomProp(2.7500, 39.0980);
  defaults_["KR"]=AtomProp(2.0200, 83.8000);
  defaults_["LA"]=AtomProp(2.0000, 138.9060);
  defaults_["LI"]=AtomProp(1.8200, 6.9410);
  defaults_["LR"]=AtomProp(1.5000, 103.0000);
  defaults_["LU"]=AtomProp(2.0000, 174.9670);
  defaults_["MD"]=AtomProp(2.0000, 258.0000);
  defaults_["MG"]=AtomProp(1.7300, 24.3050);
  defaults_["MN"]=AtomProp(2.0000, 54.9380);
  defaults_["MO"]=AtomProp(2.0000, 95.9400);
  defaults_["MT"]=AtomProp(2.0000, 268.0000);
  defaults_["N"]=AtomProp(1.5500, 14.0070);
  defaults_["NA"]=AtomProp(2.2700, 22.9910);
  defaults_["NB"]=AtomProp(2.0000, 92.9060);
  defaults_["ND"]=AtomProp(2.0000, 144.2400);
  defaults_["NE"]=AtomProp(1.5400, 20.1800);
  defaults_["NI"]=AtomProp(1.6300, 58.6930);
  defaults_["NO"]=AtomProp(2.0000, 259.0000);
  defaults_["NP"]=AtomProp(2.0000, 237.0000);
  defaults_["O"]=AtomProp(1.5200, 15.9990);
  defaults_["OS"]=AtomProp(2.0000, 190.2300);
  defaults_["P"]=AtomProp(1.8000, 30.9740);
  defaults_["PA"]=AtomProp(2.0000, 231.0360);
  defaults_["PB"]=AtomProp(2.0200, 207.2000);
  defaults_["PD"]=AtomProp(1.6300, 106.4200);
  defaults_["PM"]=AtomProp(2.0000, 145.0000);
  defaults_["PO"]=AtomProp(2.0000, 210.0000);
  defaults_["PR"]=AtomProp(2.0000, 140.9080);
  defaults_["PT"]=AtomProp(1.7200, 195.0780);
  defaults_["PU"]=AtomProp(2.0000, 244.0000);
  defaults_["RA"]=AtomProp(2.0000, 226.0000);
  defaults_["RB"]=AtomProp(2.0000, 85.4680);
  defaults_["RE"]=AtomProp(2.0000, 186.2070);
  defaults_["RF"]=AtomProp(2.0000, 261.0000);
  defaults_["RH"]=AtomProp(2.0000, 102.9060);
  defaults_["RN"]=AtomProp(2.0000, 222.0000);
  defaults_["RU"]=AtomProp(2.0000, 101.0700);
  defaults_["S"]=AtomProp(1.8000, 32.0660);
  defaults_["SB"]=AtomProp(2.0000, 121.7600);
  defaults_["SC"]=AtomProp(2.0000, 44.9560);
  defaults_["SE"]=AtomProp(1.9000, 78.9600);
  defaults_["SG"]=AtomProp(2.0000, 266.0000);
  defaults_["SI"]=AtomProp(2.1000, 28.0860);
  defaults_["SM"]=AtomProp(2.0000, 150.3600);
  defaults_["SN"]=AtomProp(2.1700, 118.7100);
  defaults_["SR"]=AtomProp(2.0000, 87.6200);
  defaults_["TA"]=AtomProp(2.0000, 180.9480);
  defaults_["TB"]=AtomProp(2.0000, 158.9250);
  defaults_["TC"]=AtomProp(2.0000, 98.0000);
  defaults_["TE"]=AtomProp(2.0600, 127.6000);
  defaults_["TH"]=AtomProp(2.0000, 232.0380);
  defaults_["TI"]=AtomProp(2.0000, 47.8670);
  defaults_["TL"]=AtomProp(1.9600, 204.3830);
  defaults_["TM"]=AtomProp(2.0000, 168.9340);
  defaults_["U"]=AtomProp(1.8600, 238.0290);
  defaults_["V"]=AtomProp(2.0000, 50.9420);
  defaults_["W"]=AtomProp(2.0000, 183.8400);
  defaults_["XE"]=AtomProp(2.1600, 131.2900);
  defaults_["Y"]=AtomProp(2.0000, 88.9060);
  defaults_["YB"]=AtomProp(2.0000, 173.0400);
  defaults_["ZN"]=AtomProp(1.3900, 65.3900);
  defaults_["ZR"]=AtomProp(2.0000, 91.2240);

  defaults_[""]=AtomProp(0.0, 0.0);

}

AtomProp* AtomProp::GetDefaultProps(const String& ele)  
{
  if (defaults_.empty()) {
    AtomProp::Init();
  }
  String upper_ele=ele;
  std::transform(upper_ele.begin(),upper_ele.end(),upper_ele.begin(),toupper);  
  std::map<String, AtomProp>::iterator i=defaults_.find(upper_ele);
  if (i!=defaults_.end()) {
    return &i->second;
  }
  return &defaults_[""];
}


}}}
