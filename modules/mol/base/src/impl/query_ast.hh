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
#ifndef OST_QUERY_AST_HH
#define OST_QUERY_AST_HH

#include <boost/regex.hpp>

#include <ost/mol/module_config.hh>
#include <boost/variant.hpp>
#include <ost/geom/vec3.hh>
#include <ost/mol/view_type_fw.hh>
#include <ost/mol/property_id.hh>

#include "query_ast_fw.hh"

namespace ost { namespace mol { namespace impl {

/// holds the right and left hand operand of a within statement. The lefthand
/// operand is a distance in Angstrom, the righthand parameter may both be 
/// a lazily bound reference to a point cloud (that is only known at 
/// instantiation time) or a vector.
class DLLEXPORT_OST_MOL WithinParam {
public:
  WithinParam(const geom::Vec3& center, float radius);
  WithinParam(int ref, float radius);
  WithinParam();

  float GetRadiusSquare() const;

  bool operator==(const WithinParam& p) const;

  const geom::Vec3& GetCenter() const;
  int GetRef() const;
  bool HasValidRef() const;
private:
  geom::Vec3    center_;
  float         radius_;
  int           lazily_bound_ref_;
};

// holds either a simple string or a full regex, constructed from a glob-style string
class DLLEXPORT_OST_MOL StringOrRegexParam {
 public:
  StringOrRegexParam();
  explicit StringOrRegexParam(const String& s);
  bool Match(const String& s) const;
  bool operator==(const StringOrRegexParam&) const;
  const String& str() const {return s_;}
 private:
  bool is_regex_;
  boost::regex r_;
  String s_;
};

typedef boost::variant<int, Real, WithinParam, StringOrRegexParam> ParamType;

// AST node, used internally for building the AST tree.
class DLLEXPORT_OST_MOL Node {
public:
  Node(): parent_(NULL) { }
  virtual ~Node() { }
  virtual void Dump(int level=0) const = 0;

  Node* GetParent();

  void SetParent(Node* parent);
private:
    Node*           parent_;
};



class DLLEXPORT_OST_MOL LogicOPNode : public Node {
public:
  LogicOPNode(LogicOP op);

  ~LogicOPNode();

  //! set right-hane side node
  void SetRHS(Node* rhs);

  //! set left hand-side node
  void SetLHS(Node* lhs);

  //! Get right-hand side
  const Node* GetRHS() const;

  //! Get left-hand side argument
  const Node* GetLHS() const;


  //! Get logic operation
  LogicOP GetOP() const;

  void SetOP(LogicOP op);

  virtual void Dump(int level=0) const;

private:
  Node*           lhs_;
  Node*           rhs_;
  LogicOP         op_;
};


class DLLEXPORT_OST_MOL SelNode : public Node {
public:
  SelNode(const Prop& sel, CompOP op, const ParamType& value)
    : sel_(sel), op_(op), param_(value)
  { }

  SelNode(const SelNode& rhs)
    : sel_(rhs.sel_), op_(rhs.op_),param_(rhs.param_)
  { }

  virtual void Dump(int level=0) const;

  const Prop& GetAtomProps() const {
    return sel_;
  }
  CompOP GetCompOP() const {
    return op_;
  }
  const ParamType& GetParm() const {
    return param_;
  }
private:
  Prop     sel_;
  CompOP       op_;
  ParamType    param_;

};

}}} // ns

#endif
