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
/*
 Author: Marco Biasini
 */
#include <iostream>
#include "query_ast.hh"

namespace ost { namespace mol { namespace impl {
  

LogicOPNode::LogicOPNode(LogicOP op) 
  : lhs_(NULL), rhs_(NULL), op_(op) 
{}

LogicOPNode::~LogicOPNode() {
  if (rhs_)
    delete rhs_;
  if (lhs_)
    delete lhs_;
}

void LogicOPNode::SetRHS(Node* rhs) {
  rhs_ = rhs;
}

void LogicOPNode::SetLHS(Node* lhs) {
  lhs_ = lhs;
}
const Node* LogicOPNode::GetRHS() const {
  return rhs_;
}

const Node* LogicOPNode::GetLHS() const {
  return lhs_;
}

LogicOP LogicOPNode::GetOP() const {
  return op_;
}

void LogicOPNode::SetOP(LogicOP op) {
  op_ = op;
}

    
void LogicOPNode::Dump(int level) const{
  for(int i=0; i < level; ++i) 
    std::cout << "\t";
  String op = "";
  switch (op_) {
    case LOP_OR:
      op = "or";
      break;
    case LOP_AND:
      op = "and";
      break;
  }
  std::cout << op<< std::endl;
  if (lhs_ && rhs_) {
    lhs_->Dump(level+1);
    rhs_->Dump(level+1);
  }
}

WithinParam::WithinParam(int ref, float radius) 
  : radius_(radius), lazily_bound_ref_(ref) {
}

WithinParam::WithinParam(const geom::Vec3& center, float radius)
  : center_(center), radius_(radius), lazily_bound_ref_(-1) {
}
int WithinParam::GetRef() const {
  return lazily_bound_ref_;
}

bool WithinParam::HasValidRef() const {
  return lazily_bound_ref_!=-1;
}

WithinParam::WithinParam()
  : center_(geom::Vec3(0,0,0)), radius_(0.0), lazily_bound_ref_(-1) {
}

bool WithinParam::operator==(const WithinParam& p) const {
  if (lazily_bound_ref_>-1) {
    return lazily_bound_ref_==p.lazily_bound_ref_ && p.radius_==radius_;
  }
  return center_==p.center_;
}

float WithinParam::GetRadiusSquare() const {
  return radius_;
}

const geom::Vec3& WithinParam::GetCenter() const {
  return center_;
}

StringOrRegexParam::StringOrRegexParam():
  is_regex_(false), r_(), s_()
{}

StringOrRegexParam::StringOrRegexParam(const String& s):
  is_regex_(false),r_(),s_(s)
{
  String special("[]{}()");
  for(String::const_iterator it=s.begin();it!=s.end();++it) {
    if((*it)=='?' || (*it)=='*') {
      is_regex_=true;
      break;
    }
  }

  if(is_regex_) {
    std::ostringstream e;
    for(String::const_iterator it=s.begin();it!=s.end();++it) {
      
      if((*it)=='?' && (it==s.begin() || (*(it-1))!='\\')) {
        e << ".";
      } else if((*it)=='*' && (it==s.begin() || (*(it-1))!='\\')) {
        e << ".*";
      } else {
        if (special.find(*it)!=String::npos) {
          e << '\\';
        }
        e << *it;
      }
    }
    //std::cerr << "assembling regex [" << e.str() << "]... ";
    r_=boost::regex(e.str());
    //std::cerr << "done" << std::endl;
  }
}

bool StringOrRegexParam::Match(const String& s) const
{
  if(is_regex_) {
    //std::cerr << "matching regex to [" << s << "]... ";
    bool result = boost::regex_match(s,r_);
    //std::cerr << "done" << std::endl;
    return result;
  }
  return s==s_;
}

bool StringOrRegexParam::operator==(const StringOrRegexParam& o) const
{
  if(is_regex_) return r_==o.r_;
  return s_==o.s_;
}

Node* Node::GetParent() {
  return parent_;
}

void Node::SetParent(Node* parent) {
  parent_=parent;
}

String CompOPToString(CompOP op) {
  switch (op) {
    case COP_EQ:
      return "=";
    case COP_NEQ:
      return "!=";
    case COP_GE:
      return ">=";
    case COP_LE:
      return "<=";
    case COP_GT:
      return ">";
    case COP_LT:
      return "<";                             
  }
  // Keep the compiler happy
  return "";
}

void SelNode::Dump(int level) const{
  for(int i=0; i < level; ++i) 
    std::cout << "\t";
  String op = CompOPToString(op_);
  switch(sel_.type) {
    case Prop::STRING:
      std::cout << sel_.GetName()<< " " << op << " " << boost::get<StringOrRegexParam>(param_).str() << std::endl;
      break;
  case Prop::FLOAT:
    std::cout << sel_.GetName()<< " " << op << " " << boost::get<float>(param_) << std::endl;
    break;        
  case Prop::INT:
    std::cout << sel_.GetName()<< " " << op << " " << boost::get<int>(param_) << std::endl;
    break;      
  case Prop::VEC_DIST:
    std::cout << "{" << (boost::get<WithinParam>(param_).GetCenter())[0] << ", "
              << (boost::get<WithinParam>(param_).GetCenter())[1] << ", "
              << (boost::get<WithinParam>(param_).GetCenter())[2] <<"} <> "
              << sqrt(boost::get<WithinParam>(param_).GetRadiusSquare())
              << std::endl;
    break;
  }
}

}}} // ns
