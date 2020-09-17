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
#include "query_impl.hh"
#include <ctype.h>

#include <boost/version.hpp>
#if BOOST_VERSION>=103800
#define BOOST_SPIRIT_USE_OLD_NAMESPACE
#include <boost/spirit/include/classic_symbols.hpp>
#else
#include <boost/spirit/symbols.hpp>
#endif

#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/chain_impl.hh>
#include <ost/mol/mol.hh>

#include <memory>
#include <cstddef>

namespace ost { namespace mol {  namespace impl {

// portable solution for use of unique ptr in this file
// -> use C++11 unique_ptr if possible and old auto_ptr otherwise
#if __cplusplus >= 201103L
typedef std::unique_ptr<Node> NodeUniquePtr;
#else
typedef std::auto_ptr<Node> NodeUniquePtr;
#endif

using namespace boost::spirit;

namespace {
  
class IdentTokens : public symbols<tok::Type> {
public:  
  IdentTokens() {
    add("and", tok::And)
       ("or", tok::Or)
       ("not", tok::Not);
  }
};

} // anon namespace


QueryToken::QueryToken() 
  : type_(tok::None), range_() {
}

QueryToken::QueryToken(const Range& range, tok::Type type, 
                       const Range& value_range)
  : type_(type), range_(range) {
  if (value_range==Range()) {
    value_range_=range;
  } else {
    value_range_=value_range;
  }
}

tok::Type QueryToken::GetType() const {
  return type_;
}


  //! Returns true if token is end of file token.
bool QueryToken::IsEOF() const {
  return type_==tok::None;
}

const Range& QueryToken::GetRange() const {
  return range_;
}

const Range& QueryToken::GetValueRange() const {
  return value_range_;
}

QueryLexer::QueryLexer(const String& query_string)
  : query_string_(query_string), current_(0) 
{ }

QueryToken QueryLexer::LexQuotedStringLiteral()
{
  char quote=query_string_[current_];
  Range range(current_+1, 0);
  while(++current_ && current_<query_string_.size()) {
    if (query_string_[current_]==quote && query_string_[current_-1]!='\\') {
      current_++;
      return QueryToken(Range(range.Loc-1, current_-range.Loc+1), tok::String,
                        Range(range.Loc, current_-range.Loc-1));
    }
  }
  return QueryToken(Range(range.Loc-1, 1), tok::UnterminatedQuote, 
                    Range(range.Loc, 0));
}

QueryToken QueryLexer::LexNumericToken() {
  Range range(current_, 0);
  static String allowed_chars("_");  
  bool dot_present=false, is_string=false, last_was_dot=false;
  //TODO Better checking for duplicate -
  while(current_<query_string_.size()){
    if (isdigit(query_string_[current_]) || query_string_[current_]=='-') {
      current_++;
      last_was_dot=false;
    } else if (query_string_[current_]=='.') {
      if (dot_present) {
        is_string=true;
        break;
      }
      dot_present = true;
      last_was_dot=true;
      current_++;
    } else {
      if (last_was_dot) {
        break;
      }
      is_string=isalnum(query_string_[current_]) || 
                allowed_chars.find_first_of(query_string_[current_])!=String::npos;
      break;
    }
  }
  range.Length=current_-range.Loc;
  if (!is_string) {
    if (dot_present) {
      if (last_was_dot) {
        current_--;
        return QueryToken(Range(range.Loc, range.Length-1), 
                          tok::IntegralValue);
      }
      return QueryToken(range, tok::FloatingValue);      
    } else {
      return QueryToken(range, tok::IntegralValue);      
    }    
  }
  current_=range.Loc;
  return this->LexIdentOrStringToken();
}

bool is_ident_or_str(char c) {
  static String allowed_chars("_*?\\");
  return isalnum(c) || allowed_chars.find_first_of(c)!=String::npos;
}

QueryToken QueryLexer::LexIdentOrStringToken() {
  static IdentTokens ident_tokens;
  size_t start=current_;
  bool force_string=false;
  while (current_<query_string_.length() && 
         is_ident_or_str(query_string_[current_])) {
    if (query_string_[current_]=='*' || query_string_[current_]=='?' ||
        query_string_[current_]=='\\') {
      force_string=true;
    }
    current_++;
  }
  String ident=query_string_.substr(start, current_-start);
  if (tok::Type* t=find(ident_tokens, ident.c_str())) {
    return QueryToken(Range(start, current_-start), *t);
  }
  return QueryToken(Range(start, current_-start), 
                    force_string? tok::String : tok::Identifier);
}

QueryToken QueryLexer::LexToken() {
  //read over whitespaces
  this->EatWhitespaces();
  if (current_ < query_string_.size()) {
    char current_char=query_string_[current_];
    if (isdigit(current_char) || current_char=='-') {
      return this->LexNumericToken();
    }
    if (isalpha(current_char) || current_char=='?' || current_char=='*' || 
        current_char=='_') {
      return this->LexIdentOrStringToken();
    }
    switch (current_char) {
      case ',':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::Coma);
      case ':':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::Colon);
      case '.':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::Dot);           
      case '{':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::LeftCurlyBrace);
      case '}':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::RightCurlyBrace);
      case '(':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::LeftParen);
      case ')':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::RightParen);
      case '[':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::LeftBracket);
      case '"':
      case '\'':
        return this->LexQuotedStringLiteral();
      case ']':
        current_++;
        return QueryToken(Range(current_-1, 1), tok::RightBracket);
      case '=':
        current_++;
        if (current_<query_string_.length() && 
            query_string_[current_]=='=') {
          current_++;
          return QueryToken(Range(current_-2, 2), tok::Equal);
        } else {
          return QueryToken(Range(current_-1, 1), tok::Equal);          
        }
      case '!':
        if (current_+1<query_string_.length() && 
            query_string_[current_+1]=='=') {
          current_+=2;
          return QueryToken(Range(current_-2, 2), tok::NotEqual);
        }
        break;
      case '>':
        if (current_+1<query_string_.length() && 
            query_string_[current_+1]=='=') {
          current_+=2;
          return QueryToken(Range(current_-2, 2), tok::GreaterEqualThan);
        } else {
          current_++;
          return QueryToken(Range(current_-1, 1), tok::GreaterThan);
        }
      case '<':
        if (current_+1<query_string_.length() && 
            query_string_[current_+1]=='=') {
          current_+=2;
          return QueryToken(Range(current_-2, 2), tok::LessEqualThan);
        } else if (current_+1<query_string_.length() &&
                   query_string_[current_+1]=='>') {
          current_+=2;
          return QueryToken(Range(current_-2, 2), tok::Within);          
        } else {
          current_++;
          return QueryToken(Range(current_-1, 1), tok::LessThan);
        }
      default:
        return QueryToken(Range(current_, 1), tok::UnknownChar);
    }
  }  
  // return empty EOF token.
  return QueryToken();
}
QueryToken QueryLexer::NextToken() {
  tok_=this->LexToken();
  return tok_;
}

QueryToken QueryLexer::CurrentToken() const {
  return tok_;
}

void QueryLexer::EatWhitespaces() {
  while (current_<query_string_.length() && isspace(query_string_[current_])) {
    ++current_;
  }
}  

namespace {
}

bool QueryImpl::IsAlwaysUndef(const Node* ast,
                              Prop::Level target_level) {
  const LogicOPNode* lop_node = dynamic_cast<const LogicOPNode*>(ast);
  if (lop_node) {
    bool lhs = this->IsAlwaysUndef(lop_node->GetLHS(), target_level);
    bool rhs = this->IsAlwaysUndef(lop_node->GetRHS(), target_level);    
    return lhs && rhs;
  }  
  else {
    const SelNode* sel_node = dynamic_cast<const SelNode*>(ast);
    assert(sel_node);
    if ((int)target_level > (int)sel_node->GetAtomProps().level) {
      return true;
    }
  }
  return false;
}

void QueryImpl::ExtractSelStmts(const Node* ast) {
  const LogicOPNode* lop_node = dynamic_cast<const LogicOPNode*>(ast);
  if (lop_node) {
    this->ExtractSelStmts(lop_node->GetLHS());
    this->ExtractSelStmts(lop_node->GetRHS());
  }
  else {
    const SelNode* sel_node = dynamic_cast<const SelNode*>(ast);
    assert(sel_node);
    SelStmt val(sel_node->GetAtomProps().id,
                sel_node->GetCompOP(),
                sel_node->GetParm());
    std::vector<SelStmt>::iterator i = std::find(sel_values_.begin(),
                                                 sel_values_.end(),val);
    if (i==sel_values_.end()) {
      size_t idx = sel_values_.size();
      sel_values_.push_back(val);
      indices_[(int)sel_node->GetAtomProps().level].insert(idx);      
    }
  }  
}

bool QueryImpl::IsValid() const 
{
  return !has_error_;
}

bool QueryImpl::IsEmpty() const
{
  return empty_optimize_;
}

QueryImpl::QueryImpl(Node* ast)
  : query_string_(), has_error_(false), num_gen_prop_(0) {
  sel_stacks_.resize(3);
  indices_.resize(3);
  // Add an undef SelStmt at index 0.
  SelStmt sel_value;
  sel_values_.push_back(sel_value);
  // ast_root->Dump();      
  // Get all selection statements
  this->ExtractSelStmts(ast);
  this->ASTToSelStack(ast, Prop::CHAIN,
                      sel_stacks_[(int)Prop::CHAIN]);
  this->ASTToSelStack(ast, Prop::RESIDUE,
                      sel_stacks_[(int)Prop::RESIDUE]);
  this->ASTToSelStack(ast, Prop::ATOM,
                      sel_stacks_[(int)Prop::ATOM]);        
  delete ast;
  empty_optimize_=false;
}

void QueryImpl::ASTToSelStack(const Node* src_ast,
                              Prop::Level target_level,
                              SelStack& stack) {
  if (this->IsAlwaysUndef(src_ast,target_level)) {
    // Push index 0 which always points to the UNDEF SelStmt
    SelItem undef;
    undef.type = VALUE;
    undef.value = 0;
    stack.push_back(undef);
    return;
  }
  const LogicOPNode* lop_n=dynamic_cast<const LogicOPNode*>(src_ast);
  if (lop_n) {
    SelItem logic;
    logic.type = LOGIC;
    logic.value = (size_t)lop_n->GetOP();
    stack.push_back(logic);    
    this->ASTToSelStack(lop_n->GetLHS(),target_level, stack);
    this->ASTToSelStack(lop_n->GetRHS(),target_level, stack);      
  }
  else {
    const SelNode* sel_node = dynamic_cast<const SelNode*>(src_ast);
    assert((int)target_level <= (int)sel_node->GetAtomProps().level);
    SelStmt val(sel_node->GetAtomProps().id,sel_node->GetCompOP(),
                sel_node->GetParm());
    std::vector<SelStmt>::iterator i = std::find(sel_values_.begin(),
                                                 sel_values_.end(),val);
    assert(i!=sel_values_.end());
    SelItem value;
    value.type = VALUE;
    value.value = i - sel_values_.begin();
    stack.push_back(value);    
  }
}

QueryImpl::QueryImpl(const String& query_string) 
  : query_string_(query_string), has_error_(false), num_gen_prop_(0) {
  if (query_string.length()>0) {
    sel_stacks_.resize(3);
    indices_.resize(3);
    // Add an undef SelStmt at index 0.
    SelStmt sel_value;
    sel_values_.push_back(sel_value);    
    Node* ast_root = this->BuildAST();
    if (ast_root) {
      // Get all selection statements
      this->ExtractSelStmts(ast_root);
      this->ASTToSelStack(ast_root,Prop::CHAIN,
                          sel_stacks_[(int)Prop::CHAIN]);
      this->ASTToSelStack(ast_root,Prop::RESIDUE,
                          sel_stacks_[(int)Prop::RESIDUE]);
      this->ASTToSelStack(ast_root,Prop::ATOM,
                          sel_stacks_[(int)Prop::ATOM]);        
      delete ast_root;
      empty_optimize_=false;
      has_error_=false;
    } else {
      has_error_=true;
    }  
  } else {
    empty_optimize_=true;
  }
}


const QueryErrorDesc& QueryImpl::GetErrorDescription() const {
  return error_desc_;
}


bool is_literal_true(const String& value_string) {
  return value_string=="true" || value_string == "True" || 
         value_string == "TRUE";
}

bool is_literal_false(const String& value_string) {
  return value_string=="false" || value_string == "False" || 
         value_string == "FALSE";
}

bool QueryImpl::ParseValue(const Prop& sel, const QueryToken& op,
                           QueryLexer& lexer, ParamType& value) {
  QueryToken last=lexer.CurrentToken();                               
  QueryToken v=lexer.NextToken();
  if (!this->ExpectNotEnd(v, "value"))
    return false;
  if (v.GetType()==tok::None) {
    error_desc_.msg="unexpected end of query string after '"+
                     query_string_.substr(v.GetRange().Loc, 
                                          v.GetRange().Length)+"'";
    error_desc_.range=v.GetRange();    
  }
  String value_string=query_string_.substr(v.GetValueRange().Loc, 
                                           v.GetValueRange().Length);
  switch (v.GetType()) {
    case tok::Identifier:
      if (sel.type==Prop::INT || sel.id>=Prop::CUSTOM) {
        // todo. Add check to test that the comparison operator is only one of
        // = and !=. The others don't make too much sense.
        if (is_literal_true(value_string)) {
          // these two if/else branches are unfortunate. but we need to cast
          // the bool to either int or float, depending whether the result
          // is to be used for generic properties (which are always floats)
          // or for integral properties.
          if (sel.type==Prop::INT) {
            value=ParamType(int(1));
          } else {
            value=ParamType(Real(1));
          }
          break;
        } else if (is_literal_false(value_string)) {   
          if (sel.type==Prop::INT) {
            value=ParamType(int(0));
          } else {
            value=ParamType(Real(0));
          }
          break;
        }
      }
      // yes, not having a break here is on purpose
    case tok::String:      
      if (sel.type!=Prop::STRING) {
        if (sel.id>=Prop::CUSTOM) {
          // BZDNG-204: issue specific warning when trying to use a string value 
          //     for a generic property.
          error_desc_.msg="only numeric generic properties can be used in queries";
          error_desc_.range=v.GetValueRange();
        } else {
          error_desc_.msg="'"+sel.GetName()+"' takes "+sel.GetTypeName()+
                         " argument, but string literal given";
          error_desc_.range=v.GetValueRange();
        }
        return false;
      } else {
        value=StringOrRegexParam(value_string);
      }

      break;
    case tok::FloatingValue:
      if (sel.type==Prop::INT) {
        error_desc_.msg="'"+sel.GetName()+"' takes "+sel.GetTypeName()+
                       " argument, but floating point value given";
        error_desc_.range=v.GetRange();
        return false;
      } else if (sel.type==Prop::STRING) {
        value=StringOrRegexParam(value_string);
      } else
        value=ParamType(float(atof(value_string.c_str())));              
      break;      
    case tok::IntegralValue:
      if (sel.type==Prop::STRING) {
        value=StringOrRegexParam(value_string);
      } else {
        if (sel.type==Prop::INT) {
          value=ParamType(atoi(value_string.c_str()));
        } else {
          value=ParamType(float(atof(value_string.c_str())));
        } 
      }      

      break;
    default:
      error_desc_.msg="expected value after '"+
        query_string_.substr(last.GetRange().Loc, 
                             last.GetRange().Length)+"', but '"+
        query_string_.substr(v.GetRange().Loc, 
                             v.GetRange().Length)+"' found";
        return false;
    }                  
    return true;
}

namespace {
LogicOP Invert(LogicOP op) {
  switch (op) {
    case LOP_AND:
      return LOP_OR;
    case LOP_OR:
      return LOP_AND;
  }
  assert(0);
  return op;
}

Range TotalRange(const Range& r1, const Range& r2) {
  Range r;
  r.Loc=r1.Loc;
  r.Length=r2.End()-r.Loc;
  return r;
}

CompOP GetCompOP(tok::Type op, bool invert) {
  CompOP comp_op=COP_EQ;
  if (!invert) {
    switch (op) {
      case tok::Equal:
        comp_op=COP_EQ;
        break;
      case tok::NotEqual:
        comp_op=COP_NEQ;
        break;
      case tok::LessEqualThan:
        comp_op=COP_LE;
        break;
      case tok::GreaterEqualThan:
        comp_op=COP_GE;
        break;
      case tok::LessThan:
        comp_op=COP_LT;
        break;                
      case tok::GreaterThan:
        comp_op=COP_GT;
        break;                                
      default:
       assert(0 && "should never get here.");
    }    
  } else {
    switch (op) {
      case tok::Equal:
        comp_op=COP_NEQ;
        break;
      case tok::NotEqual:
        comp_op=COP_EQ;
        break;
      case tok::LessEqualThan:
        comp_op=COP_GT;
        break;
      case tok::GreaterEqualThan:
        comp_op=COP_LT;
        break;
      case tok::LessThan:
        comp_op=COP_GE;
        break;                
      case tok::GreaterThan:
        comp_op=COP_LE;
        break;                                
      default:
       assert(0 && "should never get here.");
    }        
  }
  return comp_op;  
}

} // end anonymous namespace

Node* QueryImpl::ParseValueOrRange(const Prop& s, const QueryToken& op, 
                                   QueryLexer& lexer) {
  ParamType v;
  if (!this->ParseValue(s, op, lexer, v))
    return NULL;
                           
  // do look-ahead.
  QueryToken n=lexer.NextToken();
  if (n.GetType()==tok::Colon) {
    ParamType v2;
    if (!this->ParseValue(s, op, lexer, v2))
      return NULL;
    if (s.type==Prop::STRING) {
      error_desc_.msg="use of range with '"+s.GetName()+"' that takes "
                      "String argument is not allowed.";
      error_desc_.range=TotalRange(op.GetRange(), 
                                   lexer.CurrentToken().GetRange());
      return NULL;      
    }
    if (!(op.GetType()==tok::NotEqual || op.GetType()==tok::Equal)) {
      String o=query_string_.substr(op.GetValueRange().Loc, 
                                    op.GetRange().Length);
      error_desc_.msg="use of range with '"+o+
                      "' is not allowed.";
      error_desc_.range=TotalRange(op.GetRange(), 
                                   lexer.CurrentToken().GetRange());
      return NULL;
    } 
    // advance to next token. this makes sure that the caller of this method 
    // does not have to care whether we had a value or value-range.
    lexer.NextToken();    
    CompOP left_op, right_op;     
    LogicOP lop;
    CompOP o=GetCompOP(op.GetType(), inversion_stack_.back());
    if (o==COP_EQ) {
      left_op=COP_GE;
      right_op=COP_LE;
      lop=LOP_AND;
    } else {
      left_op=COP_LT;
      right_op=COP_GT;
      lop=LOP_OR;
    }
    SelNode* sn1=new SelNode(s, left_op, v);
    SelNode* sn2=new SelNode(s, right_op, v2);
    return this->Concatenate(sn1, sn2, lop);
  } else {
    return new SelNode(s, GetCompOP(op.GetType(), inversion_stack_.back()), v);
  }
}

Node* QueryImpl::ParsePropValueExpr(QueryLexer& lexer) {
  QueryToken sname=lexer.CurrentToken();
  String s=query_string_.substr(sname.GetRange().Loc,sname.GetRange().Length);
  Prop property=PropertyFromString(s);
  QueryToken op=lexer.NextToken();
  
  // this block deals with the cname.rnum.aname shortcut.
  if (op.GetType()==tok::Dot) {
    QueryToken cname=sname;
    QueryToken rnum=lexer.NextToken();
    if (!this->Expect(tok::IntegralValue, "residue number", rnum)) {
      return NULL;
    }
    if (!this->Expect(tok::Dot, "'.'", lexer.NextToken())) {
      return NULL;
    }
    QueryToken aname=lexer.NextToken();
    if (!this->Expect(tok::Identifier, "atom name", aname)) {
      return NULL;
    }
    LogicOP lop=inversion_stack_.back() ? LOP_OR : LOP_AND;
    CompOP cop=inversion_stack_.back() ? COP_NEQ : COP_EQ;
    ParamType cname_val(StringOrRegexParam(query_string_.substr(cname.GetValueRange().Loc,
                                                                cname.GetValueRange().Length).c_str()));
    Prop cname_prop(Prop::CNAME, Prop::STRING, Prop::CHAIN);
    SelNode* cname_node=new SelNode(cname_prop, cop, cname_val);
    ParamType aname_val(StringOrRegexParam(query_string_.substr(aname.GetValueRange().Loc,
                                                                aname.GetValueRange().Length).c_str()));
    Prop aname_prop(Prop::ANAME, Prop::STRING, Prop::ATOM);
    SelNode* aname_node=new SelNode(aname_prop, cop, aname_val);
    ParamType rnum_val(atoi(query_string_.substr(rnum.GetValueRange().Loc,
                            rnum.GetValueRange().Length).c_str()));
    Prop rnum_prop(Prop::RNUM, Prop::INT, Prop::RESIDUE);                            
    SelNode* rnum_node=new SelNode(rnum_prop, cop, rnum_val);
    LogicOPNode* and_one=new LogicOPNode(lop);
    LogicOPNode* and_two=new LogicOPNode(lop);  
    and_one->SetLHS(and_two);
    and_one->SetRHS(cname_node);
    and_two->SetLHS(rnum_node);
    and_two->SetRHS(aname_node);
    lexer.NextToken();
    return and_one;
  }
  if (property.id>=Prop::CUSTOM) {
    property.id=(Prop::ID)(Prop::CUSTOM+num_gen_prop_++);
    EntityPropertyMapper epm=EntityPropertyMapper(s.substr(2), property.level);
    GenProp gen_prop(epm);
    if (op.GetType()==tok::Colon) {
      op=lexer.NextToken();
      String value_str=query_string_.substr(op.GetValueRange().Loc,
                                            op.GetValueRange().Length);
      if (!this->ExpectNotEnd(op, "value")) {
        return NULL;
      }
      bool numeric = tok::FloatingValue==op.GetType() ||
                     tok::IntegralValue==op.GetType() ||
                     (tok::Identifier==op.GetType() &&
                     (is_literal_true(value_str) ||
                      is_literal_false(value_str)));
      if (!numeric) {
        // BZDNG-204: issue specific warning when trying to use a string value 
        //     for a generic property. 
        if (op.GetType()==tok::String || op.GetType()==tok::Identifier) {
          error_desc_.msg="only numeric generic properties are supported in queries";
          error_desc_.range=op.GetValueRange();
          return NULL;
        }
        this->ExpectNumeric(op);
        return NULL;
      }
      if (is_literal_true(value_str)) {
        gen_prop.default_val=1.0;
      } else if (is_literal_false(value_str)) {
        gen_prop.default_val=0.0;
      } else {
        gen_prop.default_val=atof(value_str.c_str());
      }
      gen_prop.has_default=true;
      op=lexer.NextToken();
    }
    gen_prop_list_.push_back(gen_prop);
  }
  if (property.id==Prop::UNDEF) {
    error_desc_.msg="'"+s+"' is not a valid property";
    error_desc_.range=sname.GetRange();
    return NULL;
  }
  if (!this->ExpectComparisonOperator(op))
    return NULL;  
  // in case of String properties, the comparison operators <= < >, >= are
  // not allowed.
  if (property.type==Prop::STRING) {
    if (!(op.GetType()==tok::NotEqual || op.GetType()==tok::Equal)) {
      String o=query_string_.substr(op.GetRange().Loc, op.GetRange().Length);
      error_desc_.msg="use of '"+s+"' with '"+o+
                      "' that takes String argument is not allowed";
      error_desc_.range=op.GetRange();
      return NULL;
    }
  }
  NodeUniquePtr root_value;
  do {
    Node* sel_node=this->ParseValueOrRange(property, op, lexer);
    if (sel_node==NULL) {
      return NULL;
    }
    // do some tree-rearrangement foo
    if (root_value.get()==NULL) {
      root_value.reset(sel_node);
    } else {
      LogicOP lop=op.GetType()==tok::NotEqual ? LOP_AND : LOP_OR;
      if (inversion_stack_.back())
        lop=Invert(lop);      
      Node* new_root=this->Concatenate(root_value.get(), sel_node, lop);
      // first release then set root
      root_value.release(); root_value.reset(new_root);
    }
    if (lexer.CurrentToken().GetType()==tok::Coma) {
      if (!(op.GetType()==tok::Equal || op.GetType()==tok::NotEqual)) {
        String o=query_string_.substr(op.GetRange().Loc, op.GetRange().Length);
        error_desc_.msg="use of ',' with '"+o+"' does not make any sense.";
        error_desc_.range=lexer.CurrentToken().GetRange();
        return NULL;
      }
    }
  } while (lexer.CurrentToken().GetType()==tok::Coma);
  return root_value.release();
}

const String& QueryImpl::GetQueryString() const {
  return query_string_;
}

Node* QueryImpl::ParseParenSubExpr(QueryLexer& lexer) {
  QueryToken paren_token=lexer.CurrentToken();
  lexer.NextToken();
  Node* n=this->ParseSubExpr(lexer, true);
  if (n) {

    QueryToken end_paren=lexer.CurrentToken();
    if (end_paren.GetType()!=tok::RightParen) {
      delete n;      
      error_desc_.msg="unmatched parenthesis in query string";
      error_desc_.range=paren_token.GetRange();
      return NULL;
    } else {
      lexer.NextToken();          
    }
  }
  return n;
}

Node* QueryImpl::ParseBracketSubExpr(QueryLexer& lexer) {
  QueryToken bracket_token=lexer.CurrentToken();
  lexer.NextToken();
  Node* n=this->ParseSubExpr(lexer, true);
  if (n) {
    QueryToken end_bracket=lexer.CurrentToken();
    if (end_bracket.GetType()!=tok::RightBracket) {
      delete n;      
      error_desc_.msg="unmatched bracket in query string";
      error_desc_.range=bracket_token.GetRange();
      return NULL;
    } else {
      lexer.NextToken();          
    }
  }
  return n;
}

Node* QueryImpl::ParseSubExpr(QueryLexer& lexer, bool paren) {
  // Use unique ptr to protect against memory leaks. In case of early exit, i.e.
  // exceptions, the memory associated with root_node and its child nodes will
  // be freed.
  NodeUniquePtr root_node;
  Node* new_node=NULL;
  QueryToken t=lexer.CurrentToken();
  LogicOP logic_op=LOP_AND;
  while (!t.IsEOF()) {
    if (t.GetType()==tok::Not) {
      inversion_stack_.push_back(!inversion_stack_.back());
      t=lexer.NextToken();      
    } else {
      inversion_stack_.push_back(inversion_stack_.back());
    }
    if (!this->ExpectNotEnd(t, "predicate"))
      return NULL;
    switch (t.GetType()) {
      case tok::Identifier:
        new_node=this->ParsePropValueExpr(lexer);
        if (new_node==NULL) {
          return NULL;
        }
        break;
      case tok::LeftParen:
        new_node=this->ParseParenSubExpr(lexer);
        if (new_node==NULL) {
          return NULL;
        }
        break;
      case tok::IntegralValue:
      case tok::FloatingValue:
        new_node=this->ParseWithinExpr(lexer);
        break;  
      default:
        error_desc_.msg="unexpected '"+
          query_string_.substr(t.GetRange().Loc, 
                               t.GetRange().Length)+"' in query string";
        error_desc_.range=t.GetRange();
        return NULL;
    }
    if (!new_node) {
      return NULL;
    }
    // transform AST to concatenate the sub expressions.
    if (root_node.get()==NULL) {
      root_node.reset(new_node);
    } else {
      Node* new_root=this->Concatenate(root_node.get(), new_node, logic_op);
      // release first, then set. otherwise root node is freed.
      root_node.release(); root_node.reset(new_root);
    }
    inversion_stack_.pop_back();    
    t=lexer.CurrentToken();    
    if (t.GetType()==tok::RightParen||t.GetType()==tok::RightBracket) {
      
      if (paren) {
        return root_node.release();
      }
      error_desc_.msg="unexpected '"+query_string_.substr(t.GetRange().Loc, 
        t.GetRange().Length)+"' in query string";
      error_desc_.range=t.GetRange();
      return NULL;
    }
    if (t.IsEOF()) {
      return root_node.release();
    }
    if (!this->ExpectLogicalOperator(t))
      return NULL;
    switch(t.GetType()) {
      case tok::Or:
        logic_op=LOP_OR;
        break;
      case tok::And:
        logic_op=LOP_AND;
        break;
      default:
        assert(0 && "what the heck?");
    }
    if (inversion_stack_.back())
      logic_op=Invert(logic_op);
    t=lexer.NextToken();
    if (!this->ExpectNotEnd(t, "predicate"))
      return NULL;
  }
  return root_node.release();  
}

bool QueryImpl::ExpectLogicalOperator(const QueryToken& token) {
  if (!(token.GetType()==tok::And || token.GetType()==tok::Or)) {
    if (token.GetType()==tok::UnterminatedQuote) {
      error_desc_.msg="unterminated quote in query string";
      error_desc_.range=token.GetRange();      
    } else {
      error_desc_.msg="expected logical operator, but '"+
                      query_string_.substr(token.GetRange().Loc, 
                                           token.GetRange().Length)+"' found";
      error_desc_.range=token.GetRange();      
    }

    return false;
  }  
  return true;
}

bool QueryImpl::ExpectNotEnd(const QueryToken& token, const String& s) {
  if (token.IsEOF()) {
    error_desc_.msg="expected "+s+", but query string ended "
                    "unexpectedly";
    error_desc_.range=Range(query_string_.length(),0);    
    return false;
  } else if (token.GetType()==tok::UnterminatedQuote) {
    error_desc_.msg="unterminated quote in query string";
    error_desc_.range=token.GetRange();    
    return false;
  }
  return true;
}

bool QueryImpl::ExpectComparisonOperator(const QueryToken& token) {
  if (!(token.GetType()==tok::NotEqual || token.GetType()==tok::Equal || 
      token.GetType()==tok::GreaterEqualThan || 
      token.GetType()==tok::LessEqualThan ||
      token.GetType()==tok::LessThan || token.GetType()==tok::GreaterThan)) {
      if (token.IsEOF()) {
        error_desc_.msg="expected comparison operator, but query string ended "
                        "unexpectedly";
        error_desc_.range=Range(query_string_.length(),0);
      } else if (token.GetType()==tok::UnterminatedQuote) {
        error_desc_.msg="unterminated quote in query string";
        error_desc_.range=token.GetRange();    
      } else {
        error_desc_.msg="expected comparison operator, but '"+
                        query_string_.substr(token.GetRange().Loc, 
                                             token.GetRange().Length)+"' found";                                             
        error_desc_.range=token.GetRange();
      }
      return false;
  }
  return true;
}

bool QueryImpl::Expect(tok::Type type, const String& s,
                       const QueryToken& token) {
  if (type==token.GetType())
    return true;
  if (token.IsEOF()) {
    error_desc_.msg="expected "+s+", but query string ended "
                    "unexpectedly";
    error_desc_.range=Range(query_string_.length(), 0);    
  }   else if (token.GetType()==tok::UnterminatedQuote) {
      error_desc_.msg="unterminated quote in query string";
      error_desc_.range=token.GetRange();    
  } else {
    error_desc_.msg="expected "+s+", but '"+
      query_string_.substr(token.GetRange().Loc,
                           token.GetRange().Length)+"' found";
    error_desc_.range=token.GetRange();    
  }

  return false;
}

QueryState QueryImpl::CreateQueryState(const EntityHandle& h) const 
{
  if (has_error_) {
    throw QueryError(query_string_, error_desc_);
  }
  return QueryState(*this, h);
}

QueryState QueryImpl::CreateQueryState(const EntityView& v) const 
{
  if (has_error_) {
    throw QueryError(query_string_, error_desc_);
  }
  return QueryState(*this, v);
}

bool QueryImpl::ExpectNumeric(const QueryToken& token) 
{
  tok::Type t=token.GetType();
  if (tok::FloatingValue==t || tok::IntegralValue==t)
    return true;
  if (token.IsEOF()) {
    error_desc_.msg="expected numeric value, but query string ended "
                    "unexpectedly";
    error_desc_.range=Range(query_string_.length(), 0);    
  } else if (token.GetType()==tok::UnterminatedQuote) {
        error_desc_.msg="unterminated quote in query string";
        error_desc_.range=token.GetRange();    
    } else {
    error_desc_.msg="expected numeric value, but '"+
      query_string_.substr(token.GetRange().Loc,
                           token.GetRange().Length)+"' found";
    error_desc_.range=token.GetRange();    
  }
  return false;
}


Node* QueryImpl::Concatenate(Node* lhs, Node* rhs, LogicOP logical_op) {
  assert(lhs && "lhs is NULL");
  assert(rhs && "rhs is NULL");
  LogicOPNode* lop_node=new LogicOPNode(logical_op);
  lop_node->SetLHS(lhs);
  lop_node->SetRHS(rhs);
  return lop_node;
}

bool QueryImpl::ParsePoint(QueryLexer& lexer, geom::Vec3& point) {
  QueryToken f1=lexer.NextToken();
  if (!this->ExpectNumeric(f1))
    return false;
  if (!this->Expect(tok::Coma, "','", lexer.NextToken())) 
    return false;
  QueryToken f2=lexer.NextToken();  
  if (!this->ExpectNumeric(f2))
    return false; 
  if (!this->Expect(tok::Coma, "','", lexer.NextToken())) 
    return false;
  QueryToken f3=lexer.NextToken();  
  if (!this->ExpectNumeric(f3))
    return false; 
  if (!this->Expect(tok::RightCurlyBrace, "'}'", lexer.NextToken()))
    return false;

  float fv1=atof(query_string_.substr(f1.GetRange().Loc,
                                      f1.GetRange().Length).c_str());
  float fv2=atof(query_string_.substr(f2.GetRange().Loc,
                                      f2.GetRange().Length).c_str());
  float fv3=atof(query_string_.substr(f3.GetRange().Loc,
                                      f3.GetRange().Length).c_str());                                                                            
  point=geom::Vec3(fv1, fv2, fv3);
  return true;                                   
}

Node* QueryImpl::ParseWithinExpr(QueryLexer& lexer) {
  QueryToken radius=lexer.CurrentToken();
  if (!this->ExpectNumeric(radius))
    return NULL;  
  float rv=atof(query_string_.substr(radius.GetRange().Loc,
                                     radius.GetRange().Length).c_str());    
  QueryToken within_op=lexer.NextToken();
  if (!this->Expect(tok::Within, "'<>'", within_op))
    return NULL;
  if (!this->ExpectNotEnd(lexer.NextToken(), 
                          "right hand operand of within statement"))
    return NULL;
  QueryToken ct=lexer.CurrentToken();
  if (ct.GetType()==tok::LeftCurlyBrace) {
    Node* points=NULL;
    while (true) {
      geom::Vec3 point;
      if (this->ParsePoint(lexer, point)) {
        ParamType pt(WithinParam(point, rv*rv));
        CompOP comp_op=COP_LE;
        if (inversion_stack_.back())
          comp_op=COP_GE;
        SelNode* within_node=new SelNode(Prop(Prop::WITHIN, Prop::VEC_DIST,
                                              Prop::ATOM), 
                                         comp_op, pt);
        if (points) {
          points=this->Concatenate(points, within_node, LOP_OR);
        } else {
          points=within_node;
        }
        QueryToken nt=lexer.NextToken();
        if (nt.GetType()==tok::Coma) {
          nt=lexer.NextToken();
          if (!this->Expect(tok::LeftCurlyBrace, "'{'", nt)) {
            delete points;
            return NULL;
          }
          continue;
        }
        return points;
      } else {
        return NULL;
      }
    }
  } else if (ct.GetType()==tok::LeftBracket) {
    // push false onto inversion stack to make sure we have a proper start for
    // the bracketed sub expression.
    inversion_stack_.push_back(false);
    Node* rhs=this->ParseBracketSubExpr(lexer);
    if (!rhs) {
      return NULL;      
    }
    // build new query from parsed AST. inside the WithinParam we store the 
    // index of the query for later access.
    QueryImplP q(new QueryImpl(rhs));
    bracketed_expr_.push_back(q);
    q->gen_prop_list_ = gen_prop_list_;
    ParamType pt(WithinParam(bracketed_expr_.size()-1, rv*rv));    
    inversion_stack_.pop_back();    
    CompOP comp_op= COP_LE;
    if (inversion_stack_.back()) {
      comp_op=COP_GE;
    }
      
    SelNode* within_node=new SelNode(Prop(Prop::WITHIN, Prop::VEC_DIST,
                                          Prop::ATOM), 
                                     comp_op, pt);
    return within_node;
  } else {
    error_desc_.msg="expected right hand operand of within statement, but '"+
                    query_string_.substr(ct.GetRange().Loc,
                                         ct.GetRange().Length)+"' found";
    error_desc_.range=ct.GetRange();
    return NULL;
  }
}

Node* QueryImpl::BuildAST() {
  QueryLexer lexer(query_string_);
  QueryToken t;
  t=lexer.NextToken();
  if (t.IsEOF()) {
    return NULL;
  }
  inversion_stack_.push_back(false);
  Node* node=this->ParseSubExpr(lexer);
  if (node) {
    assert(inversion_stack_.size()==1);
  }
  return node;
}

}}} // ns
