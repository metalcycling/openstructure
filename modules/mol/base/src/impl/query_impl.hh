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
#ifndef OST_QUERY_IMPL_HH
#define OST_QUERY_IMPL_HH

#include <vector>
#include <string>
#include <set>

#include <boost/logic/tribool.hpp>

#include <ost/mol/module_config.hh>
#include "query_impl_fw.hh"
#include "query_ast.hh"
#include <ost/range.hh>
#include <ost/mol/query_error.hh>
#include "atom_impl_fw.hh"
#include "residue_impl_fw.hh"
#include "chain_impl_fw.hh"
#include <ost/mol/mol.hh>
#include <ost/mol/query_state.hh>
#include <ost/mol/entity_property_mapper.hh>

namespace ost { namespace mol { namespace impl {

namespace tok {
  /// \internal
  typedef enum {
    IntegralValue,
    FloatingValue,
    Colon,
    Coma,
    Dot, 
    And,
    Or, 
    Not,
    Equal,
    NotEqual,
    GreaterEqualThan,
    GreaterThan,
    LessThan,
    LessEqualThan,
    Identifier,
    Within,
    LeftCurlyBrace,
    RightCurlyBrace,
    LeftParen,
    RightParen,
    LeftBracket,
    RightBracket,
    None,
    UnknownChar,
    UnterminatedQuote,
    String
  }
  Type;
}

/// \internal
class QueryToken {
public:
  QueryToken();
  QueryToken(const Range& range, tok::Type type, 
             const Range& value_range=Range());
  
  tok::Type GetType() const;
  //! Returns true if token is end of file (or end of line) token.
  bool IsEOF() const;
  //! Get character range of this token.
  const Range& GetRange() const;
  /// Get value range of this token. This is usually identical to GetRange()
  /// except when dealing with quoted strings.
  const Range& GetValueRange() const;
private:
  tok::Type      type_;
  Range          range_;
  Range          value_range_;
};

/// \internal
class QueryLexer {
public:
  QueryLexer(const String& query_string);
  QueryToken NextToken();
  QueryToken CurrentToken() const;
private:  
  //! Reads over whitespace characters directly on the righ side of the current 
  //  lexer position.
  void EatWhitespaces();
  //! Lex next token
  QueryToken LexToken();
  QueryToken LexQuotedStringLiteral();
  QueryToken LexNumericToken();
  QueryToken LexIdentOrStringToken();
  
  String   query_string_;
  size_t   current_;  
  QueryToken tok_;
};

/// \internal
struct SelStmt {
  Prop::ID                   sel_id;
  CompOP                       comp_op;
  ParamType                    param;    
  SelStmt()
    : sel_id(Prop::UNDEF),
      comp_op(COP_EQ),
      param(0) { }
  SelStmt(Prop::ID id ,CompOP op, const ParamType& parm) 
    : sel_id(id),
      comp_op(op),
      param(parm) { }
      
  bool operator ==(const SelStmt& rhs) {
    return sel_id == rhs.sel_id && comp_op == rhs.comp_op && param == rhs.param;
  }
};
typedef enum {
  LOGIC, VALUE
} SelItemType;

/// \internal
struct SelItem {
  SelItemType  type;
  size_t       value;
};

/// \brief Struct to store info about GenericProperties parsed from the
/// \internal
///        query string
struct GenProp {
  EntityPropertyMapper mapper;
  float                default_val;
  bool                 has_default;
  GenProp(EntityPropertyMapper epm) :
         mapper(epm), default_val(0.0), has_default(false) {}
};

/// \internal
typedef std::vector<SelItem> SelStack;

/// \internal
typedef std::vector<SelStmt> SelStmts;
  
/// \brief query statement implementation 
/// \internal
class QueryImpl {
public:
  friend class mol::QueryState;
  /// Construct new QueryImpl object.
  QueryImpl(const String& query_string="");
  /// Construct new QueryImpl from AST tree. Used internally to build 
  /// representation of bracketed sub expressions.
  QueryImpl(Node* ast);
  /// Find out if the query String is valid
  bool IsValid() const;
  
  /// Check if query was empty, ie matching all
  bool IsEmpty() const;

  /// Retrieve query String
  const String& GetQueryString() const;
  

  const QueryErrorDesc& GetErrorDescription() const;
  
  mol::QueryState CreateQueryState(const EntityHandle& h) const;
  mol::QueryState CreateQueryState(const EntityView& v) const;
private:
  
  /// Builds an abstract syntax tree from the query String.
  Node* BuildAST();  
  /// Transform AST into dual stack representation
  void ASTToSelStack(const Node* src_ast,Prop::Level target_level,
                     SelStack& stack);
  
  /// Determine whether the AST starting at ast_node has a defined result when 
  /// evaluated at target_level.
  bool IsAlwaysUndef(const Node* ast_node, Prop::Level target_level);
  
  /// extract all sel statements from the AST
  void ExtractSelStmts(const Node* ast);
  
  /// Parse a single predicate
  Node* ParsePropValueExpr(QueryLexer& lexer);
  
  /// Parse sub-expression. Will return the parsed sub-expression or NULL on
  /// failure.
  Node* ParseSubExpr(QueryLexer& lexer, bool paren=false);

  /// Parse sub-expression that starts with a '('. Will return the parsed 
  /// sub-expression on success and NULL on failure.
  Node* ParseParenSubExpr(QueryLexer& lexer);  
  
  /// Parse sub-expression that starts with a '['. This may occur on the right
  /// hand-side of a within statement. Will return the parsed sub-expression 
  /// on success and NULL on failure.  
  Node* ParseBracketSubExpr(QueryLexer& lexer);
  
  /// Parse within statement. Will return the property or NULL if the statement
  ///  was not correct.
  Node* ParseWithinExpr(QueryLexer& lexer);
  
  /// \brief Parse point.
  /// 
  /// Current token must point to the opening '{'. On exit, the current token 
  /// will point to the closing '}'.
  /// 
  /// \param lexer query lexer
  /// \param point upon success will contain the coordinates
  /// \return true on success false when a parsing error occurs.
  bool ParsePoint(QueryLexer& lexer, geom::Vec3& point);
  
  /// Check if token has the wanted type. Will set error_desc_ appropriately
  ///  and return false if this is not the case.
  bool Expect(tok::Type type, const String& s, const QueryToken& token);
  
  /// Check if token has numeric type (int or float). Will set error_desc_ 
  ///  appropriately and return false if this is not the case.  
  bool ExpectNumeric(const QueryToken& token);
  
  /// Check if token is a logical operator ('or' and 'and'). Will set 
  ///  error_desc_ appropriately and return false if this is not the case.
  bool ExpectLogicalOperator(const QueryToken& token);
  
  /// Check if token is a comparison operator. Will return false and set
  ///  error_desc_ appropriately if this is not the case.
  bool ExpectComparisonOperator(const QueryToken& token);
  
  /// Check if the end of the query String has not been reached yet. 
  /// 
  /// \param s is the error that should be emitted in case the end of the
  /// \param token the current token
  /// 
  /// query String is reached.
  bool ExpectNotEnd(const QueryToken& token, const String& s);
  
  /// Parse value. On success, this method will return true and value will be 
  /// set to the parsed value.
  bool ParseValue(const Prop& sel, const QueryToken& op,
                  QueryLexer& lexer, ParamType& value);

  /// Parse value or value range. On success the method will return a SelNode*
  ///  or LogicOPNode*. If the method fails, NULL will be returned.
  Node* ParseValueOrRange(const Prop& s, const QueryToken& op, 
                          QueryLexer& lexer);
  /// Concatenate two nodes together using LogicOP
  Node* Concatenate(Node* lhs, Node* rhs, LogicOP logical_op);
  
  String                            query_string_;
  bool                              has_error_;
  bool                              empty_optimize_;
  QueryErrorDesc                    error_desc_;
  std::vector<GenProp>              gen_prop_list_;
  int                               num_gen_prop_;
protected:  
  std::vector<SelStmt>              sel_values_;
  std::vector<SelStack>             sel_stacks_;
  std::vector<std::set<size_t> >    indices_;
  std::vector<bool>                 inversion_stack_;
  std::vector<QueryImplP>           bracketed_expr_;
};
 
}}} // ns

#endif
