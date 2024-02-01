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
 *  Authors: Marco Biasini, Juergen Haas
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/mol/query.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/entity_view.hh>

using namespace ost;
using namespace ost::mol;

EntityHandle make_query_test_entity()
{
  EntityHandle eh = CreateEntity();
  XCSEditor e=eh.EditXCS();
  ChainHandle chain = e.InsertChain("A");
  chain.SetFloatProp("testprop_c", 1.0);
  ResidueHandle res = e.AppendResidue(chain, "MET");
  res.SetFloatProp("testpropr", 1.0);
  AtomHandle ah=e.InsertAtom(res, "N",geom::Vec3(21.609,35.384,56.705), "N");
  ah.SetFloatProp("testpropa", 1.0);
  e.InsertAtom(res, "CA",geom::Vec3(20.601,35.494,57.793), "C");
  e.InsertAtom(res, "C",geom::Vec3(19.654,34.300,57.789), "C");
  e.InsertAtom(res, "O",geom::Vec3(18.447,34.456,57.595), "O");
  e.InsertAtom(res, "CB",geom::Vec3(19.789,36.783,57.639), "C");
  e.InsertAtom(res, "CG",geom::Vec3(20.629,38.055,57.606), "C");
  e.InsertAtom(res, "SD",geom::Vec3(21.638,38.325,59.084), "S");
  e.InsertAtom(res, "CE",geom::Vec3(23.233,37.697,58.529), "C");

  res = e.AppendResidue(chain, "ARG");
  e.InsertAtom(res, "N",geom::Vec3(20.202,33.112,58.011), "N");
  e.InsertAtom(res, "CA",geom::Vec3(19.396,31.903,58.033), "C");
  e.InsertAtom(res, "C",geom::Vec3(18.608,31.739,59.328), "C");
  e.InsertAtom(res, "O",geom::Vec3(17.651,30.965,59.381), "O");
  e.InsertAtom(res, "CB",geom::Vec3(20.284,30.681,57.801), "C");
  e.InsertAtom(res, "CG",geom::Vec3(20.665,30.488,56.342), "C");
  e.InsertAtom(res, "CD",geom::Vec3(21.557,29.281,56.154), "C");
  e.InsertAtom(res, "NE",geom::Vec3(22.931,29.557,56.551), "N");
  e.InsertAtom(res, "CZ",geom::Vec3(23.901,28.653,56.528), "C");
  e.InsertAtom(res, "NH1",geom::Vec3(23.640,27.417,56.130), "N");
  e.InsertAtom(res, "NH2",geom::Vec3(25.132,28.980,56.893), "N");

  res = e.AppendResidue(chain, "LEU");
  e.InsertAtom(res, "N",geom::Vec3(19.003,32.473,60.366), "N");
  e.InsertAtom(res, "CA",geom::Vec3(18.330,32.402,61.664), "C");
  e.InsertAtom(res, "C",geom::Vec3(17.884,33.787,62.117), "C");
  e.InsertAtom(res, "O",geom::Vec3(17.853,34.091,63.308), "O");
  e.InsertAtom(res, "CB",geom::Vec3(19.269,31.793,62.710), "C");
  e.InsertAtom(res, "CG",geom::Vec3(19.695,30.340,62.501), "C");
  e.InsertAtom(res, "CD1",geom::Vec3(20.585,29.897,63.648), "C");
  e.InsertAtom(res, "CD2",geom::Vec3(18.461,29.459,62.420), "C");
  return eh;
}

void ensure_counts(EntityHandle e, const String& qs, int cc, int rc, int ac) {
  EntityView v;
  BOOST_CHECK_NO_THROW(v=e.Select(qs));
  BOOST_CHECK_MESSAGE(v.GetChainCount()==cc,
                      "wrong chain count " << v.GetChainCount()
                      << " for query string " << qs);
  BOOST_CHECK_MESSAGE(v.GetResidueCount()==rc,
                      "wrong residue count " << v.GetResidueCount() <<
                      " for query string " << qs);
  BOOST_CHECK_MESSAGE(v.GetAtomCount()==ac,
                      "wrong atom count " << v.GetAtomCount() <<
                      " for query string " << qs);
}


void ensure_counts_v(EntityView src, const String& qs,
                     int cc, int rc, int ac) {
  EntityView v;
  BOOST_CHECK_NO_THROW(v=src.Select(qs));
  BOOST_CHECK_MESSAGE(v.GetChainCount()==cc,
                      "wrong chain count " << v.GetChainCount()
                      << " for query String " << qs);
  BOOST_CHECK_MESSAGE(v.GetResidueCount()==rc,
                      "wrong residue count " << v.GetResidueCount() <<
                      " for query String " << qs);
  BOOST_CHECK_MESSAGE(v.GetAtomCount()==ac,
                      "wrong atom count " << v.GetAtomCount() <<
                      " for query String " << qs);
}

BOOST_AUTO_TEST_SUITE( mol_base );

BOOST_AUTO_TEST_CASE(test_query_parse_properties) 
{
  BOOST_CHECK(Query("rname=Ala").IsValid());
  BOOST_CHECK(Query("ele=C").IsValid());
  BOOST_CHECK(Query("aname=X").IsValid());
  BOOST_CHECK(Query("rnum=3").IsValid());
  BOOST_CHECK(Query("cname=A").IsValid());
  BOOST_CHECK(Query("occ=50").IsValid());
  BOOST_CHECK(Query("occ=50").IsValid());
  BOOST_CHECK(Query("rtype=H").IsValid());
  BOOST_CHECK(Query("x=3").IsValid());
  BOOST_CHECK(Query("y=4").IsValid());
  BOOST_CHECK(Query("z=6").IsValid());
  BOOST_CHECK(Query("aindex=1").IsValid());
  BOOST_CHECK(Query("gatest=7").IsValid());
  BOOST_CHECK(Query("grtest=8").IsValid());
  BOOST_CHECK(Query("gctest=9").IsValid());
  BOOST_CHECK(Query("gatest:1=7").IsValid());
  BOOST_CHECK(Query("grtest:2=8").IsValid());
  BOOST_CHECK(Query("gctest:3.0=9").IsValid());
  BOOST_CHECK(Query("anita=3").IsValid()==false);
  BOOST_CHECK(Query("gc*test=3").IsValid()==false);
  BOOST_CHECK(Query("gc?test=3").IsValid()==false);
}

BOOST_AUTO_TEST_CASE(test_query_parse_value_type) 
{
  BOOST_CHECK(Query("rnum<=3").IsValid());
  BOOST_CHECK(Query("rnum>=3").IsValid());
  BOOST_CHECK(Query("rnum<3").IsValid());
  BOOST_CHECK(Query("rnum>3").IsValid());
  BOOST_CHECK(Query("rnum!=3").IsValid());

  BOOST_CHECK(Query("rnum=3:10").IsValid());
  BOOST_CHECK(Query("rnum!=3:10").IsValid());
  BOOST_CHECK(Query("rnum!=3:10").IsValid());
  BOOST_CHECK(Query("rnum>=30:40").IsValid()==false);
  BOOST_CHECK(Query("rnum>30:40").IsValid()==false);
  BOOST_CHECK(Query("rnum<=30:40").IsValid()==false);
  BOOST_CHECK(Query("rnum<30:40").IsValid()==false);

  BOOST_CHECK(Query("rnum=1,2").IsValid());
  BOOST_CHECK(Query("rnum=1,2,3,4").IsValid());
  BOOST_CHECK(Query("rnum>1,2,3,4").IsValid()==false);

  BOOST_CHECK(Query("rnum=WTF").IsValid()==false);
  BOOST_CHECK(Query("rnum=WTF").IsValid()==false);
  BOOST_CHECK(Query("rnum=3.0").IsValid()==false);
  BOOST_CHECK(Query("ele>=XXX").IsValid()==false);

  BOOST_CHECK(Query("aindex=1,2").IsValid());
  BOOST_CHECK(Query("aindex=1:10,12:20").IsValid());
  BOOST_CHECK(Query("aindex>7").IsValid());
}

BOOST_AUTO_TEST_CASE(test_query_parse_logical_op) 
{
  BOOST_CHECK(Query("rnum=3 and rnum=5").IsValid());
  BOOST_CHECK(Query("rnum=3 or rnum=5").IsValid());
  BOOST_CHECK(Query("not rnum=3").IsValid());
  BOOST_CHECK(Query("not rnum=3").IsValid());
  BOOST_CHECK(Query("not rnum=3 and not rnum=5").IsValid());
  BOOST_CHECK(Query("not not rnum=3").IsValid()==false);
}

BOOST_AUTO_TEST_CASE(test_query_parse_unexpected) 
{
  BOOST_CHECK(Query("cname").IsValid()==false);
  BOOST_CHECK(Query("cname=").IsValid()==false);
  BOOST_CHECK(Query("not").IsValid()==false);
  BOOST_CHECK(Query("(rnum=3").IsValid()==false);
  BOOST_CHECK(Query("(rnum=3").IsValid()==false);
  BOOST_CHECK(Query("(rnum=3").IsValid()==false);
  BOOST_CHECK(Query("(not").IsValid()==false);
  BOOST_CHECK(Query("()").IsValid()==false);
  BOOST_CHECK(Query("(rnum)").IsValid()==false);
  BOOST_CHECK(Query("(rnum=)").IsValid()==false);
  BOOST_CHECK(Query("gatest:=5").IsValid()==false);
  BOOST_CHECK(Query("gatest:5.0=").IsValid()==false);
  BOOST_CHECK(Query("gatest:=").IsValid()==false);
}

BOOST_AUTO_TEST_CASE(test_query_parse_within) 
{
  BOOST_CHECK(Query("5 <> {0,0,0.0}").IsValid());
  BOOST_CHECK(Query("5 <> {0}").IsValid()==false);
  BOOST_CHECK(Query("5 <> {0").IsValid()==false);
  BOOST_CHECK(Query("5 <> {0,}").IsValid()==false);
  BOOST_CHECK(Query("5 <> {0,0.}").IsValid()==false);
}

BOOST_AUTO_TEST_CASE(test_query_eval) 
{
  EntityHandle e=make_query_test_entity();

  ensure_counts(e, "", 1, 3, 27);
  ensure_counts(e, "cname=B", 0, 0, 0);
  ensure_counts(e, "cname=A", 1, 3, 27);
  ensure_counts(e, "aname=CA", 1, 3, 3);
  ensure_counts(e, "aname=SD", 1, 1, 1);
  ensure_counts(e, "aindex=1,3,99", 1, 1, 2);
  ensure_counts(e, "rnum=1:2", 1, 2, 19);
  ensure_counts(e, "rnum=1,2", 1, 2, 19);
  ensure_counts(e, "rnum>3", 0, 0, 0);
  ensure_counts(e, "rnum<1", 0, 0, 0);
  ensure_counts(e, "ele=S", 1, 1, 1);
  ensure_counts(e, "rname=GLY", 0, 0, 0);
  ensure_counts(e, "not rname=GLY", 1, 3, 27);
  ensure_counts(e, "not rname!=GLY", 0, 0, 0);
  ensure_counts(e, "aname=CA and not aname=CA", 0, 0, 0);
  ensure_counts(e, "rtype=H", 0, 0, 0);
  ensure_counts(e, "rtype=C", 1, 3, 27);
  ensure_counts(e, "not (aname=CA and not aname=CA)", 1, 3, 27);
  ensure_counts(e, "3 <> {21.5,35,57.0}", 1, 2, 5);
  ensure_counts(e, "not 3 <> {21.5,35,57.0}", 1, 3, 22);
  ensure_counts(e, "3 <> {21.5,35,57} and not 0.5 <> {21.5,35,57} ", 1, 2, 4);
  ensure_counts(e, "not 0.5 <> [rnum=3]", 1, 2, 19);
  ensure_counts(e, "1 <> {0,0,0}", 0, 0, 0);
  ensure_counts(e, "gatestpropa:0=1", 1, 1, 1);
  ensure_counts(e, "gatestpropa:1.0=1", 1, 3, 27);
  ensure_counts(e, "gatestpropa:2.0=2", 1, 3, 26);
  ensure_counts(e, "grtestpropr:0=1", 1, 1, 8);
  ensure_counts(e, "grtestpropr:1.0=1", 1, 3, 27);
  ensure_counts(e, "grtestpropr:2.0=2", 1, 2, 19);
  ensure_counts(e, "gctestprop_c:0=1", 1, 3, 27);
  ensure_counts(e, "gctestprop_c:1.0=1", 1, 3, 27);
  ensure_counts(e, "gctestprop_c:2.0=2", 0, 0, 0);
  ensure_counts(e, "ganotthere:false=true", 0, 0, 0);
  ensure_counts(e, "ganotthere:false=false", 1, 3, 27);
}

BOOST_AUTO_TEST_CASE(test_query_eval_on_view) 
{
  EntityHandle e=make_query_test_entity();
  EntityView v=e.CreateFullView();
  ensure_counts_v(v, "", 1, 3, 27);
  ensure_counts_v(v, "cname=B", 0, 0, 0);
  ensure_counts_v(v, "cname=A", 1, 3, 27);
  ensure_counts_v(v, "aname=CA", 1, 3, 3);
  ensure_counts_v(v, "aname=SD", 1, 1, 1);
  ensure_counts_v(v, "ele=S", 1, 1, 1);
  ensure_counts_v(v, "rname=GLY", 0, 0, 0);
  ensure_counts_v(v, "not rname=GLY", 1, 3, 27);
  ensure_counts_v(v, "not rname!=GLY", 0, 0, 0);
  ensure_counts_v(v, "aname=CA and not aname=CA", 0, 0, 0);
  ensure_counts_v(v, "rtype=H", 0, 0, 0);
  ensure_counts_v(v, "rtype=C", 1, 3, 27);
  ensure_counts_v(v, "not (aname=CA and not aname=CA)", 1, 3, 27);
  ensure_counts_v(v, "rnum=1:2", 1, 2, 19);
  ensure_counts_v(v, "rnum=1,2", 1, 2, 19);
  ensure_counts_v(v, "rnum>3", 0, 0, 0);
  ensure_counts_v(v, "rnum<1", 0, 0, 0);
  ensure_counts_v(v, "gatestpropa:0=1", 1, 1, 1);
  ensure_counts_v(v, "gatestpropa:1.0=1", 1, 3, 27);
  ensure_counts_v(v, "gatestpropa:2.0=2", 1, 3, 26);
  ensure_counts_v(v, "grtestpropr:0=1", 1, 1, 8);
  ensure_counts_v(v, "grtestpropr:1.0=1", 1, 3, 27);
  ensure_counts_v(v, "grtestpropr:2.0=2", 1, 2, 19);
  ensure_counts_v(v, "gctestprop_c:0=1", 1, 3, 27);
  ensure_counts_v(v, "gctestprop_c:1.0=1", 1, 3, 27);
  ensure_counts_v(v, "gctestprop_c:2.0=2", 0, 0, 0);
}

BOOST_AUTO_TEST_CASE(test_query_throw) 
{
  EntityHandle e=make_query_test_entity();
  BOOST_CHECK_NO_THROW(e.Select("gatestpropa:0=1"));
  BOOST_CHECK_NO_THROW(e.Select("gatestpropa:1=1"));
  BOOST_CHECK_THROW(e.Select("gatestpropa=1"), ost::Error);
  BOOST_CHECK_NO_THROW(e.Select("grtestpropr:0=1"));
  BOOST_CHECK_NO_THROW(e.Select("grtestpropr:1=1"));
  BOOST_CHECK_THROW(e.Select("grtestpropr=1"), ost::Error);
  BOOST_CHECK_NO_THROW(e.Select("gctestprop_c:0=1"));
  BOOST_CHECK_NO_THROW(e.Select("gctestprop_c:1=1"));
  BOOST_CHECK_THROW(e.Select("ganotsetprop=1"), ost::Error);
  BOOST_CHECK_THROW(e.Select("grnotsetprop=1"), ost::Error);
  BOOST_CHECK_THROW(e.Select("gcnotsetprop=1"), ost::Error);
  BOOST_CHECK_NO_THROW(e.Select("ganotsetprop:0=1"));
  BOOST_CHECK_NO_THROW(e.Select("grnotsetprop:0=1"));
  BOOST_CHECK_NO_THROW(e.Select("gcnotsetprop:0=1"));
}

BOOST_AUTO_TEST_CASE(test_glob)
{
  EntityHandle e=make_query_test_entity();
  ensure_counts(e, "rname=MET and aname=C*", 1, 1, 5);
  ensure_counts(e, "rname=ARG and aname=N?1", 1, 1, 1);
  ensure_counts(e, "rname=ARG and aname=NH?", 1, 1, 2);
  ensure_counts(e, "rname=ARG and aname=\"*2\"", 1, 1, 1);
  ensure_counts(e, "rname=ARG and aname=*2", 1, 1, 1);
  ensure_counts(e, "rname=ARG and aname=N?", 1, 1, 1);
  ensure_counts(e, "rname=LEU and aname=\"?D?\"", 1, 1, 2);
  ensure_counts(e, "rname=LEU and aname=?D?", 1, 1, 2);
}

BOOST_AUTO_TEST_CASE(test_quoting)
{
  // possible letters taken from mmCIF dictionary
  // -> http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_atom_site.label_asym_id.html
  // not ok for us: '\' alone
  String single_letters = " ][_,.;:\"&<>()/{}'`~!@#$\%|+-"
                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefgh"
                          "ijklmnopqrstuvwxyz";
  // put into vector
  std::vector<String> chain_names;
  for (size_t i = 0; i < single_letters.length(); ++i) {
    chain_names.push_back(String(1, single_letters[i]));
  }
  // add some multi letter names and empty chain name too
  // not ok: mixing ' and ", '\' at end
  chain_names.push_back("\\ ][_,.;:&<>()/{}'`~!@#$\%|+-");
  chain_names.push_back("ABCDEFGHIJKLMNOPQRSTU'VWXYZ0123456789abcdefgh");
  chain_names.push_back("ijklmno\"pqrstuvwxyz");
  chain_names.push_back("");
  // setup entity
  EntityHandle ent = CreateEntity();
  XCSEditor edi = ent.EditXCS();
  for (size_t i = 0; i < chain_names.size(); ++i) {
    edi.InsertChain(chain_names[i]);
  }
  // test quoting
  String query_all = "cname=";
  for (size_t i = 0; i < chain_names.size(); ++i) {
    const String quoted_name = QueryQuoteName(chain_names[i]);
    BOOST_CHECK_EQUAL(ent.Select("cname=" + quoted_name).GetChainCount(), 1);
    query_all += quoted_name;
    if (i != chain_names.size() - 1) query_all += ",";
  }
  BOOST_CHECK_EQUAL(ent.Select(query_all).GetChainCount(),
                    int(chain_names.size()));
  // note: quoting * keeps it as wild card!
  BOOST_CHECK_EQUAL(ent.Select("cname=" + QueryQuoteName("*")).GetChainCount(),
                    int(chain_names.size()));
  // expected failures (mixing ' and ", '\' at end)
  BOOST_CHECK_THROW(QueryQuoteName("'\""), Error);
  BOOST_CHECK_THROW(QueryQuoteName("\\"), Error);
  BOOST_CHECK_THROW(QueryQuoteName("a\\"), Error);
}

BOOST_AUTO_TEST_SUITE_END();
