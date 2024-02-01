#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/mol/xcs_editor.hh>
#include <ost/mol/bond_handle.hh>

#include "helper.hh"

using namespace ost::mol;

namespace ost { namespace conop {

ResidueHandle make_cytosine(ChainHandle chain) 
{
  XCSEditor e=chain.GetEntity().EditXCS();  
  ResidueHandle res=e.AppendResidue(chain, "DC");

  e.InsertAtom(res, "P",   geom::Vec3(21.412, 34.391, 37.142));
  e.InsertAtom(res, "OP1", geom::Vec3(22.938, 34.599, 36.988));
  e.InsertAtom(res, "OP2", geom::Vec3(20.690, 35.640, 37.689));
  e.InsertAtom(res, "O5'", geom::Vec3(21.215, 33.299, 38.339));
  e.InsertAtom(res, "C5'", geom::Vec3(20.524, 33.660, 39.548));
  e.InsertAtom(res, "C4'", geom::Vec3(19.064, 33.285, 39.452));
  e.InsertAtom(res, "O4'", geom::Vec3(18.246, 34.411, 39.034));
  e.InsertAtom(res, "C3'", geom::Vec3(18.778, 32.150, 38.469));
  e.InsertAtom(res, "O3'", geom::Vec3(17.930, 31.199, 39.081));
  e.InsertAtom(res, "C2'", geom::Vec3(18.046, 32.827, 37.330));
  e.InsertAtom(res, "C1'", geom::Vec3(17.326, 33.955, 38.048));
  e.InsertAtom(res, "N1",  geom::Vec3(16.891, 35.100, 37.196));
  e.InsertAtom(res, "C2",  geom::Vec3(15.709, 35.782, 37.543));
  e.InsertAtom(res, "O2",  geom::Vec3(15.137, 35.489, 38.610));
  e.InsertAtom(res, "N3",  geom::Vec3(15.226, 36.742, 36.710));
  e.InsertAtom(res, "C4",  geom::Vec3(15.878, 37.041, 35.582));
  e.InsertAtom(res, "N4",  geom::Vec3(15.324, 37.943, 34.762));
  e.InsertAtom(res, "C5",  geom::Vec3(17.116, 36.415, 35.238));
  e.InsertAtom(res, "C6",  geom::Vec3(17.584, 35.464, 36.067));

  return res;
}

ResidueHandle make_uracil1(ChainHandle chain) 
{
  XCSEditor e=chain.GetEntity().EditXCS();
  ResidueHandle res = e.AppendResidue(chain, "U");

  e.InsertAtom(res, "P",   geom::Vec3(18.533, 30.238, 40.226));
  e.InsertAtom(res, "OP1", geom::Vec3(20.012, 30.098, 40.003));
  e.InsertAtom(res, "OP2", geom::Vec3(17.680, 29.017, 40.258));
  e.InsertAtom(res, "O5'", geom::Vec3(18.237, 31.056, 41.570));
  e.InsertAtom(res, "C5'", geom::Vec3(17.104, 31.923, 41.642));
  e.InsertAtom(res, "C4'", geom::Vec3(16.332, 31.729, 42.937));
  e.InsertAtom(res, "O4'", geom::Vec3(15.180, 32.616, 42.841));
  e.InsertAtom(res, "C3'", geom::Vec3(15.785, 30.319, 43.194));
  e.InsertAtom(res, "O3'", geom::Vec3(15.960, 29.834, 44.527));
  e.InsertAtom(res, "C2'", geom::Vec3(14.311, 30.371, 42.766));
  e.InsertAtom(res, "O2'", geom::Vec3(13.442, 29.708, 43.677));
  e.InsertAtom(res, "C1'", geom::Vec3(13.972, 31.863, 42.765));
  e.InsertAtom(res, "N1",  geom::Vec3(13.282, 32.212, 41.499));
  e.InsertAtom(res, "C2",  geom::Vec3(12.072, 32.896, 41.568));
  e.InsertAtom(res, "O2",  geom::Vec3(11.510, 33.142, 42.626));
  e.InsertAtom(res, "N3",  geom::Vec3(11.535, 33.257, 40.349));
  e.InsertAtom(res, "C4",  geom::Vec3(12.048, 32.976, 39.096));
  e.InsertAtom(res, "O4",  geom::Vec3(11.490, 33.447, 38.087));
  e.InsertAtom(res, "C5",  geom::Vec3(13.268, 32.207, 39.106));
  e.InsertAtom(res, "C6",  geom::Vec3(13.831, 31.872, 40.261));

  return res;
}

ResidueHandle make_uracil2(ChainHandle chain) 
{
  XCSEditor e=chain.GetEntity().EditXCS();
  ResidueHandle res = e.AppendResidue(chain, "U");

  e.InsertAtom(res, "P",   geom::Vec3(16.249, 28.254, 44.759));
  e.InsertAtom(res, "OP1", geom::Vec3(16.654, 28.055, 46.181));
  e.InsertAtom(res, "OP2", geom::Vec3(17.115, 27.731, 43.656));
  e.InsertAtom(res, "O5'", geom::Vec3(14.826, 27.572, 44.574));
  e.InsertAtom(res, "C5'", geom::Vec3(14.711, 26.206, 44.216));
  e.InsertAtom(res, "C4'", geom::Vec3(13.279, 25.889, 43.887));
  e.InsertAtom(res, "O4'", geom::Vec3(12.832, 26.793, 42.838));
  e.InsertAtom(res, "C3'", geom::Vec3(13.155, 24.434, 43.329));
  e.InsertAtom(res, "O3'", geom::Vec3(12.269, 23.625, 44.098));
  e.InsertAtom(res, "C2'", geom::Vec3(12.871, 24.595, 41.875));
  e.InsertAtom(res, "O2'", geom::Vec3(11.811, 23.752, 41.462));
  e.InsertAtom(res, "C1'", geom::Vec3(12.424, 26.056, 41.694));
  e.InsertAtom(res, "N1",  geom::Vec3(13.030, 26.692, 40.497));
  e.InsertAtom(res, "C2",  geom::Vec3(12.517, 26.365, 39.228));
  e.InsertAtom(res, "O2",  geom::Vec3(11.579, 25.594, 39.068));
  e.InsertAtom(res, "N3",  geom::Vec3(13.141, 26.987, 38.161));
  e.InsertAtom(res, "C4",  geom::Vec3(14.197, 27.888, 38.210));
  e.InsertAtom(res, "O4",  geom::Vec3(14.627, 28.368, 37.156));
  e.InsertAtom(res, "C5",  geom::Vec3(14.671, 28.189, 39.542));
  e.InsertAtom(res, "C6",  geom::Vec3(14.087, 27.597, 40.612));

  return res;
}

ResidueHandle make_defective_uracil2(ChainHandle chain) 
{
  XCSEditor e=chain.GetEntity().EditXCS();
  ResidueHandle res = e.AppendResidue(chain, "U");

  e.InsertAtom(res, "P",   geom::Vec3(16.249, 28.254, 44.759));
  e.InsertAtom(res, "OP1", geom::Vec3(16.654, 28.055, 46.181));
  e.InsertAtom(res, "OP2", geom::Vec3(17.115, 27.731, 43.656));
  e.InsertAtom(res, "O5'", geom::Vec3(14.826, 27.572, 44.574));
  e.InsertAtom(res, "C5'", geom::Vec3(14.711, 26.206, 44.216));
  e.InsertAtom(res, "C4'", geom::Vec3(13.279, 25.889, 43.887));
  e.InsertAtom(res, "O4'", geom::Vec3(12.832, 26.793, 42.838));
  e.InsertAtom(res, "C3'", geom::Vec3(13.155, 24.434, 43.329));
  e.InsertAtom(res, "O3'", geom::Vec3(12.269, 23.625, 44.098));
  e.InsertAtom(res, "C2'", geom::Vec3(12.871, 24.595, 41.875));
  e.InsertAtom(res, "O2'", geom::Vec3(11.811, 23.752, 41.462));
  e.InsertAtom(res, "C1'", geom::Vec3(12.424, 26.056, 41.694));
  e.InsertAtom(res, "N1",  geom::Vec3(13.030, 26.692, 40.497));
  e.InsertAtom(res, "C2",  geom::Vec3(12.517, 26.365, 39.228));
  e.InsertAtom(res, "O2",  geom::Vec3(11.579, 25.594, 39.068));
  e.InsertAtom(res, "N3",  geom::Vec3(13.141, 26.987, 38.161));
  e.InsertAtom(res, "C4",  geom::Vec3(14.197, 27.888, 38.210));
  e.InsertAtom(res, "O4",  geom::Vec3(14.627, 28.368, 37.156));
  e.InsertAtom(res, "C5",  geom::Vec3(14.671, 28.189, 39.542));
  e.InsertAtom(res, "C6",  geom::Vec3(14.087, 27.597, 80.612));

  return res;
}
ResidueHandle make_1zk(ChainHandle chain)
{
  XCSEditor edi=chain.GetEntity().EditXCS();
  ResidueHandle r=edi.AppendResidue(chain, "1ZK");
  edi.InsertAtom(r, "C", geom::Vec3(3.946, 3.520, 10.861),"C");
  edi.InsertAtom(r, "C1", geom::Vec3(6.308, 3.824, 8.036),"C");
  edi.InsertAtom(r, "CA", geom::Vec3(4.856, 4.282, 9.931),"C");
  edi.InsertAtom(r, "C2", geom::Vec3(6.770, 5.099, 8.019),"C");
  edi.InsertAtom(r, "C3", geom::Vec3(7.805, 5.484, 7.142),"C");
  edi.InsertAtom(r, "C4", geom::Vec3(8.363, 4.590, 6.284),"C");
  edi.InsertAtom(r, "C5", geom::Vec3(8.469, 2.345, 5.318),"C");
  edi.InsertAtom(r, "C6", geom::Vec3(7.986, 1.072, 5.322),"C");
  edi.InsertAtom(r, "C4A", geom::Vec3(7.885, 3.277, 6.208),"C");
  edi.InsertAtom(r, "C7", geom::Vec3(6.949, 0.682, 6.190),"C");
  edi.InsertAtom(r, "C8", geom::Vec3(6.398, 1.548, 7.074),"C");
  edi.InsertAtom(r, "C9", geom::Vec3(4.018, 0.269, 12.183),"C");
  edi.InsertAtom(r, "C1A", geom::Vec3(6.824, 2.882, 7.117),"C");
  edi.InsertAtom(r, "O", geom::Vec3(3.224, 4.213, 11.565),"O");
  edi.InsertAtom(r, "O1", geom::Vec3(5.331, 3.476, 8.904),"O");
  edi.InsertAtom(r, "O3", geom::Vec3(5.748, -4.044, 14.612),"O");
  edi.InsertAtom(r, "N", geom::Vec3(3.855, 2.201, 10.790),"N");
  edi.InsertAtom(r, "N1", geom::Vec3(4.814, 0.719, 13.173),"N");
  edi.InsertAtom(r, "CA1", geom::Vec3(3.077, 1.346, 11.644),"C");
  edi.InsertAtom(r, "O2", geom::Vec3(4.203, -0.810, 11.651),"O");
  edi.InsertAtom(r, "O4", geom::Vec3(7.311, -5.667, 18.880),"O");
  edi.InsertAtom(r, "CB", geom::Vec3(1.856, 0.712, 11.049),"C");
  edi.InsertAtom(r, "CG", geom::Vec3(1.015, 1.845, 10.511),"C");
  edi.InsertAtom(r, "ND1", geom::Vec3(1.467, 2.439, 9.321),"N");
  edi.InsertAtom(r, "N2", geom::Vec3(6.478, -3.958, 16.751),"N");
  edi.InsertAtom(r, "CD2", geom::Vec3(-0.105, 2.479, 10.887),"C");
  edi.InsertAtom(r, "CE1", geom::Vec3(0.638, 3.428, 9.002),"C");
  edi.InsertAtom(r, "NE2", geom::Vec3(-0.372, 3.461, 9.881),"N");
  edi.InsertAtom(r, "N3", geom::Vec3(6.871, -7.136, 17.332),"N");
  edi.InsertAtom(r, "CA2", geom::Vec3(5.881, -0.001, 13.808),"C");
  edi.InsertAtom(r, "CB1", geom::Vec3(7.140, 0.860, 13.743),"C");
  edi.InsertAtom(r, "CG1", geom::Vec3(7.503, 1.324, 12.299),"C");
  edi.InsertAtom(r, "C21", geom::Vec3(5.126, -8.557, 18.179),"C");
  edi.InsertAtom(r, "CD1", geom::Vec3(8.185, 0.142, 11.507),"C");
  edi.InsertAtom(r, "CD21", geom::Vec3(8.420, 2.537, 12.325),"C");
  edi.InsertAtom(r, "CE11", geom::Vec3(8.381, 0.689, 10.066),"C");
  edi.InsertAtom(r, "CE2", geom::Vec3(8.907, 2.979, 10.922),"C");
  edi.InsertAtom(r, "CZ", geom::Vec3(9.409, 1.807, 10.075),"C");
  edi.InsertAtom(r, "CH", geom::Vec3(5.592, -0.511, 15.204),"C");
  edi.InsertAtom(r, "OH", geom::Vec3(5.225, 0.377, 16.238),"O");
  edi.InsertAtom(r, "CB11", geom::Vec3(4.426, -1.543, 15.170),"C");
  edi.InsertAtom(r, "CA'", geom::Vec3(4.451, -2.730, 16.152),"C");
  edi.InsertAtom(r, "CB'", geom::Vec3(3.124, -3.441, 16.281),"C");
  edi.InsertAtom(r, "CG11", geom::Vec3(2.553, -3.986, 14.933),"C");
  edi.InsertAtom(r, "C31", geom::Vec3(4.413, -7.811, 19.117),"C");
  edi.InsertAtom(r, "CG2", geom::Vec3(3.204, -4.586, 17.345),"C");
  edi.InsertAtom(r, "OB1", geom::Vec3(3.249, -0.875, 15.134),"O");
  edi.InsertAtom(r, "CC", geom::Vec3(5.603, -3.655, 15.782),"C");
  edi.InsertAtom(r, "CA3", geom::Vec3(7.592, -4.867, 16.603),"C");
  edi.InsertAtom(r, "CD", geom::Vec3(7.274, -5.947, 17.691),"C");
  edi.InsertAtom(r, "CB2", geom::Vec3(8.986, -4.351, 16.803),"C");
  edi.InsertAtom(r, "CG12", geom::Vec3(9.488, -3.108, 16.016),"C");
  edi.InsertAtom(r, "CG21", geom::Vec3(10.000, -5.461, 16.472),"C");
  edi.InsertAtom(r, "CD11", geom::Vec3(9.099, -3.257, 14.571),"C");
  edi.InsertAtom(r, "CM", geom::Vec3(6.587, -8.286, 18.106),"C");
  edi.InsertAtom(r, "C41", geom::Vec3(3.045, -7.980, 19.287),"C");
  edi.InsertAtom(r, "C51", geom::Vec3(2.423, -8.911, 18.456),"C");
  edi.InsertAtom(r, "C61", geom::Vec3(3.164, -9.631, 17.518),"C");
  edi.InsertAtom(r, "N11", geom::Vec3(4.497, -9.459, 17.386),"N");
  return r;
}

void verify_nucleotide_connectivity(const ResidueHandle& res) 
{
  BOOST_CHECK(BondExists(res.FindAtom("P"),
                         res.FindAtom("OP1")));
  BOOST_CHECK(BondExists(res.FindAtom("P"),
                         res.FindAtom("OP2")));
  BOOST_CHECK(BondExists(res.FindAtom("P"),
                         res.FindAtom("O5'")));
  BOOST_CHECK(BondExists(res.FindAtom("O5'"),
                         res.FindAtom("C5'")));
  BOOST_CHECK(BondExists(res.FindAtom("C5'"),
                         res.FindAtom("C4'")));
  BOOST_CHECK(BondExists(res.FindAtom("C4'"),
                         res.FindAtom("O4'")));
  BOOST_CHECK(BondExists(res.FindAtom("C4'"),
                         res.FindAtom("C3'")));
  BOOST_CHECK(BondExists(res.FindAtom("C3'"),
                         res.FindAtom("O3'")));
  BOOST_CHECK(BondExists(res.FindAtom("C3'"),
                         res.FindAtom("C2'")));
  BOOST_CHECK(BondExists(res.FindAtom("C2'"),
                         res.FindAtom("C1'")));
  BOOST_CHECK(BondExists(res.FindAtom("C1'"),
                         res.FindAtom("O4'")));

  if (res.GetKey()=="DC") {
    BOOST_CHECK(BondExists(res.FindAtom("C1'"),
                           res.FindAtom("N1")));
    BOOST_CHECK(BondExists(res.FindAtom("N1"),
                           res.FindAtom("C2")));
    BOOST_CHECK(BondExists(res.FindAtom("C2"),
                           res.FindAtom("O2")));
    BOOST_CHECK(BondExists(res.FindAtom("C2"),
                           res.FindAtom("N3")));
    BOOST_CHECK(BondExists(res.FindAtom("N3"),
                           res.FindAtom("C4")));
    BOOST_CHECK(BondExists(res.FindAtom("C4"),
                           res.FindAtom("N4")));
    BOOST_CHECK(BondExists(res.FindAtom("C4"),
                           res.FindAtom("C5")));
    BOOST_CHECK(BondExists(res.FindAtom("C5"),
                           res.FindAtom("C6")));
    BOOST_CHECK(BondExists(res.FindAtom("C6"),
                           res.FindAtom("N1")));
    // TODO: Check that no other atoms are connected!
  }

  if (res.GetKey()=="U") {
    BOOST_CHECK(BondExists(res.FindAtom("O2'"),
                           res.FindAtom("C2'")));
    BOOST_CHECK(BondExists(res.FindAtom("C1'"),
                           res.FindAtom("N1")));
    BOOST_CHECK(BondExists(res.FindAtom("N1"),
                           res.FindAtom("C2")));
    BOOST_CHECK(BondExists(res.FindAtom("C2"),
                           res.FindAtom("O2")));
    BOOST_CHECK(BondExists(res.FindAtom("C2"),
                           res.FindAtom("N3")));
    BOOST_CHECK(BondExists(res.FindAtom("N3"),
                           res.FindAtom("C4")));
    BOOST_CHECK(BondExists(res.FindAtom("C4"),
                           res.FindAtom("O4")));
    BOOST_CHECK(BondExists(res.FindAtom("C4"),
                           res.FindAtom("C5")));
    BOOST_CHECK(BondExists(res.FindAtom("C5"),
                           res.FindAtom("C6")));
    BOOST_CHECK(BondExists(res.FindAtom("C6"),
                           res.FindAtom("N1")));
    // TODO: Check that no other atoms are connected!
  }
}
void verify_1zk_connectivity(const ResidueHandle& r1) 
{
  BOOST_CHECK(BondExists(r1.FindAtom("C"), r1.FindAtom("CA")));
  BOOST_CHECK(BondExists(r1.FindAtom("C"), r1.FindAtom("O")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA"), r1.FindAtom("O1")));
  BOOST_CHECK(BondExists(r1.FindAtom("C1"), r1.FindAtom("C2")));
  BOOST_CHECK(BondExists(r1.FindAtom("C1"), r1.FindAtom("C1A")));
  BOOST_CHECK(BondExists(r1.FindAtom("C1"), r1.FindAtom("O1")));
  BOOST_CHECK(BondExists(r1.FindAtom("C2"), r1.FindAtom("C3")));
  BOOST_CHECK(BondExists(r1.FindAtom("C3"), r1.FindAtom("C4")));
  BOOST_CHECK(BondExists(r1.FindAtom("C4"), r1.FindAtom("C4A")));
  BOOST_CHECK(BondExists(r1.FindAtom("C4A"), r1.FindAtom("C5")));
  BOOST_CHECK(BondExists(r1.FindAtom("C4A"), r1.FindAtom("C1A")));
  BOOST_CHECK(BondExists(r1.FindAtom("C5"), r1.FindAtom("C6")));
  BOOST_CHECK(BondExists(r1.FindAtom("C6"), r1.FindAtom("C7")));
  BOOST_CHECK(BondExists(r1.FindAtom("C7"), r1.FindAtom("C8")));
  BOOST_CHECK(BondExists(r1.FindAtom("C8"), r1.FindAtom("C1A")));
  BOOST_CHECK(BondExists(r1.FindAtom("N"), r1.FindAtom("CA1")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA1"), r1.FindAtom("C9")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA1"), r1.FindAtom("CB")));
  BOOST_CHECK(BondExists(r1.FindAtom("C9"), r1.FindAtom("O2")));
  BOOST_CHECK(BondExists(r1.FindAtom("CB"), r1.FindAtom("CG")));
  BOOST_CHECK(BondExists(r1.FindAtom("CG"), r1.FindAtom("ND1")));
  BOOST_CHECK(BondExists(r1.FindAtom("CG"), r1.FindAtom("CD2")));
  BOOST_CHECK(BondExists(r1.FindAtom("ND1"), r1.FindAtom("CE1")));
  BOOST_CHECK(BondExists(r1.FindAtom("CD2"), r1.FindAtom("NE2")));
  BOOST_CHECK(BondExists(r1.FindAtom("CE1"), r1.FindAtom("NE2")));
  BOOST_CHECK(BondExists(r1.FindAtom("N1"), r1.FindAtom("CA2")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA2"), r1.FindAtom("CB1")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA2"), r1.FindAtom("CH")));
  BOOST_CHECK(BondExists(r1.FindAtom("CB1"), r1.FindAtom("CG1")));
  BOOST_CHECK(BondExists(r1.FindAtom("CG1"), r1.FindAtom("CD1")));
  BOOST_CHECK(BondExists(r1.FindAtom("CG1"), r1.FindAtom("CD21")));
  BOOST_CHECK(BondExists(r1.FindAtom("CD1"), r1.FindAtom("CE11")));
  BOOST_CHECK(BondExists(r1.FindAtom("CD21"), r1.FindAtom("CE2")));
  BOOST_CHECK(BondExists(r1.FindAtom("CE11"), r1.FindAtom("CZ")));
  BOOST_CHECK(BondExists(r1.FindAtom("CE2"), r1.FindAtom("CZ")));
  BOOST_CHECK(BondExists(r1.FindAtom("CH"), r1.FindAtom("OH")));
  BOOST_CHECK(BondExists(r1.FindAtom("CH"), r1.FindAtom("CB11")));
  BOOST_CHECK(BondExists(r1.FindAtom("CB11"), r1.FindAtom("CA'")));
  BOOST_CHECK(BondExists(r1.FindAtom("CB11"), r1.FindAtom("OB1")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA'"), r1.FindAtom("CB'")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA'"), r1.FindAtom("CC")));
  BOOST_CHECK(BondExists(r1.FindAtom("CB'"), r1.FindAtom("CG11")));
  BOOST_CHECK(BondExists(r1.FindAtom("CB'"), r1.FindAtom("CG2")));
  BOOST_CHECK(BondExists(r1.FindAtom("CC"), r1.FindAtom("O3")));
  BOOST_CHECK(BondExists(r1.FindAtom("N2"), r1.FindAtom("CA3")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA3"), r1.FindAtom("CD")));
  BOOST_CHECK(BondExists(r1.FindAtom("CA3"), r1.FindAtom("CB2")));
  BOOST_CHECK(BondExists(r1.FindAtom("CD"), r1.FindAtom("O4")));
  BOOST_CHECK(BondExists(r1.FindAtom("CB2"), r1.FindAtom("CG12")));
  BOOST_CHECK(BondExists(r1.FindAtom("CB2"), r1.FindAtom("CG21")));
  BOOST_CHECK(BondExists(r1.FindAtom("CG12"), r1.FindAtom("CD11")));
  BOOST_CHECK(BondExists(r1.FindAtom("N3"), r1.FindAtom("CM")));
  BOOST_CHECK(BondExists(r1.FindAtom("CM"), r1.FindAtom("C21")));
  BOOST_CHECK(BondExists(r1.FindAtom("C21"), r1.FindAtom("C31")));
  BOOST_CHECK(BondExists(r1.FindAtom("C21"), r1.FindAtom("N11")));
  BOOST_CHECK(BondExists(r1.FindAtom("C31"), r1.FindAtom("C41")));
  BOOST_CHECK(BondExists(r1.FindAtom("C41"), r1.FindAtom("C51")));
  BOOST_CHECK(BondExists(r1.FindAtom("C51"), r1.FindAtom("C61")));
  BOOST_CHECK(BondExists(r1.FindAtom("C61"), r1.FindAtom("N11")));
  BOOST_CHECK(BondExists(r1.FindAtom("C"), r1.FindAtom("N")));
  BOOST_CHECK(BondExists(r1.FindAtom("C9"), r1.FindAtom("N1")));
  BOOST_CHECK(BondExists(r1.FindAtom("CC"), r1.FindAtom("N2")));
  BOOST_CHECK(BondExists(r1.FindAtom("CD"), r1.FindAtom("N3")));
}

ResidueHandle make_arg(ChainHandle chain) 
{
  XCSEditor e=chain.GetEntity().EditXCS();
  ResidueHandle res = e.AppendResidue(chain, "ARG");
  e.InsertAtom(res, "N",geom::Vec3(20.202,33.112,58.011));
  e.InsertAtom(res, "CA",geom::Vec3(19.396,31.903,58.033));
  e.InsertAtom(res, "C",geom::Vec3(18.608,31.739,59.328));
  e.InsertAtom(res, "O",geom::Vec3(17.651,30.965,59.381));
  e.InsertAtom(res, "CB",geom::Vec3(20.284,30.681,57.801));
  e.InsertAtom(res, "CG",geom::Vec3(20.665,30.488,56.342));
  e.InsertAtom(res, "CD",geom::Vec3(21.557,29.281,56.154));
  e.InsertAtom(res, "NE",geom::Vec3(22.931,29.557,56.551));
  e.InsertAtom(res, "CZ",geom::Vec3(23.901,28.653,56.528));
  e.InsertAtom(res, "NH1",geom::Vec3(23.640,27.417,56.130));
  e.InsertAtom(res, "NH2",geom::Vec3(25.132,28.980,56.893));
  return res;
}

ResidueHandle make_leu(ChainHandle chain) 
{
  XCSEditor e=chain.GetEntity().EditXCS();  
  ResidueHandle res=e.AppendResidue(chain, "LEU");

  e.InsertAtom(res, "N", geom::Vec3(19.003,32.473,60.366));
  e.InsertAtom(res, "CA", geom::Vec3(18.330,32.402,61.664));
  e.InsertAtom(res, "C", geom::Vec3(17.884,33.787,62.117));
  e.InsertAtom(res, "O", geom::Vec3(17.853,34.091,63.308));
  e.InsertAtom(res, "CB", geom::Vec3(19.269,31.793,62.710));
  e.InsertAtom(res, "CG", geom::Vec3(19.695,30.340,62.501));
  e.InsertAtom(res, "CD1", geom::Vec3(20.585,29.897,63.648));
  e.InsertAtom(res, "CD2", geom::Vec3(18.461,29.459,62.420));
  return res;
}

ResidueHandle make_defective_leu(ChainHandle chain) 
{
  XCSEditor e=chain.GetEntity().EditXCS();  
  ResidueHandle res=e.AppendResidue(chain, "LEU");

  e.InsertAtom(res, "N", geom::Vec3(19.003,32.473,60.366));
  e.InsertAtom(res, "CA", geom::Vec3(18.330,32.402,61.664));
  e.InsertAtom(res, "C", geom::Vec3(17.884,33.787,62.117));
  e.InsertAtom(res, "O", geom::Vec3(17.853,34.091,63.308));
  e.InsertAtom(res, "CB", geom::Vec3(19.269,31.793,102.710));
  e.InsertAtom(res, "CG", geom::Vec3(19.695,30.340,62.501));
  e.InsertAtom(res, "CD1", geom::Vec3(20.585,29.897,63.648));
  e.InsertAtom(res, "CD2", geom::Vec3(18.461,29.459,62.420));
  return res;
}

}}

