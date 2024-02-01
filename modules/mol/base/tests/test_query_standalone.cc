#include <ost/mol/mol.hh>

using namespace ost;
using namespace ost::mol;

int main()
{
  EntityHandle eh = CreateEntity();
  XCSEditor e=eh.EditXCS();
  ChainHandle chain = e.InsertChain("A");
  ResidueHandle res = e.AppendResidue(chain, "MET");
  e.InsertAtom(res, "N",geom::Vec3(21.609,35.384,56.705), "N");
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

  EntityView v=eh.Select("rname=MET and aname=C*");
  std::cerr << v.GetAtomCount() << std::endl;
  v=eh.Select("rname=ARG and aname=N?1");
  std::cerr << v.GetAtomCount() << std::endl;
  v=eh.Select("rname=ARG and aname=NH?");
  std::cerr << v.GetAtomCount() << std::endl;
  v=eh.Select("rname=ARG and aname=\"*2\"");
  std::cerr << v.GetAtomCount() << std::endl;
  v=eh.Select("rname=ARG and aname=N?");
  std::cerr << v.GetAtomCount() << std::endl;
}
