#include <sstream>
#include <cmath>

#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>

#include <ost/log.hh>
#include <ost/platform.hh>
#include <ost/img/alg/dft.hh>
#include <ost/mol/mol.hh>

#include "entity_to_density.hh"

namespace ost { namespace mol { namespace alg {

namespace detail {

struct AtomScatteringProps {
 String element;
  Real atomic_weight;
  Real c;
  Real a1;
  Real a2;
  Real a3;
  Real a4;
  Real b1;
  Real b2;
  Real b3;
  Real b4;
};

typedef std::vector<AtomScatteringProps> AtomScatteringPropsTable;

void FillAtomScatteringPropsTable(AtomScatteringPropsTable& table)
{
  AtomScatteringProps atom_props;

  String ost_root=GetSharedDataPath();
  String filename = "/atom_scattering_properties.txt";
  String fullpath = ost_root+filename;
  boost::filesystem::path loc(fullpath);
  boost::filesystem::ifstream infile(loc);
  if (!infile) {
    String msg;
    std::stringstream msgstr(msg);
    msgstr << "Couldn't find " << fullpath << std::endl;
    throw ost::Error(msg);
  }
  String line;
  std::getline(infile, line);
  std::getline(infile, line);
  std::getline(infile, line);
  std::getline(infile, line);
  std::getline(infile, line);
  while (std::getline(infile, line))
  {
   std::stringstream line_stream(line);
   line_stream >> atom_props.element
               >> atom_props.atomic_weight
               >> atom_props.c
               >> atom_props.a1
               >> atom_props.a2
               >> atom_props.a3
               >> atom_props.a4
               >> atom_props.b1
               >> atom_props.b2
               >> atom_props.b3
               >> atom_props.b4;
    table.push_back(atom_props);
  }
};


Real ScatteringFactor (Real frequency, const AtomScatteringProps& props,
                         Real source_wavelength)
{
   Real scattering_factor =
             props.a1*exp(-props.b1*frequency*frequency) +
             props.a2*exp(-props.b2*frequency*frequency) +
             props.a3*exp(-props.b3*frequency*frequency) +
             props.a4*exp(-props.b4*frequency*frequency) +
             props.c;

   return scattering_factor;
}


class EntityToDensityHelperBase
{

public:

  EntityToDensityHelperBase (const ost::mol::EntityView entity_view,
                             Real falloff_start,
                             Real falloff_end,
                             Real source_wavelength,
                             geom::Vec3 map_start,
                             geom::Vec3 map_end
                             ):
      entity_view_(entity_view),
      falloff_start_frequency_(1.0/falloff_start),
      falloff_end_frequency_(1.0/falloff_end),
      source_wavelength_(source_wavelength),
      map_start_(map_start),map_end_(map_end)
  {
    if (scattering_props_table_loaded_ == false)
    {
      FillAtomScatteringPropsTable(scatt_props_table_);
      scattering_props_table_loaded_ = true;
    }
  }

  void VisitState(img::ComplexHalfFrequencyImageState& is) const
  {

    geom::Vec3 frequency_sampling =
             is.GetSampling().GetFrequencySampling();

    Real minus_2_pi = -2.0*M_PI;

    uint x_limit = ceil(falloff_end_frequency_ / frequency_sampling[0]);
    uint y_limit = ceil(falloff_end_frequency_ / frequency_sampling[1]);
    uint z_limit = ceil(falloff_end_frequency_ / frequency_sampling[2]);
     img::Extent reduced_extent = img::Extent
             (img::Point(-static_cast<int>(x_limit),-static_cast<int>(y_limit),0),
              img::Point(x_limit,y_limit,z_limit));


    for (ChainViewList::const_iterator ci = entity_view_.GetChainList().begin(),
        ce = entity_view_.GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
          re = ci->GetResidueList().end(); ri != re; ++ri) {
    for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
          ae = ri->GetAtomList().end(); ai!= ae; ++ai) {
      AtomScatteringPropsTable::iterator table_iter =
                                             scatt_props_table_.begin();
      bool found = false;
      while (found != true && table_iter!=scatt_props_table_.end())
      {
        if ( (*table_iter).element == (*ai).GetElement())
        {        
          geom::Vec3 coord = ai->GetPos();

          if (coord[0] >= map_start_[0] &&
              coord[0] <= map_end_[0] &&
              coord[1] >= map_start_[1] &&
              coord[1] <= map_end_[1] &&
              coord[2] >= map_start_[2] &&
              coord[2] <= map_end_[2])
          {

            // This part of the code assumes that the three axes
            // of the map are at right angles and the origin at 0,0,0.
            // Eventually, when maps and images will be merged,
            // it will substituted by the map's xyz2uvw method

            geom::Vec3 adjusted_coord  = coord-map_start_;

            for (img::ExtentIterator mp_it(reduced_extent);
               !mp_it.AtEnd();++mp_it)
            {
              img::Point mp_it_point = img::Point(mp_it);

              geom::Vec3 frequency_sampling =
                  is.GetSampling().GetFrequencySampling();

              geom::Vec3 mp_it_vec = geom::Vec3(
                (mp_it_point[0]*frequency_sampling[0]),
                (mp_it_point[1]*frequency_sampling[1]),
                (mp_it_point[2]*frequency_sampling[2])
              );

              Real frequency = Length(mp_it_vec);

              Real sigma = (falloff_end_frequency_ -
                              falloff_start_frequency_)/3.0;

              if (frequency <= falloff_end_frequency_)
              {
                Real falloff_term = 1.0;
                if (sigma!=0 && frequency >= falloff_start_frequency_)
                {
                  falloff_term=exp(-(frequency-falloff_start_frequency_)*
                                        (frequency-falloff_start_frequency_)/
                                        (2.0*sigma*sigma));
                }
                AtomScatteringProps scatt_props = (*table_iter);
                Real scatt_fact = ScatteringFactor(frequency,
                                                     scatt_props,
                                                     source_wavelength_);

                Real amp_term = scatt_fact*falloff_term;
                Real exp = minus_2_pi * geom::Dot(mp_it_vec,adjusted_coord);
                Real real = amp_term *std::cos(exp);
                Real imag = amp_term *std::sin(exp);
                is.Value(mp_it) = is.Value(mp_it)+Complex(real,imag);

              }
            }
          }
        }
        ++table_iter;
      }
    }}} // loop over atoms
  }

  template <typename T, class D>
  void VisitState(img::ImageStateImpl<T,D>& is) const
  {
    assert(false);
  }

  static String GetAlgorithmName() {return "EntityToDensityHelper"; }

private:

  ost::mol::alg::DensityType density_type_;
  Real limit_;
  ost::mol::EntityView entity_view_;
  Real falloff_start_frequency_;
  Real falloff_end_frequency_;
  static AtomScatteringPropsTable scatt_props_table_;
  static bool scattering_props_table_loaded_;
  Real source_wavelength_;
  geom::Vec3 map_start_;
  geom::Vec3 map_end_;

};

AtomScatteringPropsTable EntityToDensityHelperBase::scatt_props_table_ =
          AtomScatteringPropsTable();
bool EntityToDensityHelperBase::scattering_props_table_loaded_ = false;

typedef img::ImageStateConstModIPAlgorithm<EntityToDensityHelperBase>
        EntityToDensityHelper;



class EntityToDensityRosettaHelperBase
{

public:

  EntityToDensityRosettaHelperBase (const ost::mol::EntityView entity_view,
                             const ost::mol::alg::DensityType& density_type,
                             Real resolution):
      density_type_(density_type), entity_view_(entity_view),
      resolution_(resolution)
  {
    if (scattering_props_table_loaded_ == false)
    {
     FillAtomScatteringPropsTable(scatt_props_table_);
     scattering_props_table_loaded_ = true;
    }
  }

  void VisitState(img::RealSpatialImageState& is) const
  {
    geom::Vec3 map_start=is.IndexToCoord(is.GetExtent().GetStart());
    geom::Vec3 map_end=is.IndexToCoord(is.GetExtent().GetEnd());
    mol::EntityView effective_entity_view=entity_view_;
    Real k,C;

    if (density_type_ == HIGH_RESOLUTION)
    {
      k = (M_PI*M_PI)/(resolution_*resolution_);
      C = sqrt((M_PI*M_PI*M_PI)/(k*k*k));
      effective_entity_view = entity_view_;

    } else {  // ROSETTA LOW RESOLUTION

      k = (M_PI*M_PI)/((2.4+0.8*resolution_)*(2.4+0.8*resolution_));
      C = sqrt((M_PI*M_PI*M_PI)/(k*k*k));

      effective_entity_view = entity_view_.Select("aname=CA and peptide=true");
    }

    Real sigma = sqrt(1.0/(2.0*k));
    Real four_sigma_squared = 16.0*sigma*sigma;
    Real four_sigma = 4.0*sigma;

    geom::Vec3 sampling = is.GetSampling().GetPixelSampling();
    img::Extent is_extent = is.GetExtent();

    for (ChainViewList::const_iterator ci = effective_entity_view.GetChainList().begin(),
        ce = effective_entity_view.GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
        re = ci->GetResidueList().end(); ri != re; ++ri) {
    for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
        ae = ri->GetAtomList().end(); ai!= ae; ++ai) {
      AtomScatteringPropsTable::iterator table_iter =
                                             scatt_props_table_.begin();
      bool found = false;
      while (found != true && table_iter!=scatt_props_table_.end()) {
        if ((*table_iter).element == ai->GetElement()) {
          found = true;
          Real a = (*table_iter).atomic_weight;
          geom::Vec3 coord = ai->GetPos();
          if (coord[0] >= map_start[0] &&
              coord[0] <= map_end[0] &&
              coord[1] >= map_start[1] &&
              coord[1] <= map_end[1] &&
              coord[2] >= map_start[2] &&
              coord[2] <= map_end[2])
          {

            geom::Vec3 adjusted_coord  = coord-map_start;
            geom::Vec3 pixel_coord=is.CoordToIndex(coord);
            img::Point rounded_pixel_coord(round(pixel_coord[0]),
                                            round(pixel_coord[1]),
                                            round(pixel_coord[2]));

            uint x_limit = ceil(2.0*four_sigma/sampling[0])+1;
            uint y_limit = ceil(2.0*four_sigma/sampling[1])+1;
            uint z_limit = ceil(2.0*four_sigma/sampling[2])+1;

            img::Extent reduced_extent = img::Extent
                                      (img::Size(x_limit,y_limit,z_limit),
                                      rounded_pixel_coord);

            img::Extent iteration_extent=
                                  img::Overlap(is_extent,reduced_extent);

            for (img::ExtentIterator mp_it(iteration_extent);
              !mp_it.AtEnd();++mp_it)
            {

              img::Point mp_it_point = img::Point(mp_it);

              geom::Vec3 mp_it_vec = geom::Vec3(
                (mp_it_point[0]*sampling[0]),
                (mp_it_point[1]*sampling[1]),
                (mp_it_point[2]*sampling[2])
              );

              Real distance_squared = Length2(mp_it_vec-adjusted_coord);

              if (distance_squared <= four_sigma_squared)
              {
                Real exp_term = exp(-k*distance_squared);
                Real value = C * a * exp_term;
                is.Value(mp_it_point) = is.Value(mp_it_point) + value;
              }
            }
          }
        }
        ++table_iter;
      }
    }}}
  }

  template <typename T, class D>
  void VisitState(img::ImageStateImpl<T,D>& is) const
  {
    assert(false);
  }

  static String GetAlgorithmName() { return "EntityToDensityRosettaHelper"; }

private:

  ost::mol::alg::DensityType density_type_;
  ost::mol::EntityView entity_view_;
  Real resolution_;
  static AtomScatteringPropsTable scatt_props_table_;
  static bool scattering_props_table_loaded_;
};

AtomScatteringPropsTable EntityToDensityRosettaHelperBase::scatt_props_table_ =
          AtomScatteringPropsTable();
bool EntityToDensityRosettaHelperBase::scattering_props_table_loaded_ = false;

typedef img::ImageStateConstModIPAlgorithm<EntityToDensityRosettaHelperBase>
    EntityToDensityRosettaHelper;


} // namespace

void EntityToDensityScattering(const mol::EntityView& entity_view,
                                     img::MapHandle& map,
                                     Real falloff_start,
                                     Real falloff_end,
                                     bool clear_map_flag,
                                     Real source_wavelength)
{
  if(falloff_start<=0.0) throw ost::Error("Invalid falloff start");
  if(falloff_end<=0.0 || falloff_end>falloff_start)
     throw ost::Error("Invalid falloff end");

 geom ::Vec3 rs_sampl = map.GetSpatialSampling();
  geom ::Vec3 abs_orig = map.GetAbsoluteOrigin();
  geom::Vec3 map_start = geom::Vec3(abs_orig[0]+map.GetExtent().GetStart()[0]*rs_sampl[0],
                                    abs_orig[1]+map.GetExtent().GetStart()[1]*rs_sampl[1],
                                    abs_orig[2]+map.GetExtent().GetStart()[2]*rs_sampl[2]);

  geom::Vec3 map_end = geom::Vec3(abs_orig[0]+map.GetExtent().GetEnd()[0]*rs_sampl[0],
                                  abs_orig[1]+map.GetExtent().GetEnd()[1]*rs_sampl[1],
                                  abs_orig[2]+map.GetExtent().GetEnd()[2]*rs_sampl[2]);
  detail::EntityToDensityHelper e_to_d_helper(entity_view,
                                              falloff_start,
                                              falloff_end,
                                              source_wavelength, map_start,map_end);
  if (clear_map_flag==true) {
    img::MapHandle mm=img::CreateImage(img::Extent(map.GetSize(),
                                       img::Point(0,0)),
                                       img::COMPLEX,img::HALF_FREQUENCY);
    // swap newly created map into place
    mm.SetSpatialSampling(map.GetSpatialSampling());
    mm.SetAbsoluteOrigin(map.GetAbsoluteOrigin());
    map.Swap(mm);
  } else {
    map.ApplyIP(img::alg::DFT());
  }

  map.ApplyIP(e_to_d_helper);
  map.ApplyIP(img::alg::DFT());
}


void EntityToDensityRosetta(const mol::EntityView& entity_view,
                            img::MapHandle& map,
                            const DensityType& density_type,
                            Real resolution,
                            bool clear_map_flag,
                            Real source_wavelength)

{
  if(resolution <=0.0) throw ost::Error("Invalid resolution");
  if (clear_map_flag==true) {
    img::MapHandle mm=img::CreateImage(img::Extent(img::Point(0,0),
                                       map.GetSize()));
    // swap newly created map into place
    mm.SetSpatialSampling(map.GetSpatialSampling());
    mm.SetAbsoluteOrigin(map.GetAbsoluteOrigin());
    map.Swap(mm);
  }

  detail::EntityToDensityRosettaHelper e_to_d_r_helper
                                                (entity_view,density_type,
                                                resolution);
  map.ApplyIP(e_to_d_r_helper);
}


}}} // ns


