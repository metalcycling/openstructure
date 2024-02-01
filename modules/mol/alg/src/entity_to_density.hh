#ifndef OST_ENTITY_TO_DENSITY_HH
#define OST_ENTITY_TO_DENSITY_HH

#include <ost/mol/entity_view.hh>
#include <ost/img/map.hh>

#include <ost/mol/alg/module_config.hh>


namespace ost { namespace mol { namespace alg {

/// \brief type of density being created by the EntityToDensity function
enum DensityType
{
  HIGH_RESOLUTION,
  LOW_RESOLUTION
};

/// \brief create a density representation of an entity in a density map (using electron scattering factors)
///
/// This functions creates a density representation of the entity provided by
/// the user in a density map, in Fourier space using the correct scattering
/// factors for the elements involved.
///
/// The user can also choose if the density map should be cleared of its
/// previous content before creating the density representation.
///
/// The density is generated in Fourier space. In order to
/// avoid artifacts in the final density representation, the function avoids
/// sharp frequency cutoffs by applying a Gaussian falloff. The user must
/// provide the resolutions at which the cutoff should begin and end, as opposed
/// to a single resolution cutoff value.
///
/// This function will only create a density represenation of the entities
/// (or portion of entities ) that fall within the borders of the map.
/// The user must take care that this condition is verified for all
/// entities for which he wants a representation.
///
void DLLEXPORT_OST_MOL_ALG EntityToDensityScattering(const mol::EntityView& entity_view,
                                        img::MapHandle& map,
                                        Real falloff_start,
                                        Real falloff_end,
                                        bool clear_map_flag = false,
                                        Real source_wavelength = 1.5418);

/// \brief create a density representation of an entity in a density map
///
/// This function creates a density representation of the entity provided by
/// the user in a density map, also provided by the user. The user can choose
/// the type of density of the output map:
///
/// ROSETTA_HIGH_RESOLUTION gaussian spheres in real space to represent
///                         density, one per atom, see Dimaio et al.,
///                         Refinement of Protein Structures into Low-Resolution
///                         Density Maps Using Rosetta. Journal of Molecular
///                         Biology (2009) pp. 1-10
/// ROSETTA_LOW_RESOLUTION guassian spheres in real space to represent
///                        density, one per residue. See reference above. Only
///                        useful at low resolution
///
/// The user can also choose if the density map should be cleared of its
/// previous content before creating the density representation.
///
/// The user must also provide a resolution parameter.
///
/// This function will only create a density represenation of the entities
/// (or portion of entities ) that fall within the borders of the map.
/// The user must take care that this condition is verified for all
/// entities for which he wants a representation.
///
void DLLEXPORT_OST_MOL_ALG EntityToDensityRosetta(const mol::EntityView& entity_view,
                                        img::MapHandle& map,
                                        const DensityType& density_type,
                                        Real resolution,
                                        bool clear_map_flag = false,
                                        Real source_wavelength = 1.5418);
}}} // ns

#endif // OST_ENTITY_TO_DENSITY
