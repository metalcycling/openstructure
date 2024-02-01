set(OST_GEOM_HEADERS
def.hh
exc.hh
geom.hh
mat.hh
mat2.hh
mat3.hh
mat4.hh
vec.hh
vec2.hh
vec3.hh
vec4.hh
vecmat_op.hh
vecmat2_op.hh
vecmat3_op.hh
vecmat4_op.hh
circular_iterator.hh
constants.hh
composite.hh
composite2.hh
composite3.hh
composite_op.hh
composite2_op.hh
composite3_op.hh
aligned_cuboid.hh
vec_mat_predicates.hh
quat.hh
transform.hh
point_cloud.hh
module_config.hh
)

set(OST_GEOM_SOURCES
mat2.cc
mat3.cc
mat4.cc
vec3.cc
vecmat2_op.cc
vecmat3_op.cc
vecmat4_op.cc
composite2.cc
composite3.cc
composite2_op.cc
composite3_op.cc
quat.cc
transform.cc
aligned_cuboid.cc
point_cloud.cc
)

set(OST_GEOM_EXPORT_HELPERS
vector.hh
)

module(NAME geom
       HEADERS ${OST_GEOM_EXPORT_HELPERS} IN_DIR export_helper ${OST_GEOM_HEADERS}
       SOURCES ${OST_GEOM_SOURCES}
       LINK ${BOOST_LIBRARIES})
add_dependencies(ost_geom ost_config)

