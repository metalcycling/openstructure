set(OST_GEOM_PYMOD_SOURCES
  export_vec2.cc
  export_vecmat2_op.cc
  export_vec3.cc
  export_vecmat3_op.cc
  export_vec4.cc
  export_vecmat4_op.cc
  wrap_geom.cc
  export_mat2.cc
  export_mat3.cc
  export_mat4.cc
  export_composite2.cc
  export_composite3.cc
  export_composite2_op.cc
  export_composite3_op.cc
  export_quat.cc
  export_transform.cc
)

if (NOT ENABLE_STATIC)
  pymod(NAME geom CPP ${OST_GEOM_PYMOD_SOURCES} PY __init__.py)
endif()
