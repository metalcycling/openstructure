set(OST_GFX_HEADERS
render_pass.hh
bitmap_io.hh
collada_exporter.hh
color.hh
entity.hh
entity_fw.hh
gl_helper.hh
gfx.hh
gfx_fw.hh
gfx_node.hh
gfx_node_fw.hh
gfx_node_visitor.hh
gfx_object_base.hh
gfx_object.hh
gfx_object_fw.hh
gfx_prim.hh
symmetry_node.hh
gfx_test_object.hh
glwin_base.hh
input.hh
prim_list.hh
scene.hh
selection.hh
surface.hh
texture.hh
vertex_array.hh
vertex_array_helper.hh
scene_observer.hh
gl_include.hh
glext_include.hh
render_mode.hh
material.hh
gradient.hh
gradient_manager.hh
module_config.hh
primitives.hh
povray_fw.hh
povray.hh
exporter.hh
exporter_fw.hh
gost_exporter.hh
color_ops/map_handle_color_op.hh
)

set(OST_GFX_COLOR_OPS_HEADERS
color_op.hh
by_element_color_op.hh
by_chain_color_op.hh
uniform_color_op.hh
gradient_color_op.hh
entity_view_color_op.hh
basic_gradient_color_op.hh
gradient_level_color_op.hh
map_handle_color_op.hh
)

set(OST_GFX_RENDER_OPTIONS_HEADERS
render_options.hh
line_render_options.hh
custom_render_options.hh
cpk_render_options.hh
simple_render_options.hh
sline_render_options.hh
line_trace_render_options.hh
cartoon_render_options.hh
trace_render_options.hh
)


set(OST_GFX_IMPL_HEADERS
cartoon_renderer.hh
custom_renderer.hh
cpk_renderer.hh
debug_renderer.hh
simple_renderer.hh
sline_renderer.hh
trace_renderer.hh
line_trace_renderer.hh
backbone_trace.hh
entity_detail.hh
entity_renderer.hh
trace_renderer_base.hh
connect_renderer_base.hh
mapped_property.hh
entity_renderer_fw.hh
tabulated_trig.hh
fast_spheres.hh
)

set(OST_GFX_SOURCES
bitmap_io.cc
exporter.cc
collada_exporter.cc
color.cc
primitives.cc
entity.cc
symmetry_node.cc
gfx_node.cc
gfx_object.cc
gfx_prim.cc
gfx_test_object.cc
input.cc
prim_list.cc
scene.cc
selection.cc
surface.cc
gradient.cc
gradient_manager.cc
vertex_array.cc
vertex_array_helper.cc
material.cc
povray.cc
gost_exporter.cc
texture.cc
color_ops/color_op.cc
color_ops/by_element_color_op.cc
color_ops/by_chain_color_op.cc
color_ops/uniform_color_op.cc
color_ops/gradient_color_op.cc
color_ops/entity_view_color_op.cc
color_ops/basic_gradient_color_op.cc
color_ops/gradient_level_color_op.cc
impl/entity_renderer.cc
impl/entity_detail.cc
impl/cartoon_renderer.cc
impl/custom_renderer.cc
impl/cpk_renderer.cc
impl/debug_renderer.cc
impl/simple_renderer.cc
impl/sline_renderer.cc
impl/line_trace_renderer.cc
impl/backbone_trace.cc
impl/connect_renderer_base.cc
impl/trace_renderer_base.cc
impl/mapped_property.cc
impl/tabulated_trig.cc
impl/trace_renderer.cc
impl/fast_spheres.cc
render_options/render_options.cc
render_options/line_render_options.cc
render_options/custom_render_options.cc
render_options/cpk_render_options.cc
render_options/simple_render_options.cc
render_options/sline_render_options.cc
render_options/line_trace_render_options.cc
render_options/cartoon_render_options.cc
render_options/trace_render_options.cc
color_ops/map_handle_color_op.cc
)

set(OST_GFX_MAP_IMPL_HEADERS
map_iso_gen.hh
map_iso_gen_s.hh
map_iso_gen_s.hh
map_iso_spec.hh
map_octree.hh
octree_isocont.hh
)
set(OST_GFX_MAP_SOURCES
map_iso.cc
map_slab.cc
impl/map_iso_gen.cc
impl/map_iso_gen_s.cc
impl/map_iso_gen_o.cc
impl/map_octree.cc
impl/octree_isocont.cc
)
set(OST_GFX_MAP_HEADERS
map_iso.hh
map_slab.hh
map_iso_prop.hh
)  

if (USE_SHADER)
  list(APPEND OST_GFX_SOURCES shader.cc)
  list(APPEND OST_GFX_HEADERS shader.hh)
  list(APPEND OST_GFX_SOURCES impl/scene_fx.cc)
  list(APPEND OST_GFX_IMPL_HEADERS scene_fx.hh)
  if (NOT APPLE)
    set(OST_GLEW_HEADERS
       glew.h
       glxew.h
       wglew.h IN_DIR GL
    )
    list(APPEND OST_GFX_SOURCES GL/glew.c)
    include_directories("${STAGE_DIR}/include/ost/gfx")
  endif()  
endif()

set(OST_GFX_DEPENDENCIES "ost_conop;ost_seq;ost_img;ost_img_alg")

module(NAME gfx SOURCES ${OST_GFX_SOURCES} ${OST_GFX_MAP_SOURCES} 
       HEADERS ${OST_GFX_MAP_IMPL_HEADERS} ${OST_GFX_IMPL_HEADERS} IN_DIR impl
               ${OST_GFX_RENDER_OPTIONS_HEADERS} IN_DIR render_options
               ${OST_GFX_COLOR_OPS_HEADERS} IN_DIR color_ops
               ${OST_GLEW_HEADERS}
               ${OST_GFX_HEADERS} ${OST_GFX_MAP_HEADERS}
       DEPENDS_ON ${OST_GFX_DEPENDENCIES})

include_directories(${PNG_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR})

# link against OpenGL and PNG libraries
target_link_libraries(ost_gfx ${OPENGL_LIBRARIES} ${PNG_LIBRARIES})

if (USE_SHADER)
  set(SHADER_FILES
    shader/aaline_fs.glsl
    shader/amboccl_fs.glsl
    shader/basic_fs.glsl
    shader/basic_vs.glsl
    shader/beacon_fs.glsl
    shader/convolute1_fs.glsl
    shader/dumpnorm_fs.glsl
    shader/dumpnorm_vs.glsl
    shader/fast_sphere_fs.glsl
    shader/fraglight_fs.glsl
    shader/fraglight_vs.glsl
    shader/iso_fs.glsl
    shader/iso_vs.glsl
    shader/noop_vs.glsl
    shader/outline_vs.glsl
    shader/quadpp_vs.glsl
    shader/scenefx_fs.glsl
    shader/scenefx_vs.glsl
    shader/selfx_fs.glsl
    shader/selfx_vs.glsl
    shader/screenblur4_fs.glsl
    shader/test_tex_fs.glsl
    shader/material_phong.glsl
    shader/material_hemi.glsl
    shader/material_toon1.glsl
    shader/material_toon2.glsl
    shader/anaglyph_fs.glsl
  )
  copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${SHARED_DATA_PATH}/shader" "${SHADER_FILES}" 
                    "SHADER_TARGETS" ost_gfx)  
  install(FILES ${SHADER_FILES} DESTINATION "share/openstructure/shader")
endif(USE_SHADER)
set(TEXTURE_FILES
  test_texture.png
  glyph_texture.png
)
copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${SHARED_DATA_PATH}/textures" "${TEXTURE_FILES}" 
                  "TEXTURE_TARGETS" ost_gfx)  

install(FILES ${TEXTURE_FILES} DESTINATION "share/openstructure/textures")
