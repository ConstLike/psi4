list(APPEND sources
  reks_arch.cc
  reks_cassette.cc
  reks_fon_relax.cc
  reks_fon_step.cc
  reks_gvb_driver.cc
  reks_naturals.cc
  reks_objectives.cc
  reks_properties.cc
  reks_rdm.cc
  reks_reel.cc
  reks_report.cc
  reks_scf.cc
  reks_si.cc
  reks_si_types.cc
  reks_trah.cc
  reks_xc.cc
  )

# REKS_LAYOUT_HASH derived from reks_layout.json; reks_arch.cc includes the generated header
# and the blob generator stamps the same value into each .bin (fail-closed at load).
set(_reks_layout_json ${PROJECT_SOURCE_DIR}/share/psi4/reks_catalog/reks_layout.json)
set(_reks_layout_hash_h ${CMAKE_CURRENT_BINARY_DIR}/reks_layout_hash.h)
add_custom_command(
  OUTPUT ${_reks_layout_hash_h}
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_layout_hash.py --header
          ${_reks_layout_json} ${_reks_layout_hash_h}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gen_layout_hash.py
          ${_reks_layout_json}
  COMMENT "Generating reks_layout_hash.h from reks_layout.json"
  VERBATIM)
add_custom_target(reks_layout_hash DEPENDS ${_reks_layout_hash_h})

psi4_add_module(lib reks sources)
add_dependencies(reks reks_layout_hash)
target_include_directories(reks PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

if (psi4_ENABLE_PRECOMPILE_HEADERS)
    target_precompile_headers(reks REUSE_FROM psio)
endif()
