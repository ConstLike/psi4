list(APPEND sources
  cuhf.cc
  frac.cc
  hf.cc
  mom.cc
  reks.cc
  reks_arch.cc
  reks_cassette.cc
  reks_rdm.cc
  reks_scf.cc
  reks_si.cc
  reks_naturals.cc
  reks_report.cc
  reks_objectives.cc
  reks_reel.cc
  reks_fon_relax.cc
  reks_si_types.cc
  rhf.cc
  rohf.cc
  sad.cc
  uhf.cc
  )
# REKS catalog blobs are runtime data (offset-ABI reks_catalog/*.bin loaded by make_catalog
# from the data dir), no longer compiled C++. Install every blob to the data dir; a new
# variant needs only its .bin dropped into reks_catalog/ -- no CMake edit (glob install).
# Use CMAKE_INSTALL_DATADIR (win-conda remaps "share"). Keep the blobs in THIS source dir,
# NOT under share/psi4/ source, whose bulk install excludes *.bin (psi4/CMakeLists.txt).
install(DIRECTORY reks_catalog/
        DESTINATION ${CMAKE_INSTALL_DATADIR}/psi4/reks_catalog
        FILES_MATCHING PATTERN "*.bin")

# Single-source REKS catalog-blob layout hash: gen_layout_hash.py derives
# REKS_LAYOUT_HASH from reks_layout.json; reks_arch.cc includes the generated header
# and the ssr2 generator stamps the same value into each .bin (fail-closed at load).
set(_reks_layout_hash_h ${CMAKE_CURRENT_BINARY_DIR}/reks_layout_hash.h)
add_custom_command(
  OUTPUT ${_reks_layout_hash_h}
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_layout_hash.py --header
          ${CMAKE_CURRENT_SOURCE_DIR}/reks_catalog/reks_layout.json ${_reks_layout_hash_h}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gen_layout_hash.py
          ${CMAKE_CURRENT_SOURCE_DIR}/reks_catalog/reks_layout.json
  COMMENT "Generating reks_layout_hash.h from reks_layout.json"
  VERBATIM)
add_custom_target(reks_layout_hash DEPENDS ${_reks_layout_hash_h})

psi4_add_module(lib scf_solver sources)
add_dependencies(scf_solver reks_layout_hash)
target_include_directories(scf_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_GTFOCK)
    add_dependencies(scf_solver INTERFACE GTFockInt GTFockCInt)
endif()

target_link_libraries(scf_solver PUBLIC diis)

if(TARGET ECPINT::ecpint)
  target_compile_definitions(scf_solver
    PUBLIC
      USING_ecpint
    )
  target_link_libraries(scf_solver
    PUBLIC
      ECPINT::ecpint
    )
endif()

if(TARGET gauxc::gauxc)
  target_include_directories(scf_solver
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()

if(TARGET BrianQC::static_wrapper)
  target_compile_definitions(scf_solver
    PUBLIC
      USING_BrianQC
    )
  target_include_directories(scf_solver
    PUBLIC
      $<TARGET_PROPERTY:BrianQC::static_wrapper,INTERFACE_INCLUDE_DIRECTORIES>
    )
endif()

if(TARGET OpenOrbitalOptimizer::OpenOrbitalOptimizer)
  target_compile_definitions(scf_solver
    PUBLIC
      USING_OpenOrbitalOptimizer
    )
  target_link_libraries(scf_solver
    PUBLIC
      OpenOrbitalOptimizer::OpenOrbitalOptimizer
    )
endif()
