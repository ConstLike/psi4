include(TestingMacros)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../reks-h2-fon-regimes/h2_opening.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_regression_test(reks-dft "psi;quicktests;smoketests;reks;dft")
