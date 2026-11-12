include(TestingMacros)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../reks-h4-fon-regimes/h4_closed.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_regression_test(reks-guess "psi;quicktests;reks")
