include(TestingMacros)

file(COPY c2h4_opening.json c2h4_open.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_regression_test(reks-c2h4-si-properties "psi;quicktests;reks")
