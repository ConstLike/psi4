include(TestingMacros)

file(COPY h4_closed.json h4_opening.json h4_open.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_regression_test(reks-h4-fon-regimes "psi;quicktests;reks")
