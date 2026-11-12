include(TestingMacros)

file(COPY h2_closed.json h2_opening.json h2_open.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_regression_test(reks-h2-fon-regimes "psi;quicktests;smoketests;reks")
