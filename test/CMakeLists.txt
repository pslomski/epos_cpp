include_directories (
    ${CMAKE_SOURCE_DIR}/epos
    ${CMAKE_SOURCE_DIR}/step1)

add_executable(epos_test
    epos_test.cpp
    test_epos.cpp
	test_fit.cpp
    test_channel.cpp
    test_data_set_bin.cpp
    test_data_set_txt.cpp
    test_histogram.cpp
    test_impulse_fun.cpp
    test_peak_detector.cpp
    test_step1.cpp
    utils.cpp utils.hpp
    ${CMAKE_SOURCE_DIR}/step1/step1.cpp
	)

target_link_libraries(epos_test epos nric)