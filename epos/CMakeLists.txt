add_library(epos
	epos.cpp epos.hpp
	histogram.cpp histogram.hpp
	peak_detector_CeBr.cpp peak_detector_CeBr.hpp
    impulse_fun.cpp impulse_fun.hpp
	)

target_link_libraries(epos nric)