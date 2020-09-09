add_executable(step1
    main.cpp
	step1.cpp step1.hpp
	)

target_link_libraries(step1
    epos
    nric)

include_directories(step1)
