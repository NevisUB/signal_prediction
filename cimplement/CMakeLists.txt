cmake_minimum_required (VERSION 2.6)
project (miniGen)
# The version number.
set (miniGen_VERSION_MAJOR 2)
set (miniGen_VERSION_MINOR 0)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
    else ()
	    set (CMAKE_CXX_STANDARD 11)
endif ()

#set (CMAKE_CXX_FLAGS "  -I/grid/fermiapp/products/larsoft/root/v6_08_06d/Linux64bit+2.6-2.12-e14-nu-prof/include -L/grid/fermiapp/products/larsoft/root/v6_08_06d/Linux64bit+2.6-2.12-e14-nu-prof/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic")
#execute_process(COMMAND "root-config --glibs --cflags" OUTPUT_VARIABLE CMAKE_CXX_FLAGS)


find_package(ROOT REQUIRED)# COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread)
include(${ROOT_USE_FILE})

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/)
#ROOT_GENERATE_DICTIONARY(G__Event MCEventWeight.h LINKDEF MCEventWeightLinkDef.h)

#add_subdirectory(inc)
#add_subdirectory(tinyxml)
add_subdirectory(src)
#add_subdirectory(examples)

