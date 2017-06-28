# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (miniGen miniGen.cxx miniClass.cxx)
add_executable (miniSVD testSVD.cxx miniSVD.cxx)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
include_directories ( ${CMAKE_SOURCE_DIR})
target_link_libraries(miniGen ${ROOT_LIBRARIES})
target_link_libraries(miniSVD ${ROOT_LIBRARIES})

#ROOT_GENERATE_DICTIONARY(H__vec mymap.h LINKDEF mymapLinkDef.h)

