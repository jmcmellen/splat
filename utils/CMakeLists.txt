add_executable(bearing bearing.c)
target_link_libraries(bearing m)

add_executable(citydecoder citydecoder.c)

add_executable(fontdata fontdata.c)
target_link_libraries(fontdata z)

add_executable(srtm2sdf srtm2sdf.c)
target_link_libraries(srtm2sdf bz2)

add_executable(usgs2sdf usgs2sdf.c)

add_custom_target(srtm2sdf-hd ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink srtm2sdf srtm2sdf-hd
                  DEPENDS srtm2sdf
                  COMMENT "Generating srtm2sdf-hd symbolic link")

install(TARGETS bearing DESTINATION bin)
install(TARGETS citydecoder DESTINATION bin)
install(TARGETS fontdata DESTINATION bin)
install(TARGETS srtm2sdf DESTINATION bin)
install(TARGETS usgs2sdf DESTINATION bin)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/srtm2sdf-hd" DESTINATION bin)
