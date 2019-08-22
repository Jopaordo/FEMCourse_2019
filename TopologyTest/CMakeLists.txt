########### next target ###############

file(GLOB HEADERS *.h)

file(GLOB SOURCE *.cpp)

add_executable(TopologyTest main.cpp ${SOURCE} ${HEADERS})

target_link_libraries(TopologyTest pz)