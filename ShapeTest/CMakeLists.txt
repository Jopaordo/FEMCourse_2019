########### next target ###############

file(GLOB HEADERS *.h)

file(GLOB SOURCE *.cpp)

add_executable(ShapeTest main.cpp ${SOURCE} ${HEADERS})

target_link_libraries(ShapeTest pz)