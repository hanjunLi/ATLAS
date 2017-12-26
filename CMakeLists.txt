cmake_minimum_required(VERSION 2.8)
project(MPCHonestMajorityNoTriples)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -maes -mpclmul -mbmi2" )

SET (BOOST_ROOT $ENV{HOME}/boost_1_64_0/)
SET (BOOST_INCLUDEDIR $ENV{HOME}/boost_1_64_0/)
SET (BOOST_LIBRARYDIR $ENV{HOME}/boost_1_64_0/stage/lib)
find_package(Boost COMPONENTS system thread REQUIRED)

INCLUDE_DIRECTORIES($ENV{HOME}/ntl-9.10.0/include $ENV{HOME} ${BOOST_INCLUDEDIR}
        $ENV{HOME}/libscapi/lib/OTExtensionBristol $ENV{HOME}/libscapi/lib/)

set(SOURCE_FILES main.cpp
        ProtocolParty.cpp ProtocolParty.h ProtocolTimer.cpp ProtocolTimer.h
        ZpKaratsubaElement.cpp ZpKaratsubaElement.h
        HashEncrypt.cpp HashEncrypt.h HashAbstract.h TemplateField.cpp)


link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib ${BOOST_LIBRARYDIR})


add_executable(MPCHonestMajorityNoTriples ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(MPCHonestMajorityNoTriples $ENV{HOME}/ntl-9.10.0/src/ntl.a
        $ENV{HOME}/libscapi/libscapi.a gmp gmpxx
        OTExtensionBristol
        $ENV{HOME}/libscapi/install/lib/libjsoncpp.a
        $ENV{HOME}/libscapi/install/lib/libsimpleot.a
        boost_system boost_thread pthread crypto dl ssl z stdc++fs)