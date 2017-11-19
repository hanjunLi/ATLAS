cmake_minimum_required(VERSION 2.8)
project(MPCHonestMajorityNoTriples)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -maes -mpclmul -mbmi2")

INCLUDE_DIRECTORIES($ENV{HOME}/ntl-9.10.0/include $ENV{HOME} $ENV{HOME}/libscapi/install/include
        $ENV{HOME}/libscapi/install/include/OTExtensionBristol)

set(SOURCE_FILES main.cpp
        ProtocolParty.cpp ProtocolParty.h ProtocolTimer.cpp ProtocolTimer.h
        ZpKaratsubaElement.cpp ZpKaratsubaElement.h
        HashEncrypt.cpp HashEncrypt.h HashAbstract.h TemplateField.cpp)


link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib )


add_executable(MPCHonestMajorityNoTriples ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(MPCHonestMajorityNoTriples $ENV{HOME}/ntl-9.10.0/src/ntl.a
        $ENV{HOME}/libscapi/libscapi.a gmp gmpxx
        libjsoncpp.a OTExtensionBristol libsimpleot.a
        libboost_system.a libboost_thread.a pthread crypto dl ssl z)