cmake_minimum_required(VERSION 2.8)
project(FastMPC)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -maes -mpclmul -mbmi2")

INCLUDE_DIRECTORIES($ENV{HOME}/libscapi/install/include $ENV{HOME}
  $ENV{HOME}/libscapi/install/include/OTExtensionBristol)

set(SOURCE_FILES main.cpp
        ProtocolParty.cpp ProtocolParty.h ProtocolTimer.cpp ProtocolTimer.h
        ZpKaratsubaElement.cpp ZpKaratsubaElement.h
        TemplateField.h TemplateField.cpp
        Interpolate.h Interpolate.cpp Communication.h Communication.cpp
        Dispute.h Dispute.cpp)


link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib )


add_executable(MPCHonestMajorityNoTriples ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(MPCHonestMajorityNoTriples $ENV{HOME}/libscapi/build/NTL/src/ntl.a
        $ENV{HOME}/libscapi/libscapi.a gmp gmpxx
        libboost_system.a libboost_thread.a pthread crypto dl ssl z)
