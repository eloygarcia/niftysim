#cmake_minimum_required(VERSION 2.4)

#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/Source/xmlParser)

#-----------------------------------------------------------------------------
SET(NAME xmlParser)

ADD_LIBRARY(${NAME} xmlParser.cpp)
INSTALL( FILES xmlParser.h DESTINATION include )
INSTALL(TARGETS ${NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
#-----------------------------------------------------------------------------



