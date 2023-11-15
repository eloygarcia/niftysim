// =========================================================================
// File:       tledUnitTest.h
// Purpose:    Unit test module
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifndef __UNITTEST_HXX
#define __UNITTEST_HXX
/*----------------------------------------------------------------------------------------------------*/
#include "tledHelper.h"
#include "tledMesh.h"
#include "tledTestResources.h"

#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <cassert>

#ifdef WIN32
#define drand48() (((float)rand())/((float)RAND_MAX))
#endif

#ifdef _WIN32
#define PATHSEP "\\"
#else 
#define PATHSEP "/"
#endif

#ifndef TLED_TEST_DATA_DIR
#define TLED_TEST_DATA_DIR #TLED_TEST_DATA_DIR
#endif

/**
 * \defgroup unittesthelpers General-purpose unit test routines and types
 */

/**
 * \ingroup unittesthelpers
 */
namespace tledUnitTest {
  void InitUnitTest(void);

  /** Returns unit test mesh-file paths */
  __tled_inline std::string GetMeshPath(const std::string &meshName) {
    return TLED_TEST_DATA_DIR PATHSEP + meshName;
  }

  /** Getter for unit test resource paths */
  __tled_inline std::string GetResourcePath(const std::string &resFile) {
    return TLED_TEST_DATA_DIR PATHSEP + resFile;
  }

  tledMesh LoadMSHMesh(const std::string &path, const char *type);

  /** Generates a random temporary file name */
  std::string MakeTemporaryFilePath(const std::string &prefix, const std::string &suffix);
  
  /** Removes a file after test is done */
  void RemoveFile(const std::string &path);

  /**
   * \brief Custom definition of assert
   */
#define tledUnitTestAssert(COND) \
   {\
      if(!(COND))\
      {\
         std::cerr << "Test failed " << __FILE__ << ":" << __LINE__;\
	 std::cerr << " in function " << __FUNCTION__ << std::endl;\
	 std::cerr << "Condition: " << #COND << std::endl;		   \
	 std::exit(EXIT_FAILURE);			     \
      }\
   }

#define tledUnitTestPrintSuccess \
  { \
    std::cout << "Test " << __FILE__ << " succeeded\nExiting\n"; \
  }

#define tledUnitTestDisabled \
  { \
    std::cout << "Test deactivated\n"; \
    return EXIT_SUCCESS; \
  }

  /**
   * \brief Parses a MSH file (tetrahedral meshes only).
   *
   *
   * If called with nodesDst == NULL, elsDst == NULL, or surfEls == NULL, the number of required buffer elements will still
   * be returned in the corresponding num* output argument, however, no writing is done.
   */
  bool ParseMSH(float (*p_nodesDst)[3], int &r_numNodes, int (*p_elsDst)[4], int &r_numEls, int (*p_surfEls)[3], int &r_numSurfEls, const std::string &path);
}

/**
 * \brief Simple NiftySim XML writer class
 * \ingroup unittesthelpers
 *
 *
 * Only supports tetrahedral meshes.
 */
class tledUnitTestXMLWriter {
  /**
   * \name File I/O
   * @{
   */
private:
  std::ofstream m_FileOut;
  std::string m_FilePath;

public:
  /**
   * \brief Returns a reference to the output stream. 
   *
   *
   * Not useable until StartXML has been called, and after CloseXML has been called. 
   */
  std::ofstream& GetFileStream(void) { return m_FileOut; }

  const std::string& GetFilePath(void) const { return m_FilePath; }

  /** Removes the file from the F/S */
  void CleanUp(void);
  /** @} */

  /**
   * \name File Constituent Writers
   *
   *
   * All of these routines return true on success, false otherwise.
   * @{
   */
public:
  /** Writes the XML header */
  bool StartXML(void);

  /** Writes the mesh section of the XML */
  bool WriteMesh(const float nodes[][3], const int numNodes, const int els[][4], const int numEls);

  /** Closes the XML. Has to be called before the file can be parsed */
  bool CloseXML(void);
  /** @} */

  /**
   * \name Construction
   * @{
   */
public:
  /**
   * \brief Creates a file at the specified position.
   */
  tledUnitTestXMLWriter(const char outputPath[]);

  /**
   * \brief Creates a temporary file for output.
   */
  tledUnitTestXMLWriter(void);
  /** @} */
};
/*----------------------------------------------------------------------------------------------------*/
#endif




