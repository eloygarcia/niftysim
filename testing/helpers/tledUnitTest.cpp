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

#include "tledUnitTest.h"
#include "tledMSHMeshLoader.h"

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <cstring>

#ifdef WIN32
#define srand48(x) (srand((x)))
#endif

using namespace std;

void tledUnitTest::InitUnitTest(void) {
  srand(static_cast<unsigned int>(time(NULL)));
  srand48(rand());
} /* InitUnitTest */

static inline void _skip_line(std::ifstream &r_fin) {
  char c;

  r_fin.unsetf(std::ios_base::skipws);
  do {
    r_fin >> c;
  } while (c != '\n' && !r_fin.eof() && !r_fin.fail());

  r_fin.setf(std::ios_base::skipws);
} /* _skip_line */

bool tledUnitTest::ParseMSH(float (*p_nodes)[3], int &r_nof_nodes, int (*p_elems)[4], int &r_nof_tetras, int (*p_surfelems)[3], int &r_nof_trias, const std::string &path) {
#define __parseMSHAssert(cond) if (!(cond)) goto tetramesh_c_tetramesh_c_fail
  ifstream fin(path.c_str());
  string tmpstring;

  __parseMSHAssert(fin.is_open());

  fin >> tmpstring;
  __parseMSHAssert(tmpstring == "$MeshFormat");
  fin >> tmpstring;
  __parseMSHAssert(tmpstring == "2.1" || tmpstring == "2.2");
  _skip_line(fin);
  fin >> tmpstring;
  __parseMSHAssert(tmpstring == "$EndMeshFormat");
  fin >> tmpstring;
  __parseMSHAssert(tmpstring == "$Nodes");

  __parseMSHAssert(!(fin.fail() || fin.eof()));

  {
    int vInd, tmp_vInd;

    fin >> r_nof_nodes;
    if (p_nodes != NULL) {
      for (vInd = 0; vInd < r_nof_nodes && !fin.fail() && !fin.eof(); vInd++) {
	float *p_vtx;

	p_vtx = p_nodes[vInd];
	fin >> tmp_vInd;
	__parseMSHAssert(tmp_vInd == vInd + 1);
	fin >> p_vtx[0], fin >> p_vtx[1], fin >> p_vtx[2];      
      }    
    } else {
      for (vInd = 0; vInd < r_nof_nodes && !fin.fail() && !fin.eof(); vInd++) {
	float tmpCd;
	
	fin >> tmp_vInd;
	__parseMSHAssert(tmp_vInd == vInd + 1);
	fin >> tmpCd, fin >> tmpCd, fin >> tmpCd;      
      }	
    }

    __parseMSHAssert(vInd == r_nof_nodes);
  }

  fin >> tmpstring;
  __parseMSHAssert(tmpstring == "$EndNodes");

  fin >> tmpstring;
  __parseMSHAssert(tmpstring == "$Elements");

  {
    int nof_elems, eidx, tmp;

    /*
     * Note: MSH format allows for surface and volume elements within the same element block
     * -> need to count tetras!
     */
    fin >> nof_elems;
    r_nof_tetras = 0, r_nof_trias = 0;
    for (eidx = 0; eidx < nof_elems && !fin.eof() && !fin.fail(); eidx++) {
      int elemtype, nof_labels, lidx;

      fin >> tmp;
      assert(tmp == eidx + 1);
      fin >> elemtype;
      if (elemtype == 2) {
	int *p_tria = NULL;
	
	if (p_surfelems) p_tria = p_surfelems[r_nof_trias];
	r_nof_trias += 1;

	fin >> nof_labels;
	for (lidx = 0; lidx < nof_labels; lidx++) fin >> tmp;
	if (p_surfelems) for (lidx = 0; lidx < 3; lidx++) fin >> p_tria[lidx], p_tria[lidx] -= 1;
	else for (lidx = 0; lidx < 3; lidx++) fin >> tmp;
      } else {	
	int *p_tetra = NULL;
	
	assert(elemtype == 4);
	if (p_elems != NULL) p_tetra = p_elems[r_nof_tetras];
	r_nof_tetras += 1;
	fin >> nof_labels;
	for (lidx = 0; lidx < nof_labels; lidx++) fin >> tmp;
	if (p_elems != NULL) for (lidx = 0; lidx < 4; lidx++) fin >> p_tetra[lidx], p_tetra[lidx] -= 1;
	else for (lidx = 0; lidx < 4; lidx++) fin >> tmp;
      }
    }
    
    __parseMSHAssert(eidx == r_nof_tetras + r_nof_trias);
  }

  return true;

 tetramesh_c_tetramesh_c_fail:
  cerr << "Reading MSH file " << path << " failed.\n";
  
  return false;

#undef __parseMSHAssert
} /* ParseMSH */

tledMesh tledUnitTest::LoadMSHMesh(const std::string &path, const char *type) {
  tledMSHMeshLoader loader;
  tledMesh mesh;

  loader.SetMeshType(type);
  loader.SetOutputMesh(mesh);
  loader.SetFilename(path);
  loader.Read();

  return mesh;
}

std::string tledUnitTest::MakeTemporaryFilePath(const std::string &prefix, const std::string &suffix) {
  string path;

#if defined __GNUC__ && !defined __llvm__
  ostringstream oss;
  char buffer[512];
  
  oss << "/tmp/" << prefix << "XXXXXX" << suffix;
  assert(512 > oss.str().length());
  mkstemps(strncpy(buffer, oss.str().c_str(), 512), suffix.length());
  path = buffer;
#else
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif
  /** Simply using tmpnam... */
  path = std::string(tempnam(NULL, NULL)) + suffix;
#endif

  return path;
}

void tledUnitTest::RemoveFile(const std::string &path) {
  std::remove(path.c_str());
}

/*****************************************************************************************************************/
tledUnitTestXMLWriter::tledUnitTestXMLWriter(const char outputPath[]) : m_FilePath(outputPath) {}
tledUnitTestXMLWriter::tledUnitTestXMLWriter() {
  m_FilePath = tledUnitTest::MakeTemporaryFilePath("tledUnitTestMesh_", ".xml");
}

bool tledUnitTestXMLWriter::StartXML() {
  m_FileOut.open(GetFilePath().c_str());
  
  m_FileOut << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
	    << "<Model>\n";

  return !m_FileOut.fail();
}

bool tledUnitTestXMLWriter::WriteMesh(const float nodes[][3], const int numNodes, const int els[][4], const int numEls) {
  int i;

  m_FileOut << "\t<Nodes DOF=\"3\" NumNodes=\"" << numNodes << "\">\n";    
  for (i = 0; i < numNodes && !m_FileOut.fail(); i++) {
    m_FileOut << "\t\t" << nodes[i][0] << " " << nodes[i][1] << " " << nodes[i][2] << endl;
  }
  m_FileOut << "\t</Nodes>\n";
  
  if (i < numNodes) return false;

  m_FileOut << "\t<Elements NumEls=\"" << numEls << "\" Type=\"T4\">\n";
  for (i = 0; i < numEls  && !m_FileOut.fail(); i++) {
    m_FileOut << "\t\t" << els[i][0] << " " << els[i][1] << " " << els[i][2] << " " << els[i][3] << endl;
  }
  m_FileOut << "\t</Elements>\n";
  
  return i == numEls;
}
    
bool tledUnitTestXMLWriter::CloseXML() {
  m_FileOut << "</Model>\n";
  m_FileOut.close();

  return !m_FileOut.fail();
}

void tledUnitTestXMLWriter::CleanUp() {
  tledUnitTest::RemoveFile(GetFilePath());
}
