// =========================================================================
// File:       testMSHMeshLoader.cpp
// Purpose:    Unit test for VTK mesh loader
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledModel.h"
#include "tledMeshTopology.h"

#include <cstdlib>

using namespace std;

/**
 * \brief Runs a number of tests on the given geometry
 */
static void _TestAgainstVTKLoader(const string &mshPath, const string &vtkPath) {
  tledUnitTestXMLWriter vtkXMLWriter, mshXMLWriter;
  bool wasSuccess;

  wasSuccess = true;
  wasSuccess &= vtkXMLWriter.StartXML();
  vtkXMLWriter.GetFileStream() << "\t<VTKMesh Type=\"T4\">\n"
			       << "\t\t" << vtkPath << endl
			       << "\t</VTKMesh>\n";
  wasSuccess &= vtkXMLWriter.CloseXML();
  assert(wasSuccess);

  wasSuccess &= mshXMLWriter.StartXML();
  mshXMLWriter.GetFileStream() << "\t<MSHMesh Type=\"T4\">\n"
			       << "\t\t" << mshPath << endl
			       << "\t</MSHMesh>\n";
  wasSuccess &= mshXMLWriter.CloseXML();
  assert(wasSuccess);

  if (wasSuccess) {
    tledModel mshModel(const_cast<char*>(mshXMLWriter.GetFilePath().c_str()));  
    tledMesh &r_mshMesh = *mshModel.GetMesh();
    tledModel vtkModel(const_cast<char*>(vtkXMLWriter.GetFilePath().c_str()));  
    tledMesh &r_vtkMesh = *vtkModel.GetMesh();
    float const *pc_vtkNCd, *pc_mshNCd;
    
    tledUnitTestAssert(string(r_vtkMesh.GetElType()) == "T4");
    tledUnitTestAssert(r_mshMesh.GetNumNodes() == r_vtkMesh.GetNumNodes());
    tledUnitTestAssert(r_mshMesh.GetNumEls() == r_vtkMesh.GetNumEls());
    tledUnitTestAssert(r_mshMesh.GetNodesPerEl() == r_vtkMesh.GetNodesPerEl());
    tledUnitTestAssert(std::string(r_mshMesh.GetElType()) == std::string(r_vtkMesh.GetElType()));
    tledUnitTestAssert(std::equal(r_mshMesh.GetAllElNodeInds(), r_mshMesh.GetAllElNodeInds() + r_mshMesh.GetNumEls()*r_mshMesh.GetNodesPerEl(), r_vtkMesh.GetAllElNodeInds()));

    for (pc_vtkNCd = r_vtkMesh.GetAllNodeCds(), pc_mshNCd = r_mshMesh.GetAllNodeCds(); pc_mshNCd < r_mshMesh.GetAllNodeCds() + r_mshMesh.GetNumNodes()*3; pc_mshNCd++, pc_vtkNCd++) {
      tledUnitTestAssert(fabsf(*pc_vtkNCd - *pc_mshNCd) <= 1e-4*fabsf(*pc_mshNCd));
    }
  } else {
    cerr << "Warning unit test is buggy\n";
  }
} /* _TestAgainstVTKLoader */

template <const int t_numElNodes>
static void _TestAgainstKnown(const string &mshPath, const int numNodes, const int numEls) {
  tledUnitTestXMLWriter vtkXMLWriter, mshXMLWriter;
  bool wasSuccess;

  wasSuccess = true;
  wasSuccess &= mshXMLWriter.StartXML();
  mshXMLWriter.GetFileStream() << "\t<MSHMesh Type=\"" << (t_numElNodes == 4? "T4" : "H8") << "\">\n"
			       << "\t\t" << mshPath << endl
			       << "\t</MSHMesh>\n";
  wasSuccess &= mshXMLWriter.CloseXML();
  assert(wasSuccess);

  if (wasSuccess) {
    tledModel mshModel(const_cast<char*>(mshXMLWriter.GetFilePath().c_str()));  
    tledMesh &r_mesh = *mshModel.GetMesh();
    
    tledUnitTestAssert(r_mesh.GetNodesPerEl() == t_numElNodes);
    tledUnitTestAssert(r_mesh.GetNumNodes() == numNodes);
    tledUnitTestAssert(r_mesh.GetNumEls() == numEls);
  }
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestAgainstKnown<8>(tledUnitTest::GetMeshPath("box8.msh"), 36, 12);

  _TestAgainstVTKLoader(tledUnitTest::GetMeshPath("collision_bar.msh"), tledUnitTest::GetMeshPath("collision_bar.vtk"));
  _TestAgainstVTKLoader(tledUnitTest::GetMeshPath("organic_shape.msh"), tledUnitTest::GetMeshPath("organic_shape.vtk"));
  _TestAgainstVTKLoader(tledUnitTest::GetMeshPath("sphere.msh"), tledUnitTest::GetMeshPath("sphere.vtk"));

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
