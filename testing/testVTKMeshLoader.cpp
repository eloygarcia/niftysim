// =========================================================================
// File:       testVTKMeshLoader.cpp
// Purpose:    Unit test for VTK mesh loader
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

#if defined _Visualisation_ && defined __GNUC__

#include "tledUnitTest.h"
#include "tledModel.h"
#include "tledMeshTopology.h"

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkPoints.h>

#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;

static bool _ExportVTK(const std::string &dstPath, const float (*nodes)[3], const int numNodes, const int (*els)[4], const int numEls) {
  vtkSmartPointer<vtkUnstructuredGrid> sp_vtkMesh;
  vtkSmartPointer<vtkUnstructuredGridWriter> sp_vtkMeshWriter;
  vtkSmartPointer<vtkTetra> sp_tetra;
  int const *pc_el;
  float const *pc_nodeCds;

  sp_tetra = vtkSmartPointer<vtkTetra>::New();
  sp_vtkMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  sp_vtkMesh->SetPoints(vtkSmartPointer<vtkPoints>::New());
  for (pc_nodeCds = &nodes[0][0]; pc_nodeCds < &nodes[numNodes][0]; pc_nodeCds += 3) sp_vtkMesh->GetPoints()->InsertNextPoint(pc_nodeCds);
  assert(sp_vtkMesh->GetNumberOfPoints() == numNodes);

  for (pc_el = &els[0][0]; pc_el < &els[numEls][0]; pc_el += 4) {
    vtkIdType vtkElNodeInds[4];

    copy(pc_el, pc_el + 4, vtkElNodeInds);
    sp_vtkMesh->InsertNextCell(sp_tetra->GetCellType(), 4, vtkElNodeInds);
  }
  assert(sp_vtkMesh->GetNumberOfCells() == numEls);

  sp_vtkMeshWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  sp_vtkMeshWriter->SetFileName(dstPath.c_str());
  tledVTK6CompatSetInput(sp_vtkMeshWriter, sp_vtkMesh);
  sp_vtkMeshWriter->Update();
     
  return sp_vtkMeshWriter->GetErrorCode() == 0;
}

static void _TestVTKLoader(const std::string &meshPath) {
  static const int maxNumEls = (int)2e3;
  static const int maxNumNodes = (int)1e3;
  const std::string vtkMeshPath = tledUnitTest::MakeTemporaryFilePath("vtkMeshLoaderTest", ".vtk");

  bool wasSuccess;
  int els[maxNumEls][4], numEls, numNodes, numSurfEls;
  float nodes[maxNumNodes][3];
  tledUnitTestXMLWriter stdXMLWriter, vtkXMLWriter;
  
  assert(tledUnitTest::ParseMSH(NULL, numNodes, NULL, numEls, NULL, numSurfEls, meshPath));
  assert(numNodes <= maxNumNodes && numEls <= maxNumEls);

  wasSuccess = tledUnitTest::ParseMSH(nodes, numNodes, els, numEls, NULL, numSurfEls, meshPath);
  wasSuccess &= _ExportVTK(vtkMeshPath, nodes, numNodes, els, numEls);

  wasSuccess &= stdXMLWriter.StartXML();
  wasSuccess &= stdXMLWriter.WriteMesh(nodes, numNodes, els, numEls);
  wasSuccess &= stdXMLWriter.CloseXML();

  wasSuccess &= vtkXMLWriter.StartXML();
  vtkXMLWriter.GetFileStream() << "\t<VTKMesh Type=\"T4\">\n"
			       << "\t\t" << vtkMeshPath << std::endl
			       << "\t</VTKMesh>\n";
  wasSuccess &= vtkXMLWriter.CloseXML();
  assert(wasSuccess);

  if (wasSuccess) {
    tledModel stdModel(const_cast<char*>(stdXMLWriter.GetFilePath().c_str()));  
    tledMesh &r_stdMesh = *stdModel.GetMesh();
    tledModel vtkModel(const_cast<char*>(vtkXMLWriter.GetFilePath().c_str()));  
    tledMesh &r_vtkMesh = *vtkModel.GetMesh();
    
    tledUnitTestAssert(string(r_vtkMesh.GetElType()) == "T4");
    tledUnitTestAssert(r_stdMesh.GetNumNodes() == r_vtkMesh.GetNumNodes());
    tledUnitTestAssert(r_stdMesh.GetNumEls() == r_vtkMesh.GetNumEls());
    tledUnitTestAssert(r_stdMesh.GetNodesPerEl() == r_vtkMesh.GetNodesPerEl());
    tledUnitTestAssert(std::string(r_stdMesh.GetElType()) == std::string(r_vtkMesh.GetElType()));
    tledUnitTestAssert(std::equal(r_stdMesh.GetAllElNodeInds(), r_stdMesh.GetAllElNodeInds() + r_stdMesh.GetNumEls()*r_stdMesh.GetNodesPerEl(), r_vtkMesh.GetAllElNodeInds()));
    tledUnitTestAssert(std::equal(r_stdMesh.GetAllNodeCds(), r_stdMesh.GetAllNodeCds() + r_stdMesh.GetNumNodes()*3, r_vtkMesh.GetAllNodeCds()));
  } else {
    std::cerr << "Warning unit test is buggy\n";
  }

  stdXMLWriter.CleanUp();
  vtkXMLWriter.CleanUp();
  tledUnitTest::RemoveFile(vtkMeshPath);
} /* _TestVTKLoader */

static void _TestTransform(const std::string &meshPath) {
  using namespace tledVectorArithmetic;

  static const int maxNumEls = (int)2e3;
  static const int maxNumNodes = (int)1e3;
  static const float axisRots[3] = {45.f, 22.5f, 17.5f};
  static const float translation[3][3] = {{0.f, 0.f, 0.f}, {1.f, 2.f, 3.f}, {99.f, 0.25f, 39.12f}};
  const std::string vtkMeshPath = tledUnitTest::MakeTemporaryFilePath("vtkMeshLoaderTest", ".vtk");  

  bool wasSuccess;
  int els[maxNumEls][4], numEls, numNodes, numSurfEls;
  float nodes[maxNumNodes][3], cnt[3];
  
  assert(tledUnitTest::ParseMSH(NULL, numNodes, NULL, numEls, NULL, numSurfEls, meshPath));
  assert(numNodes <= maxNumNodes && numEls <= maxNumEls);

  wasSuccess = tledUnitTest::ParseMSH(nodes, numNodes, els, numEls, NULL, numSurfEls, meshPath);
  wasSuccess &= _ExportVTK(vtkMeshPath, nodes, numNodes, els, numEls);
  assert(wasSuccess);

  std::fill(cnt, cnt + 3, float(0));
  for (int n = 0; n < numNodes; n++) Add(cnt, cnt, nodes[n]);
  ScalarDiv(cnt, float(numNodes));

  for (int t = 0; t < 3; t++) {
    tledUnitTestXMLWriter vtkXMLWriter;

    wasSuccess = vtkXMLWriter.StartXML();
    vtkXMLWriter.GetFileStream() << "\t<VTKMesh Type=\"T4\">\n"
				 << "\t\t" << vtkMeshPath << std::endl;
    if (Norm(translation[t]) > 0.f) {
      vtkXMLWriter.GetFileStream() << "\t\t<Translation>\n"
				   << "\t\t\t" << translation[t][0] << " " << translation[t][1] << " " << translation[t][2] << std::endl
				   << "\t\t</Translation>\n";
    }

    if (std::fabs(axisRots[t]) > 0.f) {
      vtkXMLWriter.GetFileStream() << "\t\t<Rotation>\n"
				   << cnt[0] << " " << cnt[1] << " " << cnt[2];
      for (int c = 0; c < t; c++) vtkXMLWriter.GetFileStream() << " 0";
      vtkXMLWriter.GetFileStream() << " " << axisRots[t];
      for (int c = t + 1; c < 3; c++) vtkXMLWriter.GetFileStream() << " 0";
      vtkXMLWriter.GetFileStream() << "\n\t\t</Rotation>\n";
    } 

    vtkXMLWriter.GetFileStream() << "\t</VTKMesh>\n";
    wasSuccess &= vtkXMLWriter.CloseXML();

    if (wasSuccess) {
      tledModel vtkModel(const_cast<char*>(vtkXMLWriter.GetFilePath().c_str()));  
      tledMesh &r_vtkMesh = *vtkModel.GetMesh();

      for (int n = 0; n < numNodes; n++) {
	float ref[3], test[3];

	Sub(ref, nodes[n], cnt);
	Sub(test, Sub(test, r_vtkMesh.GetAllNodeCds() + 3*n, cnt), translation[t]);
	ref[t] = test[t] = 0.f;
	assert(axisRots[t] >= 0.f);
	tledUnitTestAssert(std::fabs(ComputeAngle(test, ref)/tledPi*180.f - axisRots[t]) < 1e-2f*axisRots[t]);
      }
    } else {
      std::cerr << "Warning unit test is buggy, error creating transform " << t << std::endl;
    }
    vtkXMLWriter.CleanUp();
  }

  tledUnitTest::RemoveFile(vtkMeshPath);
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestVTKLoader(tledUnitTest::GetMeshPath("collision_bar.msh"));
  _TestVTKLoader(tledUnitTest::GetMeshPath("organic_shape.msh"));
  _TestVTKLoader(tledUnitTest::GetMeshPath("sphere.msh"));

  _TestTransform(tledUnitTest::GetMeshPath("collision_bar.msh"));
  _TestTransform(tledUnitTest::GetMeshPath("organic_shape.msh"));
  _TestTransform(tledUnitTest::GetMeshPath("sphere.msh"));

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */

#else
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << "Test requires VTK and Unix.\n";

  return EXIT_SUCCESS;
}
#endif
