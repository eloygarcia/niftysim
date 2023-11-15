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
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>

#include "tledUnitTest.h"
#include "tledVTKMeshLoader.h"
#include "tledVTKMeshExporter.h"
#include "tledVectorArithmetic.h"

#include <cstdlib>
#include <algorithm>

using namespace std;

/*
 * Test requires a working VTK mesh loader implementation!
 */
static void _TestWithLoader(const string &meshPath) {
  const std::string fileNameTemp = tledUnitTest::MakeTemporaryFilePath("vtkMeshExporterTest", ".vtk");

  tledVTKMeshLoader loader;
  tledVTKMeshExporter exporter;
  tledMesh mesh, reloadMesh;

  loader.SetMeshType("T4");
  loader.SetOutputMesh(mesh);
  loader.SetFilename(meshPath);
  loader.Read();

  exporter.SetFileName(fileNameTemp);
  exporter.SetMesh(mesh);
  tledUnitTestAssert(exporter.Write());
  
  loader.SetOutputMesh(reloadMesh);
  loader.SetFilename(fileNameTemp);
  loader.Read();
  tledUnitTestAssert(reloadMesh.GetNumNodes() == mesh.GetNumNodes());
  tledUnitTestAssert(reloadMesh.GetNumEls() == mesh.GetNumEls());
  tledUnitTestAssert(string(reloadMesh.GetElType()) == mesh.GetElType());

  tledUnitTestAssert(equal(reloadMesh.GetAllElNodeInds(), reloadMesh.GetAllElNodeInds() + reloadMesh.GetNodesPerEl()*reloadMesh.GetNumEls(), mesh.GetAllElNodeInds()));

  {
    using namespace tledVectorArithmetic;

    const float *nodesStart = mesh.GetAllNodeCds();
    const float *nodesEnd = mesh.GetAllNodeCds() + 3*mesh.GetNumNodes();
    const float *rNodesStart = reloadMesh.GetAllNodeCds();

    float const *pc_node, *pc_rNode;
    float diff[3];

    for (pc_node = nodesStart, pc_rNode = rNodesStart; pc_node < nodesEnd; pc_node += 3, pc_rNode += 3) {
      tledUnitTestAssert(Norm(Sub(diff, pc_node, pc_rNode)) < 1e-3*Norm(pc_node));
    }
  }  

  tledUnitTest::RemoveFile(fileNameTemp);
}

static void _TestAttributes() {
  static const int numElems = 200;  
  const std::string fileNameTemp = tledUnitTest::MakeTemporaryFilePath("vtkMeshExporterTest", ".vtk");

  tledVTKMeshExporter exporter;
  tledMesh rMesh;
  float refScalarAttrib[4*numElems], refVecAttrib[3*4*numElems];  

  rMesh.SetNumberOfElements(numElems, "T4");
  rMesh.SetNumberOfNodes(4*numElems, 3);
  for (int eInd = 0; eInd < numElems; eInd++) {
    for (int vInd = 0; vInd < 4; vInd++) {
      const int nInd = 4*eInd+vInd;

      refScalarAttrib[nInd] = 100*(drand48() - drand48());
      for (int cInd = 0; cInd < 3; cInd++) {
	rMesh.GetAllNodeCds()[3*nInd+cInd] = 100*(drand48() - drand48());
	refVecAttrib[3*nInd+cInd] = 100*(drand48() - drand48());
      }
      rMesh.GetAllElNodeInds()[nInd] = nInd;
    }
  }

  exporter.SetFileName(fileNameTemp);
  exporter.SetMesh(rMesh);
  exporter.AddNodeVectorAttribute("vec_attrib", refVecAttrib);
  exporter.AddNodeScalarAttribute("scal_attrib", refScalarAttrib);
  exporter.Write();

  {
    vtkSmartPointer<vtkUnstructuredGridReader> sp_reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    vtkSmartPointer<vtkUnstructuredGrid> sp_grid;
    vtkFloatArray *p_array;
    
    sp_reader->SetFileName(fileNameTemp.c_str());
    sp_reader->Update();
    sp_grid = sp_reader->GetOutput();

    tledUnitTestAssert(sp_grid->GetPointData()->HasArray("vec_attrib"));
    tledUnitTestAssert(sp_grid->GetPointData()->HasArray("scal_attrib"));

    p_array = dynamic_cast<vtkFloatArray*>(sp_grid->GetPointData()->GetArray("vec_attrib"));
    tledUnitTestAssert(p_array->GetNumberOfTuples() == 4*numElems);
    tledUnitTestAssert(p_array->GetNumberOfComponents() == 3);
    for (int vInd = 0; vInd < p_array->GetNumberOfTuples(); vInd++) {
      float tmpTuple[3];

      p_array->GetTupleValue(vInd, tmpTuple);
      for (int cInd = 0; cInd < 3; cInd++) {
	tledUnitTestAssert(std::fabs(tmpTuple[cInd] - refVecAttrib[cInd+3*vInd]) < 1e-5*std::fabs(refVecAttrib[cInd+3*vInd]));
      }
    }

    p_array = dynamic_cast<vtkFloatArray*>(sp_grid->GetPointData()->GetArray("scal_attrib"));
    tledUnitTestAssert(p_array->GetNumberOfTuples() == 4*numElems);
    tledUnitTestAssert(p_array->GetNumberOfComponents() == 1);
    for (int vInd = 0; vInd < p_array->GetNumberOfTuples(); vInd++) {
      float tmpTuple;

      p_array->GetTupleValue(vInd, &tmpTuple);
      tledUnitTestAssert(std::fabs(tmpTuple - refScalarAttrib[vInd]) < 1e-5*std::fabs(refScalarAttrib[vInd]));
    }
  }

  tledUnitTest::RemoveFile(fileNameTemp);
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestWithLoader(tledUnitTest::GetMeshPath("collision_bar.vtk"));
  _TestWithLoader(tledUnitTest::GetMeshPath("organic_shape.vtk"));
  _TestWithLoader(tledUnitTest::GetMeshPath("sphere.vtk"));

  _TestAttributes();

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
