// =========================================================================
// File:       testXMLDeformableContactSurfaceCreator
// Purpose:    Unit test for VTK mesh loader
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledDeformableContactSurfaceXMLExporter.h"
#include "tledXMLDeformableContactSurfaceCreator.h"
#include "tledModel.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"

#ifdef GPU_GP_CONTACT
#include "tledCUDAUnitTest.h"
#include "tledDeformableContactSurfaceGPU.h"
#endif

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>

template <const bool t_isCPU>
struct _SpecialTests {
  template <class TSurface>
  void RunTests(TSurface &r_surface, TSurface &output) {}
};

template <>
struct _SpecialTests<true> {
  template <class TSurface>
  void RunTests(TSurface &r_surface, TSurface &output) {
    using namespace tledVectorArithmetic;

    for (int n = 0; n < r_surface.GetNumberOfNodes(); n++) {
      float d[3];

      tledUnitTestAssert(Norm(Sub(d, r_surface.GetOldNodeCoordinates(n), output.GetOldNodeCoordinates(n))) <= 1e-3*Norm(r_surface.GetOldNodeCoordinates(n)));
      tledUnitTestAssert(Norm(Sub(d, r_surface.GetOldNodeCoordinates(n), output.GetOldNodeCoordinates(n))) <= 1e-3*Norm(r_surface.GetOldNodeCoordinates(n)));
    }
  }
};

template <>
struct _SpecialTests<false> {
  template <class TSurface>
  void RunTests(TSurface &r_surface, TSurface &output) {
  }
};

template <const bool t_isCPU, class TSurface>
static void _TestWriteRead(TSurface &r_surface) {
  using namespace tledVectorArithmetic;

  TSurface *p_output;

  if (dynamic_cast<tledDeformableMembraneContactSurface*>(&r_surface) != NULL) {
    p_output = static_cast<TSurface*>(tledDeformableMembraneContactSurface::CreateSurface(r_surface.GetNumberOfFacetVertices() == 3? "T3" : "Q4", !t_isCPU));
  } else {
    p_output = static_cast<TSurface*>(tledDeformableContactSurface::CreateSurface(r_surface.GetNumberOfFacetVertices() == 3? "T3" : "Q4", !t_isCPU));
  }

  {
    tledDeformableContactSurfaceXMLExporter<TSurface> exporter;
    tledXMLDeformableContactSurfaceCreator<TSurface> importer;

    exporter.SetInput(r_surface);
    exporter.Export();
    importer.SetXMLRoot(exporter.GetRootNode());
    importer.SetOutputMesh(*p_output);
    importer.Create();

    exporter.GetRootNode().deleteNodeContent();
  }

  tledUnitTestAssert(p_output->GetNumberOfNodes() == r_surface.GetNumberOfNodes());
  tledUnitTestAssert(p_output->GetNumberOfFacets() == r_surface.GetNumberOfFacets());
  tledUnitTestAssert(p_output->GetUpdateCount() == 0 && p_output->GetSaveCount() == 0);
  
  for (int n = 0; n < r_surface.GetNumberOfNodes(); n++) {
    float d[3];

    tledUnitTestAssert(Norm(Sub(d, r_surface.GetNodeCoordinates(n), p_output->GetNodeCoordinates(n))) <= 1e-3*Norm(r_surface.GetNodeCoordinates(n)));
    tledUnitTestAssert(r_surface.GetNumberOfNodeFacets(n) == p_output->GetNumberOfNodeFacets(n));
    tledUnitTestAssert(std::equal(r_surface.GetNodeFacetIndices(n), r_surface.GetNodeFacetIndices(n) + r_surface.GetNumberOfNodeFacets(n), p_output->GetNodeFacetIndices(n)));
  }

  for (int f = 0; f < r_surface.GetNumberOfFacets(); f++) {
    typedef typename TSurface::Facet __Facet;

    tledUnitTestAssert(std::equal(r_surface.GetFacet(f).NodeIndices, r_surface.GetFacet(f).NodeIndices + __Facet::NumberOfVertices, p_output->GetFacet(f).NodeIndices));
    tledUnitTestAssert(std::equal(r_surface.GetFacet(f).EdgeIndices, r_surface.GetFacet(f).EdgeIndices + __Facet::NumberOfVertices, p_output->GetFacet(f).EdgeIndices));
  }

  _SpecialTests<t_isCPU>().RunTests(r_surface, *p_output);
}

template <const bool t_doCPUTests, class TSurfaceT3, class TMembraneSurfaceT3>
static void _TestXMLExportImport(const std::string &inSimFile) {
  const tledModel model(inSimFile.c_str());
  const tledMesh &mesh = *model.GetMesh();

  if (mesh.GetElType() == std::string("T4") || mesh.GetElType() == std::string("T4ANP")) {
    if (model.GetNumberOfShellElementSets() == 0 || model.DoShellUseMeshSurface()) {
      TSurfaceT3 *p_surface = static_cast<TSurfaceT3*>(tledDeformableContactSurface::CreateSurface(mesh, !t_doCPUTests));
      
      p_surface->Init();
      _TestWriteRead<t_doCPUTests, TSurfaceT3>(*p_surface);
      delete p_surface;
    } else {
      tledSurface *p_membrane = model.GetGenericShellMesh();
      std::vector<float> t(p_membrane->GetNumberOfFacets(), 0);
      TMembraneSurfaceT3 *p_surface = static_cast<TMembraneSurfaceT3*>(tledDeformableMembraneContactSurface::CreateSurface(mesh, *p_membrane, !t_doCPUTests));
      
      p_surface->SetMembraneFacetThickenesses(&t.front());
      p_surface->Init();
      _TestWriteRead<t_doCPUTests, TMembraneSurfaceT3>(*p_surface);

      delete p_membrane;
      delete p_surface;
    }
  } else {
    tledFatalError("Incompatible contact surface.");
  }
}

static void _TestXMLExportImportCPU(const std::string &inSimFile) {
  _TestXMLExportImport<true, tledDeformableContactSurfaceT3CPU, tledDeformableMembraneContactSurfaceT3CPU>(inSimFile);
}

#ifdef GPU_GP_CONTACT
static void _TestXMLExportImportGPU(const std::string &inSimFile) {
  _TestXMLExportImport<false, tledDeformableContactSurfaceT3GPU, tledDeformableMembraneContactSurfaceT3CPU>(inSimFile);
}
#endif

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestXMLExportImportCPU(tledUnitTest::GetMeshPath("ball_ball.xml"));
  _TestXMLExportImportCPU(tledUnitTest::GetMeshPath("two_cubes_sheet.xml"));

#ifdef GPU_GP_CONTACT
  tledCUDAUnitTest::InitCUDATests();

  _TestXMLExportImportGPU(tledUnitTest::GetMeshPath("ball_ball.xml"));

  tledCUDAUnitTest::FinishCUDATests();
#endif

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
