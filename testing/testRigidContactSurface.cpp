// =========================================================================
// File:       testRigidContactSurface.cpp
// Purpose:    tledRigidContactSurface unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    August 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledUnitTest.h"
#include "tledSimulator.h"
#include "tledModel.h"
#include "tledRigidContactSurfaceCPU.h"
#include "tledMovingRigidContactSurfaceCPU.h"

#include <string>

using namespace tledVectorArithmetic;

template <class TSurface>
static TSurface& _LoadSurface(tledModel* &rp_model, tledSimulator* &rp_sim, const std::string &xmlPath, const int rigidSurfaceInd = 0) {
  rp_model = new tledModel(xmlPath.c_str());
  rp_sim = new tledSimulator(rp_model);

  assert(rigidSurfaceInd < rp_sim->GetContactManager()->GetUnstructuredContactManager().GetNumberOfRigidSurfaces());

  return rp_sim->GetContactManager()->GetUnstructuredContactManager().GetRigidSurface<TSurface>(rigidSurfaceInd);
}

tledUnitTestXMLWriter* _GenerateVTKSimXML(const std::string &solidMeshPath, const char solidType[], const std::string &rigidVTKMeshPath, const char surfaceType[], const bool doMoving = false) {
  bool wasSuccess = true;
  tledUnitTestXMLWriter *p_xmlWriter = new tledUnitTestXMLWriter;

  wasSuccess &= p_xmlWriter->StartXML();
  p_xmlWriter->GetFileStream() << "\t<MSHMesh Type=\"" << solidType << "\">\n"
			       << "\t\t" << solidMeshPath << std::endl
			       << "\t</MSHMesh>\n";
  p_xmlWriter->GetFileStream() << "\t<ElementSet Size=\"all\">\n"
			       << "\t\t<Material Type=\"NH\">\n"      
			       << "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			       << "\t\t</Material>\n"
			       << "\t\t0\n"
			       << "\t</ElementSet>\n";
  p_xmlWriter->GetFileStream() << "\t<SystemParams>\n"
			       << "\t\t<TimeStep>1e-3</TimeStep>\n"
			       << "\t\t<TotalTime>0.1</TotalTime>\n"
			       << "\t\t<DampingCoeff>0.1</DampingCoeff>\n"
			       << "\t\t<Density>1000</Density>\n"
			       << "\t</SystemParams>\n";
  p_xmlWriter->GetFileStream() << "\t<ContactSurface>\n"
			       << "\t\t<VTKSurface Type=\"" << surfaceType << "\">\n"
			       << "\t\t\t" << rigidVTKMeshPath << std::endl
			       << "\t\t</VTKSurface>\n";
  if (doMoving) {
    p_xmlWriter->GetFileStream() << "\t\t<Motion Type=\"Translation\">\n"
				 << "\t\t\t1.0 0.0  1.0\n"
				 << "\t\t</Motion>\n";
  }
  p_xmlWriter->GetFileStream() << "\t</ContactSurface>\n";
  wasSuccess &= !p_xmlWriter->GetFileStream().fail();
  wasSuccess &= p_xmlWriter->CloseXML();
  assert(wasSuccess);  

  return p_xmlWriter;
}

tledUnitTestXMLWriter* _GenerateVTKTransformXML(const std::string &solidMeshPath, const char solidType[], const std::string &rigidVTKMeshPath, const char surfaceType[], const float t[], const float s, const float cor[], const float angles[], const bool doMoving) {
  bool wasSuccess = true;
  tledUnitTestXMLWriter *p_xmlWriter = new tledUnitTestXMLWriter;

  wasSuccess &= p_xmlWriter->StartXML();
  p_xmlWriter->GetFileStream() << "\t<MSHMesh Type=\"" << solidType << "\">\n"
			       << "\t\t" << solidMeshPath << std::endl
			       << "\t</MSHMesh>\n";
  p_xmlWriter->GetFileStream() << "\t<ElementSet Size=\"all\">\n"
			       << "\t\t<Material Type=\"NH\">\n"      
			       << "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			       << "\t\t</Material>\n"
			       << "\t\t0\n"
			       << "\t</ElementSet>\n";
  p_xmlWriter->GetFileStream() << "\t<SystemParams>\n"
			       << "\t\t<TimeStep>1e-3</TimeStep>\n"
			       << "\t\t<TotalTime>0.1</TotalTime>\n"
			       << "\t\t<DampingCoeff>0.1</DampingCoeff>\n"
			       << "\t\t<Density>1000</Density>\n"
			       << "\t</SystemParams>\n";
  p_xmlWriter->GetFileStream() << "\t<ContactSurface>\n"
			       << "\t\t<VTKSurface Type=\"" << surfaceType << "\">\n"
			       << "\t\t\t" << rigidVTKMeshPath << std::endl;
  if (Norm(t) != 0.f) {
    p_xmlWriter->GetFileStream() << "\t\t\t<Translation>\n"
				 << "\t\t\t" << t[0] << " " << t[1] << " " << t[2] << std::endl
				 << "\t\t\t</Translation>\n";
  }

  if (s != 1.f) {
    p_xmlWriter->GetFileStream() << "\t\t\t<ScaleFactor>\n"
				 << "\t\t\t" << s << std::endl
				 << "\t\t\t</ScaleFactor>\n";
  }

  if (Norm(angles) != 0.f) {
    p_xmlWriter->GetFileStream() << "\t\t\t<Rotation>\n";
    std::copy(cor, cor + 3, std::ostream_iterator<float>(p_xmlWriter->GetFileStream() << "\t\t\t", " "));
    std::copy(angles, angles + 3, std::ostream_iterator<float>(p_xmlWriter->GetFileStream(), " "));
    p_xmlWriter->GetFileStream() << "\n\t\t\t</Rotation>\n";
  }

  p_xmlWriter->GetFileStream() << "\t\t</VTKSurface>\n";
  if (doMoving) {
    p_xmlWriter->GetFileStream() << "\t\t<Motion Type=\"Translation\">\n"
				 << "\t\t\t1.0 0.0  1.0\n"
				 << "\t\t</Motion>\n";
  }
  p_xmlWriter->GetFileStream() << "\t</ContactSurface>\n";
  wasSuccess &= !p_xmlWriter->GetFileStream().fail();
  wasSuccess &= p_xmlWriter->CloseXML();
  assert(wasSuccess);  

  return p_xmlWriter;
}

template <class TSurface>
static void _TestLoadTransform(const bool doMotion) {  
  const std::string solidPath = tledUnitTest::GetMeshPath("collision_bar.msh"), surfacePath = tledUnitTest::GetMeshPath("box_surface.vtk");

  tledModel *p_refModel;
  tledSimulator *p_refSim;
  TSurface const *pc_refSurface;

  {
    tledUnitTestXMLWriter *p_refWriter = _GenerateVTKSimXML(solidPath, "T4", surfacePath, "T3", doMotion);
    
    pc_refSurface = &_LoadSurface<TSurface>(p_refModel, p_refSim, p_refWriter->GetFilePath());
    p_refWriter->CleanUp();
    delete p_refWriter;
  }

  for (int angleInd = 0; angleInd < 3; angleInd++) {
    const float refAngles[3] = {45.f, 73.98f, 67.5f};
    const float angle = refAngles[angleInd];

    for (int axis = 0; axis < 3; axis++) {
      using namespace tledVectorArithmetic;

      const float cor[] = {0.25f, 0.25f, 0.25f};          
      const float t[] = {1.f, 2.2f, 3.4f};

      tledModel *p_testModel;
      tledSimulator *p_testSim;
      TSurface const *pc_testSurface;
      float angles[3];

      std::fill(angles, angles + 3, 0.f);
      angles[axis] = angle;    

      {
	tledUnitTestXMLWriter *p_writer = _GenerateVTKTransformXML(solidPath, "T4", surfacePath, "T3", t, 1.f, cor, angles, doMotion);

	pc_testSurface = &_LoadSurface<TSurface>(p_testModel, p_testSim, p_writer->GetFilePath());
	p_writer->CleanUp();
	delete p_writer;    
      }

      tledUnitTestAssert(pc_refSurface->GetNumberOfNodes() == pc_testSurface->GetNumberOfNodes());
      for (int n = 0; n < pc_refSurface->GetNumberOfNodes(); n++) {
	float cntRef[3], cntTest[3];

	Sub(cntRef, pc_refSurface->GetNodeCoordinates(n), cor);
	Sub(cntTest, Sub(cntTest, pc_testSurface->GetNodeCoordinates(n), t), cor);
	cntTest[axis] = cntRef[axis] = 0.f;

	if (!tledHelper::IsNumericallyZero(Norm(cntTest), std::fabs(angle)) && !tledHelper::IsNumericallyZero(Norm(cntRef), std::fabs(angle))) {
	  tledUnitTestAssert(std::fabs(ComputeAngle(cntRef, cntTest) - tledPi*angle/180.f)/tledPi*180.f <= tledPi*std::fabs(angle)*1e-2f);
	}
      }

      delete p_testSim;
      delete p_testModel;
    }
  }

  delete p_refSim;
  delete p_refModel;
}

template <class TSurface>
static void _TestLoadNormals(const std::string &xmlPath) {
  tledModel *p_model;
  tledSimulator *p_sim;
  TSurface &r_surface = _LoadSurface<TSurface>(p_model, p_sim, xmlPath);

  for (int f = 0; f < r_surface.GetNumberOfFacets(); f++) {
    tledUnitTestAssert(std::fabs(Norm(r_surface.GetFacetNormal(f)) - 1.f) < 1e-4f);
  }

  for (int n = 0; n < r_surface.GetNumberOfNodes(); n++) {
    tledUnitTestAssert(std::fabs(Norm(r_surface.GetNodeNormal(n)) - 1.f) < 1e-4f);
  }

  delete p_sim;
  delete p_model;
}

template <class TSurface>
static void _TestLoadEdges(const std::string &xmlPath) {
  tledModel *p_model;
  tledSimulator *p_sim;
  TSurface &r_surface = _LoadSurface<TSurface>(p_model, p_sim, xmlPath);
  
  for (int f = 0; f < r_surface.GetNumberOfFacets(); f++) {
    const int *facet = r_surface.GetFacet(f).NodeIndices;

    for (int v = 0; v < TSurface::Facet::NumberOfVertices; v++) {
      const int e0 = facet[v], e1 = facet[(v+1)%TSurface::Facet::NumberOfVertices];

      std::vector<std::pair<int, int> >::const_iterator ic_e;

      for (ic_e = r_surface.GetAllEdges().begin(); ic_e < r_surface.GetAllEdges().end(); ic_e++) {
	if ((ic_e->first == e0 && ic_e->second == e1) || (ic_e->first == e1 && ic_e->second == e0)) {
	  break;
	}
      }
      tledUnitTestAssert(ic_e < r_surface.GetAllEdges().end());
    }
  }

  delete p_sim;
  delete p_model;
}

int main(int argc, char *argv[]) {
  tledUnitTestXMLWriter *p_writer;

  tledUnitTest::InitUnitTest();
  _TestLoadTransform<tledMovingRigidContactSurfaceT3CPU>(true);
  _TestLoadTransform<tledRigidContactSurfaceT3CPU>(false);

  p_writer = _GenerateVTKSimXML(tledUnitTest::GetMeshPath("collision_bar.msh"), "T4", tledUnitTest::GetMeshPath("box_surface.vtk"), "T3");
  _TestLoadNormals<tledRigidContactSurfaceT3CPU>(p_writer->GetFilePath());
  _TestLoadEdges<tledRigidContactSurfaceT3CPU>(p_writer->GetFilePath());
  p_writer->CleanUp();
  delete p_writer;

  p_writer = _GenerateVTKSimXML(tledUnitTest::GetMeshPath("collision_bar.msh"), "T4", tledUnitTest::GetMeshPath("box_surface.vtk"), "T3", true);
  _TestLoadNormals<tledMovingRigidContactSurfaceT3CPU>(p_writer->GetFilePath());
  _TestLoadEdges<tledMovingRigidContactSurfaceT3CPU>(p_writer->GetFilePath());
  p_writer->CleanUp();
  delete p_writer;

  p_writer = _GenerateVTKSimXML(tledUnitTest::GetMeshPath("collision_bar.msh"), "T4", tledUnitTest::GetMeshPath("sphere_surface.vtk"), "T3");
  _TestLoadNormals<tledRigidContactSurfaceT3CPU>(p_writer->GetFilePath());
  _TestLoadEdges<tledRigidContactSurfaceT3CPU>(p_writer->GetFilePath());
  p_writer->CleanUp();
  delete p_writer;

  p_writer = _GenerateVTKSimXML(tledUnitTest::GetMeshPath("collision_bar.msh"), "T4", tledUnitTest::GetMeshPath("sphere_surface.vtk"), "T3", true);
  _TestLoadNormals<tledMovingRigidContactSurfaceT3CPU>(p_writer->GetFilePath());
  _TestLoadEdges<tledMovingRigidContactSurfaceT3CPU>(p_writer->GetFilePath());
  p_writer->CleanUp();
  delete p_writer;

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
