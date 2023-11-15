// =========================================================================
// File:       testMeshTopology.cpp
// Purpose:    Mesh topology module unit test
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
#include "tledModel.h"
#include "tledHelper.h"
#include "tledSubModelManager.h"
#include "tledVectorArithmetic.h"

#include <cstdlib>
#include <iterator>
#include <algorithm>
#include <limits>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

template <const int t_numElNodes>
static void _TestNormalBoundaries(const string &meshPath, const vector<int> &nodeInds, const float normal[], const float angle, const std::vector<std::string> bdyRestrictions = std::vector<std::string>()) {
  bool wasSuccess = true;
  tledUnitTestXMLWriter xmlWriter;

  wasSuccess &= xmlWriter.StartXML();
  xmlWriter.GetFileStream() << "\t<MSHMesh Type=\"" << (t_numElNodes == 4? "T4" : "H8") << "\">\n"
			    << "\t\t" << meshPath << endl
			    << "\t</MSHMesh>\n";
  xmlWriter.GetFileStream() << "\t<Constraint DOF=\"0\" Type=\"Disp\" SpecType=\"NORMAL\" LoadShape=\"STEP\">\n"
			    << "\t\t<Normal ToleranceAngle=\"" << angle << "\">\n";   
  copy(normal, normal + 2, ostream_iterator<float>(xmlWriter.GetFileStream() << "\t\t\t", "\n\t\t\t"));
  xmlWriter.GetFileStream() << normal[2] << endl
			    << "\t\t</Normal>\n"
			    << "\t\t<Magnitudes Type=\"UNIFORM\">\n"
			    << "\t\t\t" << 1 << endl
			    << "\t\t</Magnitudes>\n";
  for (std::vector<std::string>::const_iterator ic_r = bdyRestrictions.begin(); ic_r < bdyRestrictions.end(); ic_r++) {
    xmlWriter.GetFileStream() << *ic_r;
  }
  xmlWriter.GetFileStream() << "\t</Constraint>\n";
  xmlWriter.GetFileStream() << "\t<ElementSet Size=\"all\">\n"
			    << "\t\t<Material Type=\"NH\">\n"      
			    << "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			    << "\t\t</Material>\n"
			    << "\t</ElementSet>\n";
  wasSuccess &= !xmlWriter.GetFileStream().fail();
  wasSuccess &= xmlWriter.CloseXML();
  assert(wasSuccess);

  {
    tledModel model(const_cast<char*>(xmlWriter.GetFilePath().c_str()));
    vector<int> inds;
    
    tledUnitTestAssert(model.GetNumConstraints() == 1);
    tledUnitTestAssert(model.GetConstraintDOF(0) == 0);
    inds = model.GetConstraintInd(0);      
    tledUnitTestAssert(tledHelper::MakeSortedUnique(inds).size() == inds.size());

    inds = tledHelper::MakeSortedUnique(inds);
    tledUnitTestAssert(equal(inds.begin(), inds.end(), nodeInds.begin()));
  }

  xmlWriter.CleanUp();
}

static void _TestNormalBoundariesOnBar() {
  const string meshPath = tledUnitTest::GetMeshPath("collision_bar.msh");
  const float normal[] = {0, 0, 1};
  const float angle = 1;

  tledMesh mesh = tledUnitTest::LoadMSHMesh(meshPath, "T4");

  {
    vector<int> nodeInds;

    for (float const *pc_nodeCd = mesh.GetAllNodeCds() + 2; pc_nodeCd < mesh.GetAllNodeCds() + 3*mesh.GetNumNodes(); pc_nodeCd += 3) {
      if (std::fabs(*pc_nodeCd - 0.25) < 0.01) nodeInds.push_back((pc_nodeCd - mesh.GetAllNodeCds())/3);
    }
    nodeInds = tledHelper::MakeSortedUnique(nodeInds);

    _TestNormalBoundaries<4>(meshPath, nodeInds, normal, angle);
  }

  {
    float barBounds[3][2], restBounds[3][2];
    vector<int> nodeInds;
    std::vector<std::string> restrictions;
    std::ostringstream oss;

    for (int c = 0; c < 3; c++) barBounds[c][1] = -(barBounds[c][0] = std::numeric_limits<float>::max());
    for (float const *pc_n = mesh.GetAllNodeCds(); pc_n < mesh.GetAllNodeCds() + 3*mesh.GetNumNodes();) {
      for (int c = 0; c < 3; c++, pc_n++) {
	barBounds[c][0] = std::min(barBounds[c][0], *pc_n);
	barBounds[c][1] = std::max(barBounds[c][1], *pc_n);
      }
    }
    assert((barBounds[0][1] - barBounds[0][0])*(barBounds[1][1] - barBounds[1][0])*(barBounds[2][1] - barBounds[2][0]) > 1e-3);

    restBounds[0][0] = barBounds[0][0] - 1e-2f;
    restBounds[0][1] = 0.75f*barBounds[0][0] + 0.25f*barBounds[0][1];

    restBounds[1][1] = barBounds[1][1] + 1e-2f;
    restBounds[1][0] = (barBounds[1][1] + barBounds[1][0])/2 - 1e-2f;

    restBounds[2][0] = barBounds[2][0] - 1e-2f;
    restBounds[2][1] = barBounds[2][1] + 1e-2f;

    oss << "\t\t\t<RestrictTo Type=\"Box\">\n";
    for (int c = 0; c < 3; c++) for (int b = 0; b < 2; b++) {
	oss << restBounds[c][b] << " ";
      }
    oss << "\n\t\t\t</RestrictTo>\n";
    restrictions.push_back(oss.str());

    for (float const *pc_nodeCd = mesh.GetAllNodeCds(); pc_nodeCd < mesh.GetAllNodeCds() + 3*mesh.GetNumNodes(); pc_nodeCd += 3) {
      if (pc_nodeCd[0] >= restBounds[0][0] && pc_nodeCd[0] <= restBounds[0][1] && pc_nodeCd[1] >= restBounds[1][0] && pc_nodeCd[1] <= restBounds[1][1]) {
	if (std::fabs(*(pc_nodeCd + 2) - 0.25) < 0.01) nodeInds.push_back((pc_nodeCd - mesh.GetAllNodeCds())/3);
      }
    }
    nodeInds = tledHelper::MakeSortedUnique(nodeInds);

    _TestNormalBoundaries<4>(meshPath, nodeInds, normal, angle, restrictions);        
  }

  {
    float barBounds[3][2], cnt[3], r;
    vector<int> nodeInds;
    std::vector<std::string> restrictions;
    std::ostringstream oss;

    for (int c = 0; c < 3; c++) barBounds[c][1] = -(barBounds[c][0] = std::numeric_limits<float>::max());
    for (float const *pc_n = mesh.GetAllNodeCds(); pc_n < mesh.GetAllNodeCds() + 3*mesh.GetNumNodes();) {
      for (int c = 0; c < 3; c++, pc_n++) {
	barBounds[c][0] = std::min(barBounds[c][0], *pc_n);
	barBounds[c][1] = std::max(barBounds[c][1], *pc_n);
      }
    }
    cnt[2] = barBounds[2][1];
    cnt[1] = (barBounds[1][0] + barBounds[1][1])/2;
    cnt[0] = (barBounds[0][0] + barBounds[0][1])/2;
    r = std::min((barBounds[1][1] - barBounds[1][0])/2, (barBounds[0][1] - barBounds[0][0])/2);

    oss << "\t\t\t<RestrictTo Type=\"Sphere\">\n";
    oss << r << " ";
    for (int c = 0; c < 3; c++) {
	oss << cnt[c] << " ";
    }
    oss << "\n\t\t\t</RestrictTo>\n";
    restrictions.push_back(oss.str());

    for (float const *pc_nodeCd = mesh.GetAllNodeCds(); pc_nodeCd < mesh.GetAllNodeCds() + 3*mesh.GetNumNodes(); pc_nodeCd += 3) {
      using namespace tledVectorArithmetic;

      float d[3];
      
      if (Norm(Sub(d, pc_nodeCd, cnt))) {
	nodeInds.push_back((pc_nodeCd - mesh.GetAllNodeCds())/3);
      }
    }
    nodeInds = tledHelper::MakeSortedUnique(nodeInds);

    _TestNormalBoundaries<4>(meshPath, nodeInds, normal, angle, restrictions);        
  }
} 

static void _TestNormalBoundariesOnBox() {
  const string meshPath = tledUnitTest::GetMeshPath("box8.msh");
  const float normal[] = {0, -1, 0};
  const float angle = 1;

  tledMesh mesh = tledUnitTest::LoadMSHMesh(meshPath, "H8");
  vector<int> nodeInds;
  float const *pc_nodeCd;

  for (pc_nodeCd = mesh.GetAllNodeCds() + 1; pc_nodeCd < mesh.GetAllNodeCds() + 3*mesh.GetNumNodes(); pc_nodeCd += 3) {
    if (std::fabs(*pc_nodeCd + 0.5) < 0.01) nodeInds.push_back((pc_nodeCd - mesh.GetAllNodeCds())/3);
  }
  nodeInds = tledHelper::MakeSortedUnique(nodeInds);

  _TestNormalBoundaries<8>(meshPath, nodeInds, normal, angle);
} 

static void _TestConstraintAllNodes(const std::string &meshPath, const std::string &meshtype) {
  bool wasSuccess = true;
  tledUnitTestXMLWriter xmlWriter;

  wasSuccess &= xmlWriter.StartXML();
  xmlWriter.GetFileStream() << "\t<MSHMesh Type=\"" << meshtype << "\">\n"
			    << "\t\t" << meshPath << endl
			    << "\t</MSHMesh>\n";
  xmlWriter.GetFileStream() << "\t<Constraint DOF=\"0\" NumNodes=\"all\" Type=\"Disp\" SpecType=\"NODES\" LoadShape=\"STEP\">\n"
			    << "\t\t<Magnitudes Type=\"UNIFORM\">\n"
			    << "\t\t\t1\n"
			    << "\t\t</Magnitudes>\n"
			    << "\t</Constraint>\n";
  xmlWriter.GetFileStream() << "\t<ElementSet Size=\"all\">\n"
			    << "\t\t<Material Type=\"NH\">\n"      
			    << "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			    << "\t\t</Material>\n"
			    << "\t</ElementSet>\n";
  wasSuccess &= !xmlWriter.GetFileStream().fail();
  wasSuccess &= xmlWriter.CloseXML();
  assert(wasSuccess);

  {
    tledModel model(const_cast<char*>(xmlWriter.GetFilePath().c_str()));
    tledMesh &mesh = *model.GetMesh();
    vector<int> inds;
    
    tledUnitTestAssert(model.GetNumConstraints() == 1);
    tledUnitTestAssert(model.GetConstraintDOF(0) == 0);
    inds = model.GetConstraintInd(0);      
    tledUnitTestAssert((int)tledHelper::MakeSortedUnique(inds).size() == mesh.GetNumNodes());

    inds = tledHelper::MakeSortedUnique(inds);
    tledUnitTestAssert(equal(inds.begin(), inds.end(), tledSequenceGenerator::MakeSequence(0, mesh.GetNumNodes()).begin()));
  }  

  xmlWriter.CleanUp();
}

template <const int t_numFacetVtcs>
std::vector<int> _ExtractCubeFaces(const tledMesh &cube, const float normal[]) {
  tledMeshSurface<t_numFacetVtcs> surface(cube);
  std::vector<int> faces;

  /* Not ideal: using mesh topology, surface classes in very basic unit test. Disable test if problems with one of these classes. */
  for (int f = 0; f < surface.GetNumberOfFacets(); f++) {
    float tmpN[3];

    if (std::fabs(tledVectorArithmetic::Dot(surface.ComputeNormalisedFacetNormal(tmpN, f), normal) - 1.f) < 1e-3f) {
      faces.insert(faces.end(), surface.GetFacet(f).NodeIndices, surface.GetFacet(f).NodeIndices + t_numFacetVtcs);
    }
  }

  return faces;
}

static void _TestSurfaceConstraintOnBox(const std::string meshType, const std::string constType) {
  const std::string meshPath = meshType == "T4"? tledUnitTest::GetMeshPath("box.msh") : tledUnitTest::GetMeshPath("box8.msh");
  const int numFacetVtcs = meshType == "T4"? 3 : 4;
  static const float normal[] = {0.f, 0.f, -1.f};

  tledMesh cube = tledUnitTest::LoadMSHMesh(meshPath, meshType.c_str());
  std::vector<int> faces = meshType == "T4"? _ExtractCubeFaces<3>(cube, normal) : _ExtractCubeFaces<4>(cube, normal);
  tledUnitTestXMLWriter xmlWriter;

  assert(faces.size() >= 4);

  {
    bool wasSuccess = true;

    wasSuccess &= xmlWriter.StartXML();
    xmlWriter.GetFileStream() << "\t<MSHMesh Type=\"" << meshType << "\">\n"
			      << "\t\t" << meshPath << endl
			      << "\t</MSHMesh>\n";    
    xmlWriter.GetFileStream() << "\t<Constraint Type=\"" << constType << "\" NumFaces=\"" << faces.size()/numFacetVtcs << "\" LoadShape=\"STEP\">\n";
    copy(faces.begin(), faces.end(), ostream_iterator<int>(xmlWriter.GetFileStream(), "\n"));
    if (constType == "Pressure") {
      xmlWriter.GetFileStream() << "\t\t<Magnitude>\n"
				<< "\t\t\t1\n"
				<< "\t\t</Magnitude>\n"
				<< "\t</Constraint>\n";
    } else {
      xmlWriter.GetFileStream() << "\t\t<Magnitudes Type=\"UNIFORM\">\n"
				<< "\t\t\t1 0 0\n"
				<< "\t\t</Magnitudes>\n"
				<< "\t</Constraint>\n";
    }
    xmlWriter.GetFileStream() << "\t<ElementSet Size=\"all\">\n"
			      << "\t\t<Material Type=\"NH\">\n"      
			      << "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			      << "\t\t</Material>\n"
			      << "\t</ElementSet>\n";
    wasSuccess &= !xmlWriter.GetFileStream().fail();
    wasSuccess &= xmlWriter.CloseXML();
    assert(wasSuccess);
  }
  
  {
    tledModel model(const_cast<char*>(xmlWriter.GetFilePath().c_str()));
    vector<int> inds, refInds;
    
    tledUnitTestAssert(model.GetNumConstraints() == 1);
    if (constType == "Pressure") {
      inds = model.GetPressureFaceNodeInds(0);          
      tledUnitTestAssert(int(faces.size()/numFacetVtcs) == model.GetPressureNumFaces(0));
      tledUnitTestAssert(std::fabs(model.GetPressureMagnitude(0) - 1.f) < 1e-3f);
    } else {
      vector<float> mags;

      inds = model.GetTractionFaceNodeInds(0);          
      tledUnitTestAssert(int(faces.size()/numFacetVtcs) == model.GetTractionNumFaces(0));
      mags = model.GetTractionFaceTractions(0);
      tledUnitTestAssert(mags.size() == 3*faces.size()/numFacetVtcs || mags.size() == 3);
      for (vector<float>::const_iterator ic_t = mags.begin(); ic_t < mags.end(); ic_t += 3) {
	tledUnitTestAssert(std::fabs(*ic_t - 1.f) + std::fabs(*(ic_t + 1)) + std::fabs(*(ic_t + 2)) < 1e-3f);
      }
    }

    inds = tledHelper::MakeSortedUnique(inds);
    refInds = tledHelper::MakeSortedUnique(faces);
    tledUnitTestAssert(equal(inds.begin(), inds.end(), refInds.begin()));
  }  
  
  xmlWriter.CleanUp();
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestNormalBoundariesOnBar();
  _TestNormalBoundariesOnBox();

  _TestConstraintAllNodes(tledUnitTest::GetMeshPath("box8.msh"), "H8");
  _TestConstraintAllNodes(tledUnitTest::GetMeshPath("collision_bar.msh"), "T4");

  _TestSurfaceConstraintOnBox("H8", "Pressure");
  _TestSurfaceConstraintOnBox("T4", "Pressure");
  _TestSurfaceConstraintOnBox("H8", "Traction");
  _TestSurfaceConstraintOnBox("T4", "Traction");

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
