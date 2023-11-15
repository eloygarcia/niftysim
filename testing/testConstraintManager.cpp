// =========================================================================
// File:       testConstraintManager.cpp
// Purpose:    tledConstraintManager unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledUnitTest.h"
#include "tledModel.h"
#include "tledSolverCPU.h"
#include "tledConstraintManager.h"
#include "tledMesh.h"

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cassert>

using namespace std;

static void _GenRandomConstraint(vector<vector<int> > &rvv_inds, vector<vector<float> > &rvv_mags, const int numNodes) {
  const int numInds = rand()%(numNodes - 3) + 3;

  vector<int> newC;

  newC.reserve(numInds);
  for (int n = 0; n < numInds; n++) {
    newC.push_back(rand()%numNodes);	
  }
  rvv_inds.push_back(tledHelper::MakeSortedUnique(newC));

  rvv_mags.push_back(vector<float>());
  rvv_mags.back().reserve(rvv_inds.back().size());
  for (int n = 0; n < (int)rvv_inds.back().size(); n++) {
    rvv_mags.back().push_back((float)(10*(drand48() - drand48())));
  }
}

static bool _WriteConstraint(tledUnitTestXMLWriter &r_writer, const vector<int> &inds, const vector<float> &mags, const string &type, const string &dof) {
  r_writer.GetFileStream() << "\t<Constraint ";
  if (dof == "all" && rand()%2 == 0) {
    cout << "Omitting DOF\n"; 
  } else {
    r_writer.GetFileStream() << "DOF=\"" << dof << "\" ";
  }
  r_writer.GetFileStream() << "Type=\"" << type << "\" SpecType=\"NODES\" LoadShape=\"STEP\" NumNodes=\"" << inds.size() << "\">\n"
			   << "\t\t<Nodes>\n";
  copy(inds.begin(), inds.end(), ostream_iterator<int>(r_writer.GetFileStream() << "\t\t\t", "\n\t\t\t"));
  r_writer.GetFileStream() << "</Nodes>\n\t\t<Magnitudes Type=\"DIFFORM\">\n";
  copy(mags.begin(), mags.end(), ostream_iterator<float>(r_writer.GetFileStream() << "\t\t\t", "\n\t\t\t"));    
  r_writer.GetFileStream() << "\t\t</Magnitudes>\n"
			   << "\t</Constraint>\n";

  return !r_writer.GetFileStream().fail();
}

static void _TestConstraintType(const vector<vector<int> > &refInds, const vector<vector<float> > &refMags, const vector<int> &testInds, const vector<float> &testMags) {
  vector<int> refIndsAll;
  vector<float> refMagsAll;

  for (int c = 0; c < (int)refInds.size(); c++) {
    refIndsAll.insert(refIndsAll.end(), refInds[c].begin(), refInds[c].end());
    refMagsAll.insert(refMagsAll.end(), refMags[c].begin(), refMags[c].end());    
  }

  tledUnitTestAssert(refIndsAll.size() == testInds.size());
  tledUnitTestAssert(refMagsAll.size() == testMags.size());
  
  for (int n = 0; n < (int)refIndsAll.size(); n++) {
    int tInd = 0;
    
    for (; tInd < (int)testInds.size() && fabs(testMags[tInd] - refMagsAll[n]) > 1e-5f; tInd++);
    tledUnitTestAssert(tInd < (int)testInds.size() && fabs(testMags[tInd] - refMagsAll[n]) <= 1e-5f);
  }
}

static void _GenRandomConstraints(vector<vector<int> > &rvv_inds, vector<vector<float> > &rvv_mags, const int numConsts, const int numNodes) {
  if (numConsts > 0) {
    rvv_inds.reserve(numConsts);
    rvv_mags.reserve(numConsts);

    for (int c = 0; c < numConsts; c++) {
      _GenRandomConstraint(rvv_inds, rvv_mags, numNodes);
    }
  }
}

static void _TestRandomConstraints(const string &meshPath, const string &meshType) {
  vector<vector<int> > forceInds, dispInds;
  vector<vector<float> > forceMags, dispMags;
  string xmlPath;
  int numNodes = tledUnitTest::LoadMSHMesh(meshPath, meshType.c_str()).GetNumNodes();

  _GenRandomConstraints(dispInds, dispMags, rand()%5, numNodes);
  _GenRandomConstraints(forceInds, forceMags, rand()%5, numNodes);

  {
    bool wasSuccess = true;
    tledUnitTestXMLWriter xmlWriter;

    xmlPath = xmlWriter.GetFilePath();
    wasSuccess &= xmlWriter.StartXML();
    xmlWriter.GetFileStream() << "\t<MSHMesh Type=\"" << meshType << "\">\n"
			      << "\t\t" << meshPath << endl
			      << "\t</MSHMesh>\n";
    for (int c = 0; c < (int)dispInds.size(); c++) {
      wasSuccess &= _WriteConstraint(xmlWriter, dispInds[c], dispMags[c], "Disp", "0");
    }
    for (int c = 0; c < (int)forceInds.size(); c++) {
      wasSuccess &= _WriteConstraint(xmlWriter, forceInds[c], forceMags[c], "Force", "0");
    }
    xmlWriter.GetFileStream() << "\t<ElementSet Size=\"all\">\n"
			      << "\t\t<Material Type=\"NH\">\n"      
			      << "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			      << "\t\t</Material>\n"
			      << "\t\t0\n"
			      << "\t</ElementSet>\n";
    wasSuccess &= !xmlWriter.GetFileStream().fail();
    wasSuccess &= xmlWriter.CloseXML();    
    assert(wasSuccess);
  }

  {
    tledModel model(xmlPath.c_str());
    tledSolverCPU solver;
    tledConstraintManager mgr(&model, &solver);
    std::vector<int> all;

    tledUnitTestAssert(mgr.GetNumConstraints() == int(dispInds.size() + forceInds.size()));
    mgr.GetDispInd(0);
    _TestConstraintType(dispInds, dispMags, *mgr.GetDispInd(0), *mgr.GetDispVal(0, 0, 0.0, 0.0));
    mgr.GetForceInd(0);
    _TestConstraintType(forceInds, forceMags, *mgr.GetForceInd(0), *mgr.GetForceVal(0, 0, 0.0, 0.0));
  }

  tledUnitTest::RemoveFile(xmlPath);
}

static void _TestRandomAllConstraints(const string &meshPath, const string &meshType) {
  vector<vector<int> > allForceInds, allDispInds;
  vector<vector<float> > allForceMags, allDispMags;
  string xmlPath;
  int numNodes = tledUnitTest::LoadMSHMesh(meshPath, meshType.c_str()).GetNumNodes();

  _GenRandomConstraints(allDispInds, allDispMags, rand()%5, numNodes);
  _GenRandomConstraints(allForceInds, allForceMags, rand()%5, numNodes);

  {
    bool wasSuccess = true;
    tledUnitTestXMLWriter xmlWriter;

    xmlPath = xmlWriter.GetFilePath();
    wasSuccess &= xmlWriter.StartXML();
    xmlWriter.GetFileStream() << "\t<MSHMesh Type=\"" << meshType << "\">\n"
			      << "\t\t" << meshPath << endl
			      << "\t</MSHMesh>\n";
    for (int c = 0; c < (int)allDispInds.size(); c++) {
      wasSuccess &= _WriteConstraint(xmlWriter, allDispInds[c], allDispMags[c], "Disp", "all");
    }
    for (int c = 0; c < (int)allForceInds.size(); c++) {
      wasSuccess &= _WriteConstraint(xmlWriter, allForceInds[c], allForceMags[c], "Force", "all");
    }        
    xmlWriter.GetFileStream() << "\t<ElementSet Size=\"all\">\n"
			      << "\t\t<Material Type=\"NH\">\n"      
			      << "\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			      << "\t\t</Material>\n"
			      << "\t\t0\n"
			      << "\t</ElementSet>\n";
    wasSuccess &= !xmlWriter.GetFileStream().fail();
    wasSuccess &= xmlWriter.CloseXML();    
    assert(wasSuccess);
  }

  {
    tledModel model(xmlPath.c_str());
    tledSolverCPU solver;
    tledConstraintManager mgr(&model, &solver);
    std::vector<int> all;

    for (int c = 0; c < 3; c++) {
      tledUnitTestAssert(mgr.GetNumConstraints() == 3*int(allDispInds.size() + allForceInds.size()));
      mgr.GetDispInd(c);
      _TestConstraintType(allDispInds, allDispMags, *mgr.GetDispInd(c), *mgr.GetDispVal(0, 0, 0.0, 0.0));
      mgr.GetForceInd(c);
      _TestConstraintType(allForceInds, allForceMags, *mgr.GetForceInd(c), *mgr.GetForceVal(0, 0, 0.0, 0.0));      
    }
  }

  tledUnitTest::RemoveFile(xmlPath);
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestRandomConstraints(tledUnitTest::GetMeshPath("box8.msh"), "H8");
  _TestRandomAllConstraints(tledUnitTest::GetMeshPath("box8.msh"), "H8");
  _TestRandomConstraints(tledUnitTest::GetMeshPath("collision_bar.msh"), "T4");
  _TestRandomAllConstraints(tledUnitTest::GetMeshPath("collision_bar.msh"), "T4");
  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
