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
#include <vector>
#include <fstream>

using namespace std;

/*
 * Extracts a submesh.
 * Returns the submesh nodes, elements, and a node index map.
 */
static void _SplitMesh(float sMNodes[][3], int sMEls[][4], int &r_numSMNodes, vector<int> &r_glob2SMNodeMap, const float nodes[][3], const int els[][4], const int sMStartElInd, const int sMEndElInd, const int numMeshNodes, const float nodeTrans[]) {
  vector<int> sM2GlobNodeMap;
  int elInd;
  
  assert(sMStartElInd < sMEndElInd);
  for (elInd = sMStartElInd; elInd < sMEndElInd; elInd++) {
    sM2GlobNodeMap.insert(sM2GlobNodeMap.end(), els[elInd], els[elInd] + 4);    
  }
  sM2GlobNodeMap = tledHelper::MakeSortedUnique(sM2GlobNodeMap);
  r_numSMNodes = sM2GlobNodeMap.size();

  {
    vector<int>::const_iterator ic_globNodeInd;

    r_glob2SMNodeMap.insert(r_glob2SMNodeMap.end(), numMeshNodes, -1);
    for (ic_globNodeInd = sM2GlobNodeMap.begin(); ic_globNodeInd < sM2GlobNodeMap.end(); ic_globNodeInd++) {
      const int sMNodeInd = ic_globNodeInd - sM2GlobNodeMap.begin();

      r_glob2SMNodeMap[*ic_globNodeInd] = sMNodeInd;
      tledVectorArithmetic::Add(sMNodes[sMNodeInd], nodes[*ic_globNodeInd], nodeTrans);
    }
  }

  for (elInd = sMStartElInd; elInd < sMEndElInd; elInd++) {
    int vInd;

    for (vInd = 0; vInd < 4; vInd++) {
      assert(r_glob2SMNodeMap[els[elInd][vInd]] >= 0);
      sMEls[elInd-sMStartElInd][vInd] = r_glob2SMNodeMap[els[elInd][vInd]];
      assert(sMEls[elInd-sMStartElInd][vInd] >= 0 && sMEls[elInd-sMStartElInd][vInd] < r_numSMNodes);
    }
  }
}

static void _TestConstraints(const string &meshPath) {
  static const int maxNumEls = (int)2e3;
  static const int maxNumNodes = (int)1e3;
  static const int numTests = 100;
  static const int smDispMags[] = {10, -20};

  bool wasSuccess;
  int els[maxNumEls][4], numEls, numNodes, numSurfEls, testInd;
  float nodes[maxNumNodes][3];

  assert(tledUnitTest::ParseMSH(NULL, numNodes, NULL, numEls, NULL, numSurfEls, meshPath));
  assert(numNodes <= maxNumNodes && numEls <= maxNumEls);

  wasSuccess = tledUnitTest::ParseMSH(nodes, numNodes, els, numEls, NULL, numSurfEls, meshPath);
  assert(wasSuccess);

  for (testInd = 0; testInd < numTests && wasSuccess; testInd++) {
    static const float null[] = {0, 0, 0};
    const int bdyElInd = rand()%(numEls/3) + numEls/3;

    int sm0Els[maxNumEls][4], sm1Els[maxNumEls][4], numSMNodes[2];
    float sm0Nodes[maxNumNodes][3], sm1Nodes[maxNumNodes][3];
    vector<int> constNodeInds;
    tledUnitTestXMLWriter xmlWriter;
    vector<int> smNodeIndMaps[2];

    wasSuccess &= xmlWriter.StartXML();

    _SplitMesh(sm0Nodes, sm0Els, numSMNodes[0], smNodeIndMaps[0], nodes, els, 0, bdyElInd, numNodes, null);
    _SplitMesh(sm1Nodes, sm1Els, numSMNodes[1], smNodeIndMaps[1], nodes, els, bdyElInd, numEls, numNodes, null);
    xmlWriter.GetFileStream() << "\t<SubModel>\n";
    wasSuccess &= xmlWriter.WriteMesh(sm0Nodes, numSMNodes[0], sm0Els, bdyElInd);    

    constNodeInds = tledSequenceGenerator::MakeSequence(0, numSMNodes[0]);
    xmlWriter.GetFileStream() << "\t\t<Constraint DOF=\"0\" Type=\"Disp\" NumNodes=\"" << numSMNodes[0] << "\" LoadShape=\"STEP\">\n"
			      << "\t\t\t<Nodes>\n";   
    copy(constNodeInds.begin(), constNodeInds.end() - 1, ostream_iterator<int>(xmlWriter.GetFileStream() << "\t\t\t\t", "\n\t\t\t\t"));
    xmlWriter.GetFileStream() << constNodeInds.back() << endl
			      << "\t\t\t</Nodes>\n"
			      << "\t\t\t<Magnitudes Type=\"UNIFORM\">\n"
			      << "\t\t\t\t" << smDispMags[0] << endl
			      << "\t\t\t</Magnitudes>\n"
			      << "\t\t</Constraint>\n";
    xmlWriter.GetFileStream() << "\t\t<ElementSet Size=\"" << bdyElInd << "\">\n"
			      << "\t\t\t<Material Type=\"NH\">\n"      
			      << "\t\t\t\t<ElasticParams NumParams=\"2\">100 1000</ElasticParams>\n"
			      << "\t\t\t</Material>\n"
			      << "\t\t\t0\n"
			      << "\t\t</ElementSet>\n";
    xmlWriter.GetFileStream() << "\t</SubModel>\n";

    xmlWriter.GetFileStream() << "\t<SubModel>\n";
    wasSuccess &= xmlWriter.WriteMesh(sm1Nodes, numSMNodes[1], sm1Els, numEls - bdyElInd);    

    constNodeInds = tledSequenceGenerator::MakeSequence(0, numSMNodes[1]);
    xmlWriter.GetFileStream() << "\t\t<Constraint DOF=\"1\" Type=\"Disp\" NumNodes=\"" << numSMNodes[1] << "\" LoadShape=\"RAMP\">\n"
			      << "\t\t\t<Nodes>\n";   
    copy(constNodeInds.begin(), constNodeInds.end() - 1, ostream_iterator<int>(xmlWriter.GetFileStream() << "\t\t\t\t", "\n\t\t\t\t"));
    xmlWriter.GetFileStream() << constNodeInds.back() << endl
			      << "\t\t\t</Nodes>\n"
			      << "\t\t\t<Magnitudes Type=\"UNIFORM\">\n"
			      << "\t\t\t\t" << smDispMags[1] << endl
			      << "\t\t\t</Magnitudes>\n"
			      << "\t\t</Constraint>\n";
    xmlWriter.GetFileStream() << "\t\t<ElementSet Size=\"" << numEls - bdyElInd << "\">\n"
			      << "\t\t\t<Material Type=\"NH\">\n"
			      << "\t\t\t\t<ElasticParams NumParams=\"2\">1000 10000</ElasticParams>\n"
			      << "\t\t\t</Material>\n"
			      << "\t\t\t0\n"
			      << "\t\t</ElementSet>\n";
    xmlWriter.GetFileStream() << "\t</SubModel>\n";

    wasSuccess &= !xmlWriter.GetFileStream().fail();
    wasSuccess &= xmlWriter.CloseXML();
    assert(wasSuccess);

    {
      tledModel model(const_cast<char*>(xmlWriter.GetFilePath().c_str()));
      int cInd;
    
      tledUnitTestAssert(model.GetNumConstraints() == 2);
      for (cInd = 0; cInd < 2; cInd++) {
	vector<float> mags;
	vector<int> inds;
	vector<int>::iterator i_ind;

	tledUnitTestAssert(model.GetConstraintDOF(cInd) == cInd);
	tledUnitTestAssert(model.GetConstraintLoadShape(cInd) == (cInd == 0? STEP : RAMP));      
	mags = tledHelper::MakeSortedUnique(model.GetConstraintMag(cInd));      
	tledUnitTestAssert(mags.size() == 1 && mags.front() == smDispMags[cInd]); 
	inds = model.GetConstraintInd(cInd);      
	tledUnitTestAssert((int)tledHelper::MakeSortedUnique(inds).size() == numSMNodes[cInd]);
	tledUnitTestAssert(equal(inds.begin(), inds.end(), model.GetSubModelManager().GetSubMeshToMeshNodeIndexMap(cInd).begin()));
      } /* for sub domains */
    }

    xmlWriter.CleanUp();
  } /* for tests */
  assert(testInd == numTests);
}

static void _TestSplitMeshModels(const string &meshPath) {
  static const int maxNumEls = (int)2e3;
  static const int maxNumNodes = (int)1e3;
  static const int maxNumSubModels = 6;

  bool wasSuccess;
  int els[maxNumEls][4], numEls, numNodes, numSurfEls, numSMs, sMInd, lastElInd;
  float nodes[maxNumNodes][3];
  tledUnitTestXMLWriter refXMLWriter, sMXMLWriter;
  float transMag;
  

  assert(tledUnitTest::ParseMSH(NULL, numNodes, NULL, numEls, NULL, numSurfEls, meshPath));
  assert(numNodes <= maxNumNodes && numEls <= maxNumEls);

  wasSuccess = tledUnitTest::ParseMSH(nodes, numNodes, els, numEls, NULL, numSurfEls, meshPath);
  wasSuccess &= refXMLWriter.StartXML();
  wasSuccess &= refXMLWriter.WriteMesh(nodes, numNodes, els, numEls);
  wasSuccess &= refXMLWriter.CloseXML();

  wasSuccess &= sMXMLWriter.StartXML();
  transMag = drand48()*1e-2f;
  sMXMLWriter.GetFileStream() << "\t<SubModelNodeMergerDistance> " << 2*transMag << " </SubModelNodeMergerDistance>\n";
  numSMs = 2 + rand()%(maxNumSubModels - 2);
  lastElInd = 0;
  for (sMInd = 0; sMInd < numSMs && wasSuccess; sMInd++) {
    int sMElIndBounds[2], sMEls[maxNumEls][4], numSMNodes;
    float sMNodes[maxNumNodes][3], trans[3];
    vector<int> nodeIndMap;

    sMElIndBounds[0] = lastElInd;
    if (sMInd == numSMs - 1) lastElInd = sMElIndBounds[1] = numEls;
    else lastElInd = sMElIndBounds[1] = rand()%(numEls - 1 - sMElIndBounds[0] - (numSMs - 1 - sMInd)) + 1 + sMElIndBounds[0];

    trans[0] = drand48() - drand48() + (rand()%2 == 0? -1 : 1);
    trans[1] = drand48() - drand48() + (rand()%2 == 0? -1 : 1);
    trans[2] = drand48() - drand48() + (rand()%2 == 0? -1 : 1);
    tledVectorArithmetic::ScalarMul(trans, transMag/tledVectorArithmetic::Norm(trans)*drand48());    
    assert(tledVectorArithmetic::Norm(trans) < transMag);

    _SplitMesh(sMNodes, sMEls, numSMNodes, nodeIndMap, nodes, els, sMElIndBounds[0], sMElIndBounds[1], numNodes, trans);
    sMXMLWriter.GetFileStream() << "\t<SubModel>\n";
    wasSuccess &= sMXMLWriter.WriteMesh(sMNodes, numSMNodes, sMEls, sMElIndBounds[1] - sMElIndBounds[0]);
    sMXMLWriter.GetFileStream() << "\t</SubModel>\n";
    wasSuccess &= !sMXMLWriter.GetFileStream().fail();
  }
  wasSuccess &= sMXMLWriter.CloseXML();
  assert(wasSuccess);

  if (wasSuccess) {
    tledModel refModel(const_cast<char*>(refXMLWriter.GetFilePath().c_str()));  
    tledModel sMModel(const_cast<char*>(sMXMLWriter.GetFilePath().c_str()));
    int const *pc_sMNodeInd, *pc_refNodeInd;
  
    tledUnitTestAssert(sMModel.GetMesh()->GetNumNodes() <= refModel.GetMesh()->GetNumNodes());
    tledUnitTestAssert(sMModel.GetMesh()->GetNumEls() == refModel.GetMesh()->GetNumEls());
    for (pc_sMNodeInd = sMModel.GetMesh()->GetAllElNodeInds(), pc_refNodeInd = refModel.GetMesh()->GetAllElNodeInds(); pc_sMNodeInd < sMModel.GetMesh()->GetAllElNodeInds() + 4*sMModel.GetMesh()->GetNumEls(); pc_sMNodeInd++, pc_refNodeInd++) {
      float diff[3];

      tledUnitTestAssert(tledVectorArithmetic::Norm(tledVectorArithmetic::Sub(diff, sMModel.GetMesh()->GetAllNodeCds() + 3*(*pc_sMNodeInd), refModel.GetMesh()->GetAllNodeCds() + 3*(*pc_refNodeInd))));
    }
  } else {
    cerr << "Warning unit test is buggy\n";
  }

  sMXMLWriter.CleanUp();
  refXMLWriter.CleanUp();
} /* _TestSplitMeshModels */

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestSplitMeshModels(tledUnitTest::GetMeshPath("collision_bar.msh"));
  _TestSplitMeshModels(tledUnitTest::GetMeshPath("organic_shape.msh"));
  _TestSplitMeshModels(tledUnitTest::GetMeshPath("sphere.msh"));

  _TestConstraints(tledUnitTest::GetMeshPath("collision_bar.msh"));
  _TestConstraints(tledUnitTest::GetMeshPath("organic_shape.msh"));
  _TestConstraints(tledUnitTest::GetMeshPath("sphere.msh"));

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
