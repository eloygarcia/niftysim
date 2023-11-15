// =========================================================================
// File:       tledMSHMeshLoader.cpp
// Purpose:    
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
#include "tledMSHMeshLoader.h"
#include "tledHelper.h"

#include <string>
#include <vector>
#include <cassert>
#include <fstream>

static void _SkipLine(std::ifstream &r_fin) {
  char c;

  r_fin.unsetf(std::ios_base::skipws);
  do {
    r_fin >> c;
  } while (c != '\n' && !r_fin.eof() && !r_fin.fail());

  r_fin.setf(std::ios_base::skipws);
}

static bool _SkipAheadUntilKeyword(std::ifstream &r_fin, const std::string &keyword) {
  std::string fileString;

  do {
    r_fin >> fileString;
  } while (!r_fin.fail() && !r_fin.eof() && fileString != keyword);

  return fileString == keyword;
}

void tledMSHMeshLoader::ReadFile() {
  const std::string nsMeshType = this->GetMeshType();
  const int mshMeshType = (nsMeshType == "T4" || nsMeshType == "T4ANP")? 4 : (nsMeshType == "H8"? 5 : -1);

  std::ifstream fin(this->GetFilename().c_str());
  std::string tmpstring;

  if (!fin.is_open()) goto tledMSHMeshLoaderParseFail;

  if (_SkipAheadUntilKeyword(fin, "$MeshFormat")) {
    fin >> tmpstring;
  }

  if (!(*tmpstring.c_str() == '2' && *(tmpstring.c_str() + 1) == '.')) {
    tledFatalError(" Only supports ASCII MSH Version 2.X.");
  }

  if (_SkipAheadUntilKeyword(fin, "$Nodes")) {
    int numNodes;
    int nInd, tmp_nInd;
    std::vector<float> nodes;

    fin >> numNodes;
    nodes.reserve(3*numNodes);
    for (nInd = 0; nInd < numNodes && !fin.fail() && !fin.eof(); nInd++) {
      int cInd;
      float tmpCd;

      fin >> tmp_nInd;

      if (!(tmp_nInd == nInd + 1)) goto tledMSHMeshLoaderParseFail;
      for (cInd = 0; cInd < 3; cInd++) fin >> tmpCd, nodes.push_back(tmpCd);
    }    

    if ((int)nodes.size() != 3*numNodes) goto tledMSHMeshLoaderParseFail;      
    
    this->GetOutput().SetNumberOfNodes(numNodes, 3);
    std::copy(nodes.begin(), nodes.end(), this->GetOutput().GetAllNodeCds());
  } else goto tledMSHMeshLoaderParseFail;      

  if (_SkipAheadUntilKeyword(fin, "$Elements")) {
    int numElements, numRealElements, eInd;
    std::vector<int> elements;

    /*
     * Note: MSH format allows for surface and volume elements within the same element block
     * -> need to count tetras!
     */
    fin >> numElements;
    elements.reserve(numElements);
    numRealElements = 0;
    for (eInd = 0; eInd < numElements && !fin.eof() && !fin.fail(); eInd++) {
      int elemtype, numLabels, lInd, vInd, tmp;

      fin >> tmp;
      assert(tmp == eInd + 1);
      fin >> elemtype;
      if (elemtype == mshMeshType) {	
	fin >> numLabels;
	for (lInd = 0; lInd < numLabels; lInd++) fin >> tmp;
	for (vInd = 0; vInd < this->GetNumberOfElementNodes(); vInd++) {
	  int meshNodeInd;

	  fin >> meshNodeInd;
	  elements.push_back(meshNodeInd - 1);
	}
	numRealElements += 1;
      } else _SkipLine(fin);
    } /* for elements */

    if (eInd == numElements && numRealElements) {
      assert(numRealElements == (int)elements.size()/this->GetNumberOfElementNodes());
      this->GetOutput().SetNumberOfElements(numRealElements, GetMeshType());
      std::copy(elements.begin(), elements.end(), GetOutput().GetAllElNodeInds());
    } else goto tledMSHMeshLoaderParseFail;      
  } else goto tledMSHMeshLoaderParseFail;  

  return;

tledMSHMeshLoaderParseFail:
  tledLogErrorStream(tledHelper::FatalError() << "Reading MSH file " << GetFilename() << " failed.");
}
