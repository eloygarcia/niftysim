// =========================================================================
// File:       tledShellMesh.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <const int t_numFacetVertices> 
tledShellMesh<t_numFacetVertices>::tledShellMesh(XMLNode xModel) {
  const std::vector<int> facetDefs = GetXMLTextAsVector<int>(xModel.getChildNode("ShellElements"));

#ifndef _GPU_
  assert((t_numFacetVertices != 3 || std::string(xModel.getChildNode("ShellElements").getAttribute("Type")) == "T3") 
	 && (t_numFacetVertices != 4 || std::string(xModel.getChildNode("ShellElements").getAttribute("Type")) == "Q4"));
#endif

  this->SetNumberOfFacets(facetDefs.size()/t_numFacetVertices);
  for (int f = 0; f < (int)(facetDefs.size()/t_numFacetVertices); f++) {
    std::vector<int>::const_iterator ic_facetDefInds = facetDefs.begin() + t_numFacetVertices*f;

    std::copy(ic_facetDefInds, ic_facetDefInds + t_numFacetVertices, this->GetFacet(f).NodeIndices);
  }
}

