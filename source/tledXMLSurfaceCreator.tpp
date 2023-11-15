// =========================================================================
// File:       tledXMLSurfaceCreator.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
template <class TSurface>
void tledXMLSurfaceCreator<TSurface>::XMLImporter::ParseNodeData() {
  const XMLNode nodeNode = this->GetUniqueChild("Nodes", true);
  const int numNodes = this->template GetNumericAttribute<int>(nodeNode, "NumNodes");
  const std::vector<float> nodes = GetXMLTextAsVector<float>(nodeNode);

  this->GetOutput().SetNumberOfNodes(numNodes);
  std::copy(nodes.begin(), nodes.end(), this->GetOutput().GetAllNodeCoordinates());    
}

template <class TSurface>
void tledXMLSurfaceCreator<TSurface>::XMLImporter::ParseFacetData() {
  typedef typename TSurface::Facet __Facet;

  const XMLNode facetNode = this->GetUniqueChild("Elements", true);
  const int numEls = this->template GetNumericAttribute<int>(facetNode, "NumEls");

  std::vector<int> surfFacetDefs = GetXMLTextAsVector<int>(facetNode);

  if (__Facet::NumberOfVertices*numEls != (int)surfFacetDefs.size()) {
    tledLogErrorStream(tledHelper::FatalError() << "Expected " << __Facet::NumberOfVertices*numEls << " facet vertices, but " << int(surfFacetDefs.size()) << " found.");
  }
  
  this->GetOutput().SetNumberOfFacets(numEls);
  for (int f = 0; f < this->GetOutput().GetNumberOfFacets(); f++) {
    std::copy(surfFacetDefs.begin() + f*__Facet::NumberOfVertices, surfFacetDefs.begin() + (f + 1)*__Facet::NumberOfVertices, this->GetOutput().GetFacet(f).NodeIndices);
  }  
}

template <class TSurface>
void tledXMLSurfaceCreator<TSurface>::XMLImporter::Import() {
  this->ParseNodeData();
  this->ParseFacetData();
}

template <class TSurface>
void tledXMLSurfaceCreator<TSurface>::Create() {
  if (pc_RootNode == NULL) {
    tledFatalError("Have no XML spec.");
  } else {
    XMLImporter *p_parser = this->InstantiateParser();

    p_parser->SetRootNode(*pc_RootNode);
    p_parser->SetOuputObject(this->GetOutput());
    p_parser->Import();
    
    delete p_parser;
  }
}
