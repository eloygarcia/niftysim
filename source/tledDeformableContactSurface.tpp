// =========================================================================
// File:       tledDeformableContactSurface.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <const int t_numFacetVertices, class TAPI>
typename tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::CachedNormal& tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::GetNodeNormalCache(const int nodeIndex) { 
  NodeNormalData &r_nData = m_NodeNormals[nodeIndex];

  /* Required for certain CUDA SDK/gcc combinations */

  return r_nData.Normal; 
}

template <const int t_numFacetVertices, class TAPI>
const typename tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::CachedNormal& tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::GetNodeNormalCache(const int nodeIndex) const { 
  const NodeNormalData &nData = m_NodeNormals[nodeIndex];

  /* Required for certain CUDA SDK/gcc combinations */

  return nData.Normal; 
}

template <const int t_numFacetVertices, class TAPI>
void tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::Init() {  
  Superclass::Init();

  std::copy(this->GetAllNodeCoordinates(), this->GetAllNodeCoordinates() + 3*this->GetNumberOfNodes(), m_NodeCoordinates0.begin());  
  
  {
    std::vector<std::vector<int> > nodeFacets;

    nodeFacets.resize(this->GetNumberOfNodes());
    for (int fInd = 0; fInd < this->GetNumberOfFacets(); fInd++) {
      const Facet &facet = this->GetFacet(fInd);

      for (int const *pc_n = facet.NodeIndices; pc_n < facet.NodeIndices + Facet::NumberOfVertices; pc_n++) nodeFacets[*pc_n].push_back(fInd);
    } 

    m_NodeFacetInds.clear();
    m_NodeFacetInds.reserve(4*this->GetNumberOfNodes());
    m_NodeNormals.resize(this->GetNumberOfNodes());
    for (int nInd = 0; nInd < this->GetNumberOfNodes(); nInd++) {
      NodeNormalData &r_nd = m_NodeNormals[nInd];
    
      r_nd.NodeFacetsStartIndex = m_NodeFacetInds.size();
      m_NodeFacetInds.insert(m_NodeFacetInds.end(), nodeFacets[nInd].begin(), nodeFacets[nInd].end());
      r_nd.NodeFacetsEndIndex = m_NodeFacetInds.size();
    }
  }
}

template <const int t_numFacetVertices, class TAPI>
void tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::SetNumberOfNodes(const int numNodes) {
  Superclass::SetNumberOfNodes(numNodes);
  m_NodeNormals.resize(numNodes);
  m_NodeCoordinates0.resize(3*numNodes);
}

template <const int t_numFacetVertices, class TAPI>
void tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::ConstructFromSolidMesh(const tledMesh &volMesh) {
  tledContactVolumeSurfaceExtractor<tledDeformableContactSurfaceImpl> extractor(volMesh);

  extractor.SetOutputMesh(*this);
  extractor.Create();
}

template <const int t_numFacetVertices, class TAPI>
void tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::SetNodeFacetIndices(const std::vector<int> &allNodeFacetIndices, const std::vector<int> &numNodeFacets) {
  if ((int)numNodeFacets.size() != this->GetNumberOfNodes()) {
    tledFatalError("Size of numNodeFacets does not match number of nodes.");
  } else {
    int fStart = 0;

    if ((int)m_NodeNormals.size() != this->GetNumberOfNodes()) m_NodeNormals.resize(this->GetNumberOfNodes());

    m_NodeFacetInds = allNodeFacetIndices;
    for (int n = 0; n < this->GetNumberOfNodes(); n++) {
      m_NodeNormals[n].NodeFacetsStartIndex = fStart;
      fStart += numNodeFacets[n];
      m_NodeNormals[n].NodeFacetsEndIndex = fStart;
    }
  }
}

template <const int t_numFacetVertices, class TAPI>
XMLNode tledDeformableContactSurfaceImpl<t_numFacetVertices, TAPI>::ExportToXML() const {
  tledDeformableContactSurfaceXMLExporter<tledDeformableContactSurfaceImpl> exporter;

  exporter.SetInput(*this);
  exporter.Export();

  return exporter.GetRootNode();
}
