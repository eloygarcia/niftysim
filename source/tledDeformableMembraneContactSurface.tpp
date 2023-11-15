// =========================================================================
// File:       tledDeformableMembraneContactSurface.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <const int t_numFacetVertices, class TAPI>
void tledDeformableMembraneContactSurfaceImpl<t_numFacetVertices, TAPI>::ConstructFromSolidMeshAndMembrane(const tledMesh &volMesh, const tledSurface &genMembrane) { 
  typedef typename tledShellMesh<t_numFacetVertices>::Facet __MembraneFacet;

  const tledShellMesh<t_numFacetVertices> &membrane = static_cast<const tledShellMesh<t_numFacetVertices>&>(genMembrane);

  std::vector<int> membrane2ContactMap(membrane.GetNumberOfNodes(), -1), contact2MembraneMap;

  this->ConstructFromSolidMesh(volMesh);
  contact2MembraneMap.reserve(membrane.GetNumberOfFacets()*t_numFacetVertices);
  this->SetMembraneNodeBaseIndex(this->GetNumberOfNodes());

  for (typename std::vector<__MembraneFacet>::const_iterator ic_mf = membrane.GetAllFacets().begin(); ic_mf < membrane.GetAllFacets().end(); ic_mf++) {
    for (int const *pc_n = ic_mf->NodeIndices; pc_n < ic_mf->NodeIndices + t_numFacetVertices; pc_n++) contact2MembraneMap.push_back(*pc_n);
  }
  contact2MembraneMap = tledHelper::MakeSortedUnique(contact2MembraneMap);
  for (size_t n = 0; n < contact2MembraneMap.size(); n++) {
    membrane2ContactMap[contact2MembraneMap[n]] = n;
  }

  this->SetNumberOfMembraneNodes(contact2MembraneMap.size());
  this->SetNumberOfNodes(this->GetNumberOfNodes() + 2*this->GetNumberOfMembraneNodes());
  for (int n = 0; n < this->GetNumberOfMembraneNodes(); n++) {
    const float *x = membrane.GetNodeCoordinates(contact2MembraneMap[n]);

    std::copy(x, x + 3, this->GetAllNodeCoordinates() + 3*(this->GetMembraneNodeBaseIndex() + n));
    std::copy(x, x + 3, this->GetAllNodeCoordinates() + 3*(this->GetMembraneNodeBaseIndex() + this->GetNumberOfMembraneNodes() + n));
  }
  
  this->SetMembraneFacetBaseIndex(this->GetNumberOfFacets());
  this->SetNumberOfMembraneElements(membrane.GetNumberOfFacets());
  
  {
    tledSurfaceTopology<tledShellMesh<t_numFacetVertices> > surfaceTopology(membrane);  
    std::vector<int> surf2VolMap = this->GetSurface2VolumeNodeMap();
    std::vector<int> vol2SurfMap = this->GetVolume2SurfaceNodeMap();

    surfaceTopology.ComputeEdges();
    this->SetNumberOfMembraneEdges(surfaceTopology.GetEdges().size());
    this->SetMembraneEdgeBaseIndex(this->GetNumberOfEdges());

    this->GetAllEdges().reserve(this->GetMembraneEdgeBaseIndex() + 2*this->GetNumberOfMembraneEdges());    
    this->SetNumberOfFacets(this->GetNumberOfFacets() + 2*membrane.GetAllFacets().size());
    for (int side = 0; side < 2; side++) {            
      surf2VolMap.insert(surf2VolMap.end(), contact2MembraneMap.begin(), contact2MembraneMap.end());
      if (side == 0) for (std::vector<int>::const_iterator ic_n = contact2MembraneMap.begin(); ic_n < contact2MembraneMap.end(); ic_n++) {
	  if (vol2SurfMap[*ic_n] < 0) vol2SurfMap[*ic_n] = (ic_n - contact2MembraneMap.begin()) + this->GetMembraneNodeBaseIndex();
	}

      for (std::vector<std::pair<int, int> >::const_iterator ic_e = surfaceTopology.GetEdges().begin(); ic_e < surfaceTopology.GetEdges().end(); ic_e++) {
	std::pair<int, int> edge;

	assert(membrane2ContactMap[ic_e->first] >= 0 && membrane2ContactMap[ic_e->second] >= 0);
	edge.first = this->GetMembraneNodeBaseIndex() + side*this->GetNumberOfMembraneNodes() + membrane2ContactMap[ic_e->first];
	edge.second = this->GetMembraneNodeBaseIndex() + side*this->GetNumberOfMembraneNodes() + membrane2ContactMap[ic_e->second];

	assert(edge.first < this->GetNumberOfNodes() && edge.second < this->GetNumberOfNodes());
	
	this->GetAllEdges().push_back(edge);
      }      

      for (int fInd = 0; fInd < this->GetNumberOfMembraneElements(); fInd++) {
	const __MembraneFacet &mf = membrane.GetFacet(fInd);

	Facet &r_newFacet = this->GetFacet(this->GetMembraneFacetBaseIndex() + side*this->GetNumberOfMembraneElements() + fInd);

	for (int vInd = 0; vInd < t_numFacetVertices; vInd++) {
	  r_newFacet.NodeIndices[vInd] = this->GetMembraneNodeBaseIndex() + side*this->GetNumberOfMembraneNodes() + membrane2ContactMap[mf.NodeIndices[vInd]];
	  assert(r_newFacet.NodeIndices[vInd] >= this->GetMembraneNodeBaseIndex() && r_newFacet.NodeIndices[vInd] < this->GetNumberOfNodes());
	}

	for (int eInd = 0; eInd < t_numFacetVertices; eInd++) {
	  r_newFacet.EdgeIndices[eInd] = this->GetMembraneEdgeBaseIndex() + side*this->GetNumberOfMembraneEdges() + surfaceTopology.GetFacetEdges()[fInd][eInd];
	  assert(r_newFacet.EdgeIndices[eInd] >= this->GetMembraneEdgeBaseIndex() && r_newFacet.EdgeIndices[eInd] < this->GetNumberOfEdges());
	}
      }
    } /* for sides */
    this->SetVolume2SurfaceNodeMap(vol2SurfMap);
    this->SetSurface2VolumeNodeMap(surf2VolMap);
  }

  assert(this->GetNumberOfFacets() == this->GetMembraneFacetBaseIndex() + 2*this->GetNumberOfMembraneElements());
  assert(this->GetNumberOfNodes() == this->GetMembraneNodeBaseIndex() + 2*this->GetNumberOfMembraneNodes());
  assert(this->GetNumberOfEdges() == this->GetMembraneEdgeBaseIndex() + 2*this->GetNumberOfMembraneEdges());
}

template <const int t_numFacetVertices, class TAPI>
XMLNode tledDeformableMembraneContactSurfaceImpl<t_numFacetVertices, TAPI>::ExportToXML() const {
  tledDeformableMembraneContactSurfaceXMLExporter<tledDeformableMembraneContactSurfaceImpl> exporter;

  exporter.SetInput(*this);
  exporter.Export();
  
  return exporter.GetRootNode();
}
