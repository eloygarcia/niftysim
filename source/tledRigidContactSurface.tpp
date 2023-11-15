// =========================================================================
// File:       tledRigidContactSurface.tpp
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

template <const int t_numFacetVertices, class TAPI>
void tledRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::SetNumberOfFacets(const int numFacets) {
  Superclass::SetNumberOfFacets(numFacets);
}

template <const int t_numFacetVertices, class TAPI>
void tledRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::ComputeNodeNormals(float *p_dst) const {
  using namespace tledVectorArithmetic;

  for (typename std::vector<Facet>::const_iterator ic_facet = this->GetAllFacets().begin(); ic_facet < this->GetAllFacets().end(); ic_facet++) {
    const int numVtcs = Facet::NumberOfVertices;

    for (int nInd = 0; nInd < numVtcs; nInd++) {
      float angle, x[3], y[3], n[3];       

      assert(ic_facet->NodeIndices[nInd] < this->GetNumberOfNodes());
      Sub(x, this->GetNodeCoordinates(ic_facet->NodeIndices[(nInd+1)%numVtcs]), this->GetNodeCoordinates(ic_facet->NodeIndices[nInd]));
      Sub(y, this->GetNodeCoordinates(ic_facet->NodeIndices[(nInd-1+numVtcs)%numVtcs]), this->GetNodeCoordinates(ic_facet->NodeIndices[nInd]));

      Cross(n, x, y);      
      angle = ComputeAngle(x, y);
      assert(angle >= 0 && angle <= tledPi);

      Add(p_dst + 3*ic_facet->NodeIndices[nInd], p_dst + 3*ic_facet->NodeIndices[nInd], ScalarMul(x, n, angle));
    }
  } /* for facets */

  for (float *p_nn = p_dst; p_nn < p_dst + 3*this->GetNumberOfNodes(); p_nn += 3) ScalarDiv(p_nn, Norm(p_nn));
}

template <const int t_numFacetVertices, class TAPI>
void tledRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::BasicSurfaceXMLImporter::ReadBasicMeshSpec() {
  if (this->GetRootNode().nChildNode("VTKSurface")) {
    XMLNode xNode = this->GetRootNode().getChildNode("VTKSurface");
    tledVTKSurfaceLoader<tledRigidContactSurfaceImpl> loader;
    tledSurfaceLoaderSurfaceCreatorAdapter<tledRigidContactSurfaceImpl> loaderAdapter(loader);
    tledContactSurfaceCreator<tledRigidContactSurfaceImpl> creator(loaderAdapter);
    
    loader.SetFilename(xNode.getText());
    if (xNode.nChildNode("Translation")) {	
      float *p_trans;

      p_trans = GetXMLTextAsArray<float>(xNode.getChildNode("Translation"), 3);
      loader.SetTranslation(p_trans);
      
      delete[] p_trans;
    } 

    if (xNode.nChildNode("Rotation")) {
      float *p_rot;

      p_rot = GetXMLTextAsArray<float>(xNode.getChildNode("Rotation"), 6);
      loader.SetRotations(p_rot, p_rot[3], p_rot[4], p_rot[5]);

      delete[] p_rot;
    }

    if (xNode.nChildNode("ScaleFactor")) {
      loader.SetScaleFactor(this->template GetNumericElementValue<float>(xNode.getChildNode("ScaleFactor")));
    }

    creator.SetOutputMesh(this->GetOutput());
    creator.Create();
  } else {
    tledXMLSurfaceCreator<tledRigidContactSurfaceImpl> xmlLoader;
    tledContactSurfaceCreator<tledRigidContactSurfaceImpl> creator(xmlLoader);

    xmlLoader.SetXMLRoot(this->GetRootNode());
    creator.SetOutputMesh(this->GetOutput());
    creator.Create();
  }  

  if (this->GetRootNode().nChildNode("FrictionCoefficient")) {
    this->GetOutput().SetFrictionCoefficient(this->template GetNumericElementValue<float>(this->GetRootNode().getChildNode("FrictionCoefficient")));
  }
}

template <const int t_numFacetVertices, class TAPI>
void tledRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::BasicSurfaceXMLImporter::Import() {
  this->ReadBasicMeshSpec();
  this->GetOutput().Init();
}
