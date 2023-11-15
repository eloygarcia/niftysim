// =========================================================================
// File:       tledMovingRigidContactSurface.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <const int t_numFacetVertices, class TAPI> 
float* tledMovingRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::GetFinalDisplacement(float *p_dst) const {
  using namespace tledVectorArithmetic;

  if (this->HasRotation()) {
    float t[3], rc[3], rs[3];

    Add(t, this->GetRotationCentre(), this->GetTotalTranslation());
    for (int i = 0; i < 3; i++) {
      const float a = this->GetAllTotalRotationAngles()[i];

      rc[i] = std::cos(a);
      rs[i] = std::sin(a);      
    }

    for (int n = 0; n < this->GetNumberOfNodes(); n++) {
      float x[3], rTmp[3];

      Sub(x, this->GetNodeCoordinates(n), this->GetRotationCentre());

      for (int i = 0; i < 3; i++) {
	std::copy(x, x + 3, rTmp);
	x[(i+1)%3] = rc[i]*rTmp[(i+1)%3] - rs[i]*rTmp[(i+2)%3];
	x[(i+2)%3] = rc[i]*rTmp[(i+2)%3] + rs[i]*rTmp[(i+1)%3];
      }

      Add(x, x, t);
      Sub(p_dst + 3*n, x, this->GetNodeCoordinates(n));
    }
  } else {
    for (int n = 0; n < this->GetNumberOfNodes(); n++) {
      std::copy(this->GetTotalTranslation(), this->GetTotalTranslation() + 3, p_dst + 3*n);
    }
  }

  return p_dst;
}

template <const int t_numFacetVertices, class TAPI> 
void tledMovingRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::MovingSurfaceXMLImporter::ReadBasicMeshSpec() {
  Superclass::BasicSurfaceXMLImporter::ReadBasicMeshSpec();
}

template <const int t_numFacetVertices, class TAPI> 
void tledMovingRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::MovingSurfaceXMLImporter::ReadMotion() {
  for (int i = 0; i < this->GetRootNode().nChildNode("Motion"); i++) {
    XMLNode mNode = this->GetRootNode().getChildNode("Motion", i);    
    float *p_v;
    
    if (mNode.getAttribute("Type") != NULL && std::string("Translation") == mNode.getAttribute("Type")) {
      p_v = GetXMLTextAsArray<float>(mNode, 3);
      this->GetOutput().SetTranslation(p_v);
    } else if (mNode.getAttribute("Type") != NULL && std::string("Rotation") == mNode.getAttribute("Type")) {
      p_v = GetXMLTextAsArray<float>(mNode, 6);
      this->GetOutput().SetRotations(p_v, p_v + 3);
    }
    
    delete[] p_v;
  } 
 
  if (this->GetRootNode().nChildNode("Motion") == 0) {
    tledFatalError("Have moving rigid surface without \"Motion\" tag?!");    
  }
}

template <const int t_numFacetVertices, class TAPI> 
void tledMovingRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::MovingSurfaceXMLImporter::Import() {
  Superclass::BasicSurfaceXMLImporter::Import();
  this->ReadMotion();
}

