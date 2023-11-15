// =========================================================================
// File:       tledMovingRigidContactSurfaceCPU.tpp
// Purpose:    
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

template <class TSurface> 
void tledMovingRigidContactSurfaceImplCPU<TSurface>::MovingSurfaceXMLImporter::ReadBasicMeshSpec() {
  Superclass::BasicSurfaceXMLImporter::ReadBasicMeshSpec();
  this->GetOutput().SetAllNodeCoordinates0(this->GetOutput().GetAllNodeCoordinates());
  this->GetOutput().SetAllOldNodeCoordinates(this->GetOutput().GetAllNodeCoordinates());
}

template <class TSurface>
void tledMovingRigidContactSurfaceImplCPU<TSurface>::TranslateProjectionOperators() {
  for (int f = 0; f < this->GetNumberOfFacets(); f++) {
    const float *v0 = this->GetNodeCoordinates(this->GetFacet(f).NodeIndices[0]);

    float *p_R = this->GetFacetData(f).ProjectionOperator;

    p_R[3] = -(p_R[0]*v0[0] + p_R[1]*v0[1] + p_R[2]*v0[2]);
    p_R[7] = -(p_R[4]*v0[0] + p_R[5]*v0[1] + p_R[6]*v0[2]);
    p_R[11] = -(p_R[8]*v0[0] + p_R[9]*v0[1] + p_R[10]*v0[2]);
  }
}

template <class TSurface>
void tledMovingRigidContactSurfaceImplCPU<TSurface>::Update() {
  using namespace tledVectorArithmetic;

  Superclass::Update();
  if (this->HasRotation()) {
    float t[3], rc[3], rs[3];
    
    Add(t, this->GetRotationCentre(), this->GetTranslation(t));
    for (int i = 0; i < 3; i++) {
      const float a = this->GetRotationAngle(i);

      rc[i] = std::cos(a);
      rs[i] = std::sin(a);      
    }

    this->ApplyRotationsAndTranslationToNodeList(this->GetAllNodeCoordinates(), rc, rs, t);
    this->ApplyRotationsToNormals(this->GetAllNodeNormals(), this->GetAllNodeNormals0(), this->GetNumberOfNodes(), rc, rs);

    for (int f = 0; f < this->GetNumberOfFacets(); f++) {
      std::copy(this->GetFacetNormal0(f), this->GetFacetNormal0(f) + 3, this->GetFacetData(f).Normal);
      this->RotateVector(this->GetFacetData(f).Normal, rc, rs);
      this->ComputeFacetProjectionOperator(this->GetFacetData(f).ProjectionOperator, f);
    }

    if (this->GetCurrentStep() > this->GetHistoryLength()) {
      Add(t, this->GetRotationCentre(), this->GetTranslation(t, this->GetCurrentStep() - this->GetHistoryLength()));
      for (int i = 0; i < 3; i++) {
	const float a = this->GetRotationAngle(i, this->GetCurrentStep() - this->GetHistoryLength());
	rc[i] = std::cos(a);
	rs[i] = std::sin(a);      
      }
      this->ApplyRotationsAndTranslationToNodeList(this->GetAllOldNodeCoordinates(), rc, rs, t);
      this->ApplyRotationsToNormals(this->GetAllOldFacetNormals(), this->GetAllFacetNormals0(), this->GetNumberOfFacets(), rc, rs);
    }
  } else if (Norm(this->GetTotalTranslation()) > 0) {
    float t[3];

    this->ApplyTranslationToNodeList(this->GetAllNodeCoordinates(), this->GetTranslation(t));
    this->TranslateProjectionOperators();
    
    if (this->GetCurrentStep() > this->GetHistoryLength()) {
      this->ApplyTranslationToNodeList(this->GetAllOldNodeCoordinates(), this->GetTranslation(t, this->GetCurrentStep() - this->GetHistoryLength()));
    }
  }
}

template <class TBaseSurface> 
void tledMovingRigidContactSurfaceImplCPU<TBaseSurface>::ResetNodes() {
  std::copy(this->GetAllNodeCoordinates0(), this->GetAllNodeCoordinates0() + 3*this->GetNumberOfNodes(), this->GetAllNodeCoordinates());
}

template <class TBaseSurface> 
void tledMovingRigidContactSurfaceImplCPU<TBaseSurface>::InitNormals() {
  Superclass::InitNormals();

  std::copy(this->GetAllNodeNormals(), this->GetAllNodeNormals() + 3*this->GetNumberOfNodes(), this->GetAllNodeNormals0());
  for (int f = 0; f < this->GetNumberOfFacets(); f++) {
    std::copy(this->GetFacetData(f).Normal, this->GetFacetData(f).Normal + 3, this->GetAllFacetNormals0() + 3*f);
    assert(std::fabs(1 - tledVectorArithmetic::Norm(this->GetFacetNormal0(f))) < 1e-4f);
  }
  
  std::copy(this->GetAllFacetNormals0(), this->GetAllFacetNormals0() + 3*this->GetNumberOfFacets(), this->GetAllOldFacetNormals());
}
