// =========================================================================
// File:       tledMeshSurface.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledMeshSurface.h"
#include "tledVolumeSurfaceExtractor.h"

#include <cassert>

template <const int t_numFacetVertices>
class tledMeshSurface<t_numFacetVertices>::_Extractor : public tledVolumeSurfaceExtractor<tledMeshSurface<t_numFacetVertices> > {
public:
  typedef tledVolumeSurfaceExtractor<tledMeshSurface<t_numFacetVertices> > Superclass;

protected:
  virtual void InitNodes() {
    this->GetOutput().SetNodeVector(this->GetMesh().GetAllNodeCds(), this->GetMesh().GetNumNodes());
  }

public:
  _Extractor(const tledMesh &mesh) : Superclass(mesh) {}
  virtual ~_Extractor(void) {}
};

template <>
tledMeshSurface<3>::tledMeshSurface(const tledMesh &mesh) {
  _Extractor extractor(mesh);
  
  extractor.SetOutputMesh(*this);
  extractor.Create();
}

template <>
tledMeshSurface<4>::tledMeshSurface(const tledMesh &mesh) {
  _Extractor extractor(mesh);
  
  extractor.SetOutputMesh(*this);
  extractor.Create();
}
