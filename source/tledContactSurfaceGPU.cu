// =========================================================================
// File:       tledContactSurfaceGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledContactSurfaceGPU.h"
#include "tledCUDAMemoryBlock.h"

tledContactSurfaceGPU::~tledContactSurfaceGPU() {  
  if (mdp_NodeCoordinates != NULL) tledCheckCUDAErrors(cudaFree(mdp_NodeCoordinates));

  if (mhp_Surface != NULL) {
    tledCheckCUDAErrors(cudaFree(mhp_Surface->Edges));
    tledCheckCUDAErrors(cudaFreeHost(mhp_Surface));  
  }
  if (mdp_Surface != NULL) tledCheckCUDAErrors(cudaFree(mdp_Surface));

  if (mdp_FacetVertexIndices != NULL) tledCheckCUDAErrors(cudaFree(mdp_FacetVertexIndices));
  if (mdp_FacetEdgeIndices != NULL) tledCheckCUDAErrors(cudaFree(mdp_FacetEdgeIndices));
}

void tledContactSurfaceGPU::InitNodes(const int numNodes) {
  tledCUDAHelpers::AllocateDeviceMemory(mdp_NodeCoordinates, numNodes);
  this->GetHostGPUSurface().NodeCoordinates = mdp_NodeCoordinates;
  this->GetHostGPUSurface().NumberOfNodes = numNodes;
}

void tledContactSurfaceGPU::SetOnDeviceFacetVertexIndices(const int *inds, const int numInds) {
  assert(mdp_FacetVertexIndices == NULL);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_FacetVertexIndices, numInds);
  tledCUDAHelpers::CopyToDevice(mdp_FacetVertexIndices, inds, numInds);
}

void tledContactSurfaceGPU::SetOnDeviceFacetEdgeIndices(const int *inds, const int numInds) {
  assert(mdp_FacetEdgeIndices == NULL);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_FacetEdgeIndices, numInds);
  tledCUDAHelpers::CopyToDevice(mdp_FacetEdgeIndices, inds, numInds);
}

void tledContactSurfaceGPU::InitDeviceSurface() {
  this->AllocateDeviceSurface();
  this->InitNodes();
  this->InitEdges();
  this->InitFacets();
  this->CopyNodes();
  this->CopySurfaceToDevice();
}

