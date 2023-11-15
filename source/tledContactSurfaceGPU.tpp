// =========================================================================
// File:       tledContactSurfaceGPU.tpp
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
#ifndef tledContactSurfaceGPU_CU
#define tledContactSurfaceGPU_CU

#include "tledContactSurfaceGPU.h"
#include "tledCUDAMemoryBlock.h"

template <class TBaseSurface>
void tledContactSurfaceImplGPU<TBaseSurface>::CopyNodeVectorToDevice(float3 *dp_dst, const float hostV[]) {
  tledCUDAHostMemoryBlock &r_tmp = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<float3>(this->GetNumberOfNodes());
  float3 *p_out = r_tmp.GetBuffer<float3>();

  for (int n = 0; n < this->GetNumberOfNodes(); n++) {
    const float *nCd = hostV + 3*n;

    p_out[n].x = nCd[0];
    p_out[n].y = nCd[1];
    p_out[n].z = nCd[2];
  }
  tledCUDAHelpers::CopyToDevice(dp_dst, p_out, this->GetNumberOfNodes()); 
  r_tmp.ToggleActive();
}

template <class TBaseSurface>
void tledContactSurfaceImplGPU<TBaseSurface>::CopyNodeVectorToHost(float *p_dst, const float3 *dpc_gpuV) {
  tledCUDAHostMemoryBlock &r_tmp = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<float3>(this->GetNumberOfNodes());
  float3 *p_out = r_tmp.GetBuffer<float3>();

  tledCUDAHelpers::CopyFromDevice(p_out, dpc_gpuV, this->GetNumberOfNodes()); 
  for (int n = 0; n < this->GetNumberOfNodes(); n++) {
    float *p_nCd = p_dst + 3*n;

    p_nCd[0] = p_out[n].x;
    p_nCd[1] = p_out[n].y;
    p_nCd[2] = p_out[n].z;
  }
  r_tmp.ToggleActive();  
}

template <class TBaseSurface>
tledContactSurfaceImplGPU<TBaseSurface>::~tledContactSurfaceImplGPU() {
  tledCheckCUDAErrors(cudaFree(_GetHostGPUSurface().Facets));
}

template <class TBaseSurface>
void tledContactSurfaceImplGPU<TBaseSurface>::CopySurfaceToDevice() {  
  tledCUDAHelpers::CopyToDevice(reinterpret_cast<GPUSurface*>(this->GetDeviceGPUSurface()), &_GetHostGPUSurface());
}

template <class TBaseSurface>
void tledContactSurfaceImplGPU<TBaseSurface>::CopyNodes() {
  this->CopyNodeVectorToDevice(this->GetAllOnDeviceNodeCoordinates(), this->GetAllNodeCoordinates());
}

template <class TBaseSurface>
void tledContactSurfaceImplGPU<TBaseSurface>::InitNodes() {
  Superclass::InitNodes(this->GetNumberOfNodes());
}

template <class TBaseSurface>
void tledContactSurfaceImplGPU<TBaseSurface>::AllocateDeviceSurface() {
  tledCUDAHelpers::AllocateDeviceMemory(reinterpret_cast<GPUSurface*&>(this->GetDeviceGPUSurface()));
  tledCUDAHelpers::AllocateHostMemory(reinterpret_cast<GPUSurface*&>(this->GetHostGPUSurface()));  
}

template <class TBaseSurface>
void tledContactSurfaceImplGPU<TBaseSurface>::InitFacets() {
  tledCUDAHelpers::AllocateDeviceMemory(_GetHostGPUSurface().Facets, this->GetNumberOfFacets());  
  tledCUDAHelpers::CopyToDevice(_GetHostGPUSurface().Facets, this->GetAllFacets());
  _GetHostGPUSurface().NumberOfFacets = this->GetNumberOfFacets();

  {
    std::vector<int> inds;

    inds.reserve(Facet::NumberOfVertices*this->GetNumberOfFacets());
    for (typename std::vector<Facet>::const_iterator ic_f = this->GetAllFacets().begin(); ic_f < this->GetAllFacets().end(); ic_f++) {
      inds.insert(inds.end(), ic_f->NodeIndices, ic_f->NodeIndices + Facet::NumberOfVertices);
    }	   
    this->SetOnDeviceFacetVertexIndices(&inds.front(), inds.size());
  }

  {
    std::vector<int> inds;

    inds.reserve(Facet::NumberOfVertices*this->GetNumberOfFacets());
    for (typename std::vector<Facet>::const_iterator ic_f = this->GetAllFacets().begin(); ic_f < this->GetAllFacets().end(); ic_f++) {
      inds.insert(inds.end(), ic_f->EdgeIndices, ic_f->EdgeIndices + Facet::NumberOfVertices);
    }	   
    this->SetOnDeviceFacetEdgeIndices(&inds.front(), inds.size());
  }
}

template <class TBaseSurface>
void tledContactSurfaceImplGPU<TBaseSurface>::InitEdges() {
  tledCUDAHostMemoryBlock &r_gpuEdges = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<int2>(this->GetNumberOfEdges());
  int2 *dp_edges, *hp_edge = r_gpuEdges.GetBuffer<int2>();

  for (std::vector<std::pair<int, int> >::const_iterator ic_e = this->GetAllEdges().begin(); ic_e < this->GetAllEdges().end(); ic_e++) {
    int2 e;

    e.x = ic_e->first;
    e.y = ic_e->second;
    *(hp_edge++) = e;
  }

  tledCUDAHelpers::AllocateDeviceMemory(dp_edges, this->GetNumberOfEdges());
  tledCUDAHelpers::CopyToDevice(dp_edges, r_gpuEdges.GetBuffer<int2>(), this->GetNumberOfEdges());
  r_gpuEdges.ToggleActive();

  _GetHostGPUSurface().Edges = dp_edges;
  _GetHostGPUSurface().NumberOfEdges = this->GetNumberOfEdges();
}

#endif
