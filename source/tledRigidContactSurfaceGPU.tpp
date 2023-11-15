// =========================================================================
// File:       tledRigidContactSurfaceGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseSurface>
void tledRigidContactSurfaceImplGPU<TBaseSurface>::AllocateDeviceSurface() {
  tledCUDAHelpers::AllocateDeviceMemory<GPUSurface>(reinterpret_cast<GPUSurface*&>(this->GetDeviceGPUSurface()));
  tledCUDAHelpers::AllocateHostMemory<GPUSurface>(reinterpret_cast<GPUSurface*&>(this->GetHostGPUSurfacePointer()));
}

template <class TBaseSurface>
void tledRigidContactSurfaceImplGPU<TBaseSurface>::InitNodes() {  
  Superclass::InitNodes();

  {
    std::vector<float> nodeNormals(3*this->GetNumberOfNodes());
    tledCUDAHostMemoryBlock &r_hostNN = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<float3>(this->GetNumberOfNodes());  
    
    this->ComputeNodeNormals(&nodeNormals.front());
    for (int nInd = 0; nInd < this->GetNumberOfNodes(); nInd++) {
      r_hostNN.GetBuffer<float3>()[nInd] = tledCUDAHelpers::ConvertToFloat3(&nodeNormals[3*nInd]);
    }

    tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().NodeNormals, this->GetNumberOfNodes());
    tledCUDAHelpers::CopyToDevice(_GetHostSurface().NodeNormals, r_hostNN.GetBuffer<float3>(), this->GetNumberOfNodes());

    r_hostNN.ToggleActive();
  }
}

template <class TBaseSurface>
void tledRigidContactSurfaceImplGPU<TBaseSurface>::InitFacets() {  
  Superclass::InitFacets();

  {
    tledCUDAHostMemoryBlock &r_hostFN = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<float3>(this->GetNumberOfFacets());  
    
    for (int fInd = 0; fInd < this->GetNumberOfFacets(); fInd++) {
      float n[3];
      
      r_hostFN.GetBuffer<float3>()[fInd] = tledCUDAHelpers::ConvertToFloat3(this->ComputeNormalisedFacetNormal(n, fInd));
    }

    tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().FacetNormals, this->GetNumberOfFacets());
    tledCUDAHelpers::CopyToDevice(_GetHostSurface().FacetNormals, r_hostFN.GetBuffer<float3>(), this->GetNumberOfFacets());

    r_hostFN.ToggleActive();
  }

  {
    tledCUDAHostMemoryBlock &r_hostPO = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<float4>(3*this->GetNumberOfFacets());  

    for (int fInd = 0; fInd < this->GetNumberOfFacets(); fInd++) {
      float cpuProjOp[12];
      
      this->ComputeFacetProjectionOperator(cpuProjOp, fInd);
      for (int c = 0; c < 3; c++) {
	r_hostPO.GetBuffer<float4>()[3*fInd] = make_float4(cpuProjOp[4*c], cpuProjOp[4*c+1], cpuProjOp[4*c+2], cpuProjOp[4*c+2]);
      }
    }

    tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().FacetProjectionOperators, 3*this->GetNumberOfFacets());
    tledCUDAHelpers::CopyToDevice(_GetHostSurface().FacetProjectionOperators, r_hostPO.GetBuffer<float4>(), 3*this->GetNumberOfFacets());

    r_hostPO.ToggleActive();
  }
}

template <class TBaseSurface>
void tledRigidContactSurfaceImplGPU<TBaseSurface>::CopySurfaceToDevice() {  
  tledCUDAHelpers::CopyToDevice(reinterpret_cast<GPUSurface*>(this->GetDeviceGPUSurface()), &_GetHostSurface());
}

template <class TBaseSurface>
void tledRigidContactSurfaceImplGPU<TBaseSurface>::Init() {
  Superclass::Init();
  this->InitDeviceSurface();
  assert(tledCUDAHostMemoryBlock::GetNumberOfActiveBuffers() == 0);
}
