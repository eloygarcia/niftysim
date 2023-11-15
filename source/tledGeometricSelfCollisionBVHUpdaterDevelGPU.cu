// =========================================================================
// File:       tledGeometricSelfCollisionBVHUpdaterDevelGPU.cu
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
#ifndef tledGeometricSelfCollisionBVHUpdaterDevelGPU_CU
#define tledGeometricSelfCollisionBVHUpdaterDevelGPU_CU

#include "tledGeometricSelfCollisionBVHUpdaterDevelGPU_kernels.cu"

template <class TBVH>
tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::~tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>() {
  if (mdp_UpdateNodeNodeIndices != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_UpdateNodeNodeIndices));
    tledCheckCUDAErrors(cudaFree(this->GetOnDeviceUpdateNodes()));
  }
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::AllocateOnDeviceUpdateNodes() {
  tledCUDAHelpers::AllocateDeviceMemory(this->GetOnDeviceUpdateNodes(), this->GetUpdateNodes().size());
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::Init() {
  Superclass::Init();
  tledCUDAHelpers::AllocateDeviceMemory(mdp_UpdateNodeNodeIndices, this->GetUpdateNodeIndices().size());
  tledCUDAHelpers::CopyToDevice(mdp_UpdateNodeNodeIndices, this->GetUpdateNodeIndices());

  this->AllocateOnDeviceUpdateNodes();
  this->CopyUpdateNodesToGPU();
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::ComputeUpdateNodeCentroids() {
  const int blockSize = 128;

  tledGeometricSelfCollisionBVHUpdaterGPU_kernels::ComputeUpdateNodeCentroids<GPUUpdateNodeInfo, blockSize> <<<this->GetUpdateNodes().size(), blockSize>>> (this->GetOnDeviceUpdateNodes(), this->GetOnDeviceUpdateNodeIndices(), this->GetMesh().GetHostGPUSurface().NodeCoordinates);
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::ComputeMaxNonTranslationalDisplacement() {
  const int blockSize = 128;

  tledGeometricSelfCollisionBVHUpdaterGPU_kernels::ComputeUpdateNodeCentroids<GPUUpdateNodeInfo, blockSize> <<<this->GetUpdateNodes().size(), blockSize>>> (this->GetOnDeviceUpdateNodes(), this->GetOnDeviceUpdateNodeIndices(), this->GetMesh().GetHostGPUSurface().NodeCoordinates);
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::ComputeMaxNonRigidDisplacement() {
  std::abort();
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::SetOnDeviceUpdateStatuses() {
  const int blockSize = 128;

  //  tledGeometricSelfCollisionBVHUpdaterGPU_kernels::SetUpdateStatuses<GPUUpdateNodeInfo, tledGeometricSelfCollisionBVHUpdater<TBVH>, typename BounBoundingVolume::GPUBV, BoundingVolume::CanRotate> <<<this->GetUpdateNodes().size()/blockSize + (this->GetUpdateNodes().size()%blockSize > 0), blockSize>>> (this->GetOnDeviceUpdateNodes(), this->GetUpdateNodes().size(), this->GetBVH().GetOnDeviceBVs(), this->GetBVH().GetBVMaxDisplacement(), TBVH::VolinoThresholdAngle);
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::_SetUpdateStatusesCPU() {
  tledCUDAHostMemoryBlock &r_block = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<GPUUpdateNodeInfo>(this->GetUpdateNodes().size());
  GPUUpdateNodeInfo const *pc_src = r_block.GetBuffer<GPUUpdateNodeInfo>();

  tledCUDAHelpers::CopyFromDevice(r_block.GetBuffer<GPUUpdateNodeInfo>(), this->GetOnDeviceUpdateNodes());
  for (int u = 0; u < (int)this->GetUpdateNodes().size(); u++) {
    using namespace tledVectorArithmetic;

    const float maxNonRigid = this->GetBVH().GetBV(this->GetUpdateNodes()[u].BVIndex).SubtreeMinH/2*std::tan((TBVH::VolinoThresholdAngle - this->GetUpdateNodes()[u].LastUpdateConeAngle)/2)/(1 + 1e-3f);
    const float bvMaxDisp = this->GetBVH().GetBVMaxDisplacement();
    const GPUUpdateNodeInfo &gpuUN = r_block.GetBuffer<GPUUpdateNodeInfo>()[u];
    
    UpdateNodeInfo &r_un = this->GetUpdateNodes()[u];

    assert(!BoundingVolume::CanRotate);
    this->ConvertFromGPU(r_un, gpuUN);    
    if (Norm(r_un.Translation) + gpuUN.MaxNonTranslational < bvMaxDisp && gpuUN.MaxNonTranslational < maxNonRigid) {
      r_un.Status = tledGeometricSelfCollisionBVHUpdater<TBVH>::INACTIVE;
    } else if (gpuUN.MaxNonTranslational < maxNonRigid && gpuUN.MaxNonTranslational < bvMaxDisp) {
      r_un.Status = tledGeometricSelfCollisionBVHUpdater<TBVH>::TRANSLATION;
    } else r_un.Status = 0;
  }

  r_block.ToggleActive();
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::ComputeNodeUpdateStatuses() {  
  this->ComputeUpdateNodeCentroids();
  this->ComputeMaxNonTranslationalDisplacement();
  if (BoundingVolume::CanRotate) this->ComputeMaxNonRigidDisplacement();
  _SetUpdateStatusesCPU();
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::CopyUpdateNodesToGPU() {
  tledCUDAHostMemoryBlock &r_block = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<GPUUpdateNodeInfo>(this->GetUpdateNodes().size());
  GPUUpdateNodeInfo *p_dst = r_block.GetBuffer<GPUUpdateNodeInfo>();

  for (typename std::vector<UpdateNodeInfo>::iterator i_updateNode = this->GetUpdateNodes().begin(); i_updateNode < this->GetUpdateNodes().end(); i_updateNode++, p_dst++) {
    this->ConvertToGPU(*p_dst, *i_updateNode);
  }

  tledCUDAHelpers::CopyToDevice(this->GetOnDeviceUpdateNodes(), r_block.GetBuffer<GPUUpdateNodeInfo>());
  r_block.ToggleActive();
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::CopyUpdateNodesFromGPU() {
  tledCUDAHostMemoryBlock &r_block = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<GPUUpdateNodeInfo>(this->GetUpdateNodes().size());
  GPUUpdateNodeInfo const *pc_src = r_block.GetBuffer<GPUUpdateNodeInfo>();

  tledCUDAHelpers::CopyFromDevice(r_block.GetBuffer<GPUUpdateNodeInfo>(), this->GetOnDeviceUpdateNodes());
  for (typename std::vector<UpdateNodeInfo>::iterator i_updateNode = this->GetUpdateNodes().begin(); i_updateNode < this->GetUpdateNodes().end(); i_updateNode++, pc_src++) {
    this->ConvertFromGPU(*i_updateNode, *pc_src);
  }
  r_block.ToggleActive();
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::UpdateBVH() {
  if (false) {
    this->ComputeNodeUpdateStatuses();

    this->CopyUpdateNodesFromGPU();
    for (typename std::vector<UpdateNodeInfo>::iterator i_updateNode = this->GetUpdateNodes().begin(); i_updateNode < this->GetUpdateNodes().end(); i_updateNode++) {
      this->PerformSubtreeUpdate(*i_updateNode);
    }

    for (int b = 0; b < this->GetBVH().GetNumberOfBVs(); b++) this->GetBVH().GetBV(b).InitGPU(this->GetBVH().GetHostGPUBV(b));
  } else {
    for (int b = 0; b < this->GetBVH().GetNumberOfBVs(); b++) this->GetBVH().GetBV(b).InitGPU(this->GetBVH().GetHostGPUBV(b));
    Superclass::UpdateBVH();
  }

  this->GetBVH().CopyBVsToGPU();
  this->GetBVH().CopyBVHToGPU();
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::ConvertToGPU(GPUUpdateNodeInfo &r_gpuUN, const typename Superclass::UpdateNodeInfo &hostUN) const {
  for (int c = 0; c < 3; c++) r_gpuUN.Rotation[c] = tledCUDAHelpers::ConvertToFloat3(hostUN.Rotation[c]);
  r_gpuUN.Translation = tledCUDAHelpers::ConvertToFloat3(hostUN.Translation);
  r_gpuUN.CurrentCentroid = tledCUDAHelpers::ConvertToFloat3(hostUN.CurrentCentroid);
  r_gpuUN.LastUpdateCentroid = tledCUDAHelpers::ConvertToFloat3(hostUN.LastUpdateCentroid);
  r_gpuUN.BVIndex = hostUN.BVIndex;
  r_gpuUN.NodeStartIndex = hostUN.NodeStartIndex;
  r_gpuUN.NodeEndIndex = hostUN.NodeEndIndex;
  r_gpuUN.Status = hostUN.Status;
  r_gpuUN.RigidUpdateCounter = hostUN.RigidUpdateCounter;
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::ConvertFromGPU(typename Superclass::UpdateNodeInfo &r_hostUN, const GPUUpdateNodeInfo &gpuUN) const {
  for (int c = 0; c < 3; c++) tledCUDAHelpers::ConvertFromFloatN(r_hostUN.Rotation[c], gpuUN.Rotation[c]);
  tledCUDAHelpers::ConvertFromFloatN(r_hostUN.Translation, gpuUN.Translation);
  tledCUDAHelpers::ConvertFromFloatN(r_hostUN.CurrentCentroid, gpuUN.CurrentCentroid);
  tledCUDAHelpers::ConvertFromFloatN(r_hostUN.LastUpdateCentroid, gpuUN.LastUpdateCentroid);
  r_hostUN.BVIndex = gpuUN.BVIndex;
  r_hostUN.NodeStartIndex = gpuUN.NodeStartIndex;
  r_hostUN.NodeEndIndex = gpuUN.NodeEndIndex;
  r_hostUN.Status = gpuUN.Status;
  r_hostUN.RigidUpdateCounter = gpuUN.RigidUpdateCounter;
}

#endif
