// =========================================================================
// File:       tledDeformableContactSurfaceGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableContactSurfaceGPU_kernels.tpp"

template <class TBaseSurface>
tledDeformableContactSurfaceImplGPU<TBaseSurface>::~tledDeformableContactSurfaceImplGPU() {
  tledCheckCUDAErrors(cudaFree(_GetHostSurface().NodeFacetList));
  tledCheckCUDAErrors(cudaFree(_GetHostSurface().NodeFacetIndexRanges));
  tledCheckCUDAErrors(cudaFree(_GetHostSurface().NodeNormals));
  tledCheckCUDAErrors(cudaFree(_GetHostSurface().NodeMasses));
  tledCheckCUDAErrors(cudaFree(_GetHostSurface().VolumeToSurfaceNodeIndexMap));
  tledCheckCUDAErrors(cudaFree(_GetHostSurface().SurfaceToVolumeNodeIndexMap));
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::AllocateDeviceSurface() {
  tledCUDAHelpers::AllocateDeviceMemory<GPUSurface>(reinterpret_cast<GPUSurface*&>(this->GetDeviceGPUSurface()));
  tledCUDAHelpers::AllocateHostMemory<GPUSurface>(reinterpret_cast<GPUSurface*&>(this->GetHostGPUSurfacePointer()));
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::Init() {
  Superclass::Init();
  m_FacetNormals.resize(3*this->GetNumberOfFacets());
  this->InitDeviceSurface();
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::LoadFromXMLPostloadHook() {  
  Superclass::LoadFromXMLPostloadHook();
  m_FacetNormals.resize(3*this->GetNumberOfFacets());
  this->InitDeviceSurface();
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::InitNodes() {  
  Superclass::InitNodes();

  {
    std::vector<int2> nodeFacetIndexBounds;

    nodeFacetIndexBounds.reserve(this->GetNumberOfNodes());
    for (int nInd = 0; nInd < this->GetNumberOfNodes(); nInd++) {
      nodeFacetIndexBounds.push_back(make_int2(this->GetNodeNormalData(nInd).NodeFacetsStartIndex, this->GetNodeNormalData(nInd).NodeFacetsEndIndex));
    }
    tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().NodeFacetIndexRanges, nodeFacetIndexBounds.size());
    tledCUDAHelpers::CopyToDevice(_GetHostSurface().NodeFacetIndexRanges, nodeFacetIndexBounds);
  }

  tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().NodeFacetList, this->GetNumberOfAllNodeFacetIndices());
  tledCUDAHelpers::CopyToDevice(_GetHostSurface().NodeFacetList, this->GetAllNodeFacetIndices());
  tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().NodeNormals, this->GetNumberOfNodes());

  tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().NodeCoordinates0, this->GetNumberOfNodes());
  this->CopyNodeVectorToDevice(_GetHostSurface().NodeCoordinates0, &this->GetAllNodeCoordinates0().front());

  tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().VolumeToSurfaceNodeIndexMap, this->GetVolume2SurfaceNodeMap().size());
  tledCUDAHelpers::CopyToDevice(_GetHostSurface().VolumeToSurfaceNodeIndexMap, this->GetVolume2SurfaceNodeMap());
  tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().SurfaceToVolumeNodeIndexMap, this->GetNumberOfNodes());
  tledCUDAHelpers::CopyToDevice(_GetHostSurface().SurfaceToVolumeNodeIndexMap, this->GetSurface2VolumeNodeMap());
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::InitNodeMasses(const float globalMasses[]) {
  Superclass::InitNodeMasses(globalMasses);

  tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().NodeMasses, this->GetNumberOfNodes());
  tledCUDAHelpers::CopyToDevice(_GetHostSurface().NodeMasses, this->GetAllSurfaceNodeMasses());
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::InitFacets() {  
  Superclass::InitFacets();
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::CopyNodes() {
  Superclass::CopyNodes();
  this->CopyNodeVectorToDevice(_GetHostSurface().OldNodeCoordinates, this->GetAllOldNodeCoordinates());
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::CopySurfaceToDevice() {  
  tledCUDAHelpers::CopyToDevice(reinterpret_cast<GPUSurface*>(this->GetDeviceGPUSurface()), &_GetHostSurface());
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::Update(const float4 *dp_us) {
  tledLogDebugStream(tledHelper::Info() << "Performing on-device update " << this->GetUpdateCount());
  tledDeformableContactSurfaceGPU::Update(dp_us);
  tledDeformableContactSurfaceGPU::UpdateNodePositions(_GetHostSurface().NodeCoordinates, _GetHostSurface().NodeCoordinates0, dp_us, _GetHostSurface().SurfaceToVolumeNodeIndexMap, this->GetNumberOfNodes());
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::UpdateAllNormals() {
  tledCUDADeviceMemoryBlock &r_facetNormalVtxAngles = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<float3>(this->GetNumberOfFacets()*Facet::NumberOfVertices);

  tledLogDebugStream(tledHelper::Info() << "Client requests exhaustive node-normal update.");

  {
    const int blockSize = 128;
    const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(this->GetNumberOfFacets(), blockSize);

    tledLogDebugStream(tledHelper::Info() << "Computing facet normal x opening angle for " << this->GetNumberOfFacets()*Facet::NumberOfVertices << " facet vertices.");
    assert(r_facetNormalVtxAngles.IsActive());
    tledDeformableContactSurfaceGPU_kernels::UpdateAllFacetNormals<GPUSurface, Facet::NumberOfVertices> <<<numBlks, blockSize>>> (r_facetNormalVtxAngles.GetBuffer<float3>(), _GetDeviceSurface());
  }

  {
    const int blockSize = 128;
    const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(this->GetNumberOfNodes(), blockSize);

    tledLogDebugStream(tledHelper::Info() << "Computing node normal for " << this->GetNumberOfNodes() << " nodes.");
    assert(r_facetNormalVtxAngles.IsActive());
    tledDeformableContactSurfaceGPU_kernels::UpdateAllNodeNormals<GPUSurface, Facet::NumberOfVertices> <<<numBlks, blockSize>>> (_GetDeviceSurface(), r_facetNormalVtxAngles.GetBuffer<float3>());    
  }

  r_facetNormalVtxAngles.ToggleActive();
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::ComputeNodeNormals(const int *dpc_indBuffer, const int numNodes) {
  const int subBlockSize = 4;
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numNodes*subBlockSize, blockSize);

  tledLogDebugStream(tledHelper::Info() << "Client requests node normals for " << numNodes << " nodes.");
  tledDeformableContactSurfaceGPU_kernels::ComputeNodeNormals<GPUSurface, blockSize, subBlockSize, Facet::NumberOfVertices> <<<numBlks, blockSize>>> (_GetDeviceSurface(), dpc_indBuffer, numNodes);  
  tledDeviceSyncDebug;
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplGPU<TBaseSurface>::Save() {
  tledLogDebugStream(tledHelper::Info() << "Performing save " << this->GetSaveCount() << " on device.");
  tledDeformableContactSurface::Save();  
  
  if (this->GetCoordinateHistorySize() > 1) {
    const int blockSize = 256;
    const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(this->GetNumberOfNodes(), blockSize);

    float3 *dp_head, *dp_tail; 

    dp_head = this->GetOnDeviceOldNodeCoordinateBuffer() + this->GetNumberOfNodes()*(this->GetSaveCount()%this->GetCoordinateHistorySize());  
    tledLogDebugStream(tledHelper::Info() << "Copying current configuration node positions to buffer " << this->GetSaveCount()%this->GetCoordinateHistorySize());
    if (this->GetSaveCount() >= this->GetCoordinateHistorySize()) {
      tledLogDebugStream(tledHelper::Info() << "Setting old node coordinate pointer to buffer  " << (this->GetSaveCount() + 1)%this->GetCoordinateHistorySize());
      dp_tail = this->GetOnDeviceOldNodeCoordinateBuffer() + this->GetNumberOfNodes()*((this->GetSaveCount() + 1)%this->GetCoordinateHistorySize());
    } else dp_tail = this->GetOnDeviceOldNodeCoordinateBuffer();
    
    tledDeformableContactSurfaceGPU_kernels::SaveNodeCoordinates <<<numBlks, blockSize>>> (_GetDeviceSurface(), dp_head, dp_tail);
    _GetHostSurface().OldNodeCoordinates = dp_tail;
  } else {
    const int blockSize = 256;
    const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(this->GetNumberOfNodes(), blockSize);

    tledDeformableContactSurfaceGPU_kernels::SaveNodeCoordinates <<<numBlks, blockSize>>> (_GetDeviceSurface(), this->GetOnDeviceOldNodeCoordinateBuffer(), this->GetOnDeviceOldNodeCoordinateBuffer());
  }
  tledDeviceSyncDebug;
}

