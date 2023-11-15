// =========================================================================
// File:       tledDeformableRigidBVHTraverserImplGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    January 2015
// 
// Copyright (c) 2015, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableRigidBVHTraverserGPU_kernels.cu"

template <class TDBVH, class TRBVH, class TAPI>
void tledDeformableRigidBVHTraverserImplGPU<TDBVH, TRBVH, TAPI>::InitStartBVPairs() {  
  std::vector<int2> pairs = this->GenerateBroadPhaseCheckPairs(this->GetSlaveBVH().GetLeafBVIndices(), 0);

  this->SetBroadPhaseSearchSettings(&pairs.front(), this->GetSlaveBVH().GetNumberOfLeafs(), this->GetMasterBVH().GetMaxDepth() + 1);
}

template <class TDBVH, class TRBVH, class TAPI>
tledCUDADeviceMemoryBlock* tledDeformableRigidBVHTraverserImplGPU<TDBVH, TRBVH, TAPI>::InvokeNodeFacetNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numNodeFacet) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numNodeFacet, blockSize);
	
  tledCUDADeviceMemoryBlock *p_nfStage1Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<NodeFacetInitialProjection>(numNodeFacet);

  if (this->DoMaster()) {
    tledDeformableRigidBVHTraverserGPU_kernels::NodeFacetNarrowPhaseStage1RigidMaster<blockSize, MasterGPUSurface, SlaveGPUSurface, NodeFacetInitialProjection> <<<numBlks, blockSize, 0, stream>>> (p_nfStage1Buffer->template GetBuffer<NodeFacetInitialProjection>(), dp_outItemCounter, _GetMasterGPUSurface(), _GetSlaveGPUSurface(), inItems, numNodeFacet);      
  } else {
    tledDeformableRigidBVHTraverserGPU_kernels::NodeFacetNarrowPhaseStage1DeformableMaster<blockSize, SlaveGPUSurface, MasterGPUSurface, NodeFacetInitialProjection> <<<numBlks, blockSize, 0, stream>>> (p_nfStage1Buffer->template GetBuffer<NodeFacetInitialProjection>(), dp_outItemCounter, _GetSlaveGPUSurface(), _GetMasterGPUSurface(), inItems, numNodeFacet);      
  }

  return p_nfStage1Buffer;
}

template <class TDBVH, class TRBVH, class TAPI>
tledCUDADeviceMemoryBlock* tledDeformableRigidBVHTraverserImplGPU<TDBVH, TRBVH, TAPI>::InvokeEdgeEdgeNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numEdgeEdge) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numEdgeEdge, blockSize);

  tledCUDADeviceMemoryBlock *p_eeStage1Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<EdgeEdgeInitialProjection>(numEdgeEdge);

  if (this->DoMaster()) {
    tledBVHTraverserGPU_kernels::EdgeEdgeNarrowPhaseStage1<blockSize, MasterGPUSurface, SlaveGPUSurface, EdgeEdgeInitialProjection> <<<numBlks, blockSize, 0, stream>>> (p_eeStage1Buffer->template GetBuffer<EdgeEdgeInitialProjection>(), dp_outItemCounter, _GetMasterGPUSurface(), _GetSlaveGPUSurface(), inItems, numEdgeEdge);        
  } else {
    tledBVHTraverserGPU_kernels::EdgeEdgeNarrowPhaseStage1<blockSize, SlaveGPUSurface, MasterGPUSurface, EdgeEdgeInitialProjection> <<<numBlks, blockSize, 0, stream>>> (p_eeStage1Buffer->template GetBuffer<EdgeEdgeInitialProjection>(), dp_outItemCounter, _GetSlaveGPUSurface(), _GetMasterGPUSurface(), inItems, numEdgeEdge);        
  }

  return p_eeStage1Buffer;
}

template <class TDBVH, class TRBVH, class TAPI>
tledCUDADeviceMemoryBlock* tledDeformableRigidBVHTraverserImplGPU<TDBVH, TRBVH, TAPI>::InvokeNodeFacetNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const NodeFacetInitialProjection *inItems, const int numNodeFacet) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numNodeFacet, blockSize);

  tledCUDADeviceMemoryBlock *p_nfStage2Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<NodeFacetNarrowPhaseResult>(numNodeFacet);

  if (this->DoMaster()) {
    tledBVHTraverserGPU_kernels::NodeFacetNarrowPhaseStage2<blockSize, NodeFacetNarrowPhaseResult, NodeFacetInitialProjection, MasterGPUSurface, SlaveGPUSurface> <<<numBlks, blockSize, 0, stream>>> (p_nfStage2Buffer->template GetBuffer<NodeFacetNarrowPhaseResult>(), dp_outItemCounter, inItems, numNodeFacet, _GetMasterGPUSurface(), _GetSlaveGPUSurface());
  } else {
    tledBVHTraverserGPU_kernels::NodeFacetNarrowPhaseStage2<blockSize, NodeFacetNarrowPhaseResult, NodeFacetInitialProjection, SlaveGPUSurface, MasterGPUSurface> <<<numBlks, blockSize, 0, stream>>> (p_nfStage2Buffer->template GetBuffer<NodeFacetNarrowPhaseResult>(), dp_outItemCounter, inItems, numNodeFacet, _GetSlaveGPUSurface(), _GetMasterGPUSurface());
  }
  tledDeviceSyncDebug;      

  return p_nfStage2Buffer;
}

template <class TDBVH, class TRBVH, class TAPI>
tledCUDADeviceMemoryBlock* tledDeformableRigidBVHTraverserImplGPU<TDBVH, TRBVH, TAPI>::InvokeEdgeEdgeNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const EdgeEdgeInitialProjection *inItems, const int numEdgeEdge) {
  const int blockSize = 128;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numEdgeEdge, blockSize);

  tledCUDADeviceMemoryBlock *p_eeStage2Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<EdgeEdgeNarrowPhaseResult>(numEdgeEdge);

  if (this->DoMaster()) {
    tledDeformableRigidBVHTraverserGPU_kernels::EdgeEdgeNarrowPhaseStage2RigidMaster<blockSize, EdgeEdgeNarrowPhaseResult, EdgeEdgeInitialProjection, MasterGPUSurface, SlaveGPUSurface> <<<numBlocks, blockSize, 0, stream>>> (p_eeStage2Buffer->template GetBuffer<EdgeEdgeNarrowPhaseResult>(), dp_outItemCounter, inItems, numEdgeEdge, _GetMasterGPUSurface(), _GetSlaveGPUSurface(), this->GetNarrowPhaseMaxDistance());
  } else {
    tledDeformableRigidBVHTraverserGPU_kernels::EdgeEdgeNarrowPhaseStage2DeformableMaster<blockSize, EdgeEdgeNarrowPhaseResult, EdgeEdgeInitialProjection, SlaveGPUSurface, MasterGPUSurface> <<<numBlocks, blockSize, 0, stream>>> (p_eeStage2Buffer->template GetBuffer<EdgeEdgeNarrowPhaseResult>(), dp_outItemCounter, inItems, numEdgeEdge, _GetSlaveGPUSurface(), _GetMasterGPUSurface(), this->GetNarrowPhaseMaxDistance());
  }
  tledDeviceSyncDebug;      

  return p_eeStage2Buffer;
}
