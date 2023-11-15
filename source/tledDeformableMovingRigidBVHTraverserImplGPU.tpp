// =========================================================================
// File:       tledDeformableMovingRigidBVHTraverserImplGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledDeformableMovingRigidBVHTraverserGPU_kernels.tpp"

template <class TDBVH, class TRBVH, class TAPI>
tledCUDADeviceMemoryBlock* tledDeformableMovingRigidBVHTraverserImplGPU<TDBVH, TRBVH, TAPI>::InvokeNodeFacetNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numNodeFacet) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numNodeFacet, blockSize);
	
  tledCUDADeviceMemoryBlock *p_nfStage1Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<NodeFacetInitialProjection>(numNodeFacet);

  if (this->DoMaster()) {
    tledDeformableMovingRigidBVHTraverserGPU_kernels::NodeFacetNarrowPhaseStage1<blockSize> <<<numBlks, blockSize, 0, stream>>> (p_nfStage1Buffer->template GetBuffer<NodeFacetInitialProjection>(), dp_outItemCounter, _GetMasterGPUSurface(), _GetSlaveGPUSurface(), inItems, numNodeFacet);      
  } else {
    tledDeformableMovingRigidBVHTraverserGPU_kernels::NodeFacetNarrowPhaseStage1<blockSize> <<<numBlks, blockSize, 0, stream>>> (p_nfStage1Buffer->template GetBuffer<NodeFacetInitialProjection>(), dp_outItemCounter, _GetSlaveGPUSurface(), _GetMasterGPUSurface(), inItems, numNodeFacet);      
  }

  return p_nfStage1Buffer;
}

template <class TDBVH, class TRBVH, class TAPI>
tledCUDADeviceMemoryBlock* tledDeformableMovingRigidBVHTraverserImplGPU<TDBVH, TRBVH, TAPI>::InvokeEdgeEdgeNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const typename Superclass::EdgeEdgeInitialProjection *inItems, const int numEdgeEdge) {
  const int blockSize = 128;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numEdgeEdge, blockSize);

  tledCUDADeviceMemoryBlock *p_eeStage2Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<EdgeEdgeNarrowPhaseResult>(numEdgeEdge);

  if (this->DoMaster()) {
    tledDeformableMovingRigidBVHTraverserGPU_kernels::EdgeEdgeNarrowPhaseStage2<blockSize> <<<numBlocks, blockSize, 0, stream>>> (p_eeStage2Buffer->template GetBuffer<EdgeEdgeNarrowPhaseResult>(), dp_outItemCounter, inItems, numEdgeEdge, _GetMasterGPUSurface(), _GetSlaveGPUSurface(), this->GetNarrowPhaseMaxDistance());
  } else {
    tledDeformableMovingRigidBVHTraverserGPU_kernels::EdgeEdgeNarrowPhaseStage2<blockSize> <<<numBlocks, blockSize, 0, stream>>> (p_eeStage2Buffer->template GetBuffer<EdgeEdgeNarrowPhaseResult>(), dp_outItemCounter, inItems, numEdgeEdge, _GetSlaveGPUSurface(), _GetMasterGPUSurface(), this->GetNarrowPhaseMaxDistance());
  }
  tledDeviceSyncDebug;      

  return p_eeStage2Buffer;
}
