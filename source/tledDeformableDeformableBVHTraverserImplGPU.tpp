// =========================================================================
// File:       tledDeformableDeformableBVHTraverserImplGPU.tpp
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

#include "tledDeformableDeformableBVHTraverserGPU_kernels.cu"


template <class TBVH, class TAPI>
tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI>::~tledDeformableDeformableBVHTraverserImplGPU() {
}

template <class TBVH, class TAPI>
const std::vector<int>& tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI>::GetHostStartBVs() const {
  return this->GetMasterBVH().GetNonAdjacentGeometryNodes();
}

template <class TBVH, class TAPI>
void tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI>::InitStartBVPairs() {
  const int bvhOrder = TBVH::BVHOrder;

  std::vector<int2> masterPairs;
  int maxMasterDepth = -1;

  for (std::vector<int>::const_iterator ic_c = this->GetHostStartBVs().begin(); ic_c < this->GetHostStartBVs().end(); ic_c++) {
    const BoundingVolume &bv = this->GetMasterBVH().GetBV(*ic_c);

    std::vector<int> currLeafs;
    std::vector<int2> currPairs;

    for (int c0 = 0; c0 < bvhOrder - 1; c0++) for (int c1 = c0 + 1; c1 < bvhOrder; c1++) {
	if (bvhOrder == 2 || (bv.ChildIndices[c0] >= 0 && bv.ChildIndices[c1] >= 0)) {
	  currLeafs.clear();
	  this->GetMasterBVH().GetSubtreeLeafs(currLeafs, bv.ChildIndices[c1]);	  
	  currPairs = Superclass::GenerateBroadPhaseCheckPairs(currLeafs, bv.ChildIndices[c0]);
	  masterPairs.insert(masterPairs.end(), currPairs.begin(), currPairs.end());
	  tledLogDebugStream(tledHelper::Info() << "Inserting " << int(currPairs.size()) << " leafs to check against master " << bv.ChildIndices[c0]);
	  maxMasterDepth = std::max(maxMasterDepth, this->GetMasterBVH().GetSubtreeMaxDepth(bv.ChildIndices[c0]) + 1);
	}	
      }
  }

  this->SetBroadPhaseSearchSettings(&masterPairs.front(), masterPairs.size(), maxMasterDepth);
}

template <class TBVH, class TAPI>
tledCUDADeviceMemoryBlock*  tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI>::InvokeNodeFacetNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numNodeFacet) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numNodeFacet, blockSize);

  tledCUDADeviceMemoryBlock &r_nfStage1Buffer = tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<NodeFacetInitialProjection>(numNodeFacet);
  tledDeformableDeformableBVHTraverserGPU_kernels::NodeFacetNarrowPhaseStage1<GPUSurface, NodeFacetInitialProjection, blockSize> <<<numBlks, blockSize, 0, stream>>> (r_nfStage1Buffer.template GetBuffer<NodeFacetInitialProjection>(), dp_outItemCounter, _GetGPUSurface(), inItems, numNodeFacet);      

  return &r_nfStage1Buffer;
}

template <class TBVH, class TAPI>
tledCUDADeviceMemoryBlock*  tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI>::InvokeEdgeEdgeNarrowPhaseStage1(int *dp_outItemCounter, cudaStream_t stream, const int2 *inItems, const int numEdgeEdge) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numEdgeEdge, blockSize);

  tledCUDADeviceMemoryBlock *p_eeStage1Buffer = NULL;

  p_eeStage1Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<EdgeEdgeInitialProjection>(numEdgeEdge);
  tledBVHTraverserGPU_kernels::EdgeEdgeNarrowPhaseStage1<blockSize, GPUSurface, GPUSurface, EdgeEdgeInitialProjection> <<<numBlks, blockSize, 0, stream>>> (p_eeStage1Buffer->template GetBuffer<EdgeEdgeInitialProjection>(), dp_outItemCounter, _GetGPUSurface(), _GetGPUSurface(), inItems, numEdgeEdge);        

  return p_eeStage1Buffer;
}

template <class TBVH, class TAPI>
tledCUDADeviceMemoryBlock*  tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI>::InvokeNodeFacetNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const NodeFacetInitialProjection *inItems, const int numNodeFacet) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numNodeFacet, blockSize);

  tledCUDADeviceMemoryBlock *p_nfStage2Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<NodeFacetNarrowPhaseResult>(numNodeFacet);

  tledBVHTraverserGPU_kernels::NodeFacetNarrowPhaseStage2<blockSize, NodeFacetNarrowPhaseResult, NodeFacetInitialProjection, GPUSurface, GPUSurface> <<<numBlks, blockSize, 0, stream>>> (p_nfStage2Buffer->template GetBuffer<NodeFacetNarrowPhaseResult>(), dp_outItemCounter, inItems, numNodeFacet, _GetGPUSurface(), _GetGPUSurface());
  tledDeviceSyncDebug;      

  return p_nfStage2Buffer;
}

template <class TBVH, class TAPI>
tledCUDADeviceMemoryBlock*  tledDeformableDeformableBVHTraverserImplGPU<TBVH, TAPI>::InvokeEdgeEdgeNarrowPhaseStage2(int *dp_outItemCounter, cudaStream_t stream, const EdgeEdgeInitialProjection *inItems, const int numEdgeEdge) {
  const int blockSize = 128;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numEdgeEdge, blockSize);

  tledCUDADeviceMemoryBlock *p_eeStage2Buffer = &tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<EdgeEdgeNarrowPhaseResult>(numEdgeEdge);

  tledDeformableDeformableBVHTraverserGPU_kernels::EdgeEdgeNarrowPhaseStage2<blockSize, EdgeEdgeNarrowPhaseResult, EdgeEdgeInitialProjection, GPUSurface> <<<numBlks, blockSize, 0, stream>>> (p_eeStage2Buffer->template GetBuffer<EdgeEdgeNarrowPhaseResult>(), dp_outItemCounter, inItems, numEdgeEdge, _GetGPUSurface(), this->GetNarrowPhaseMaxDistance());
  tledDeviceSyncDebug;        

  return p_eeStage2Buffer;
}
