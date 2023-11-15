// =========================================================================
// File:       tledBottomUpBVHUpdaterGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBottomUpBVHUpdaterGPU_TPP
#define tledBottomUpBVHUpdaterGPU_TPP

#include "tledBottomUpBVHUpdaterGPU.h"

#include <thrust/sort.h>
#include <thrust/unique.h>

#include "tledBottomUpBVHUpdaterGPU_kernels.tpp"
#include "tledAABB_kernels.cu"

template <class TBVH>
void tledBottomUpBVHUpdaterGPU<TBVH>::_ProgressLevel(std::vector<bool> &r_isUpdated, std::vector<int> &r_nextLevel, const std::vector<int> &levels, const int2 currBounds) {
  r_nextLevel.clear();
  for (std::vector<int>::const_iterator ic_b = levels.begin() + currBounds.x; ic_b < levels.begin() + currBounds.y; ic_b++) {
    assert(*ic_b >= 0 && *ic_b < int(r_isUpdated.size()));
    assert(!r_isUpdated[*ic_b]);
    r_isUpdated[*ic_b] = true;
    r_nextLevel.push_back(this->GetBVH().GetBV(*ic_b).ParentIndex);
  }
  r_nextLevel = tledHelper::MakeSortedUnique(r_nextLevel);
}

template <class TBVH>
std::vector<int> tledBottomUpBVHUpdaterGPU<TBVH>::_FilterNonUpdateable(const std::vector<int> &potNextLevel, const std::vector<bool> &isUpdated) {
    std::vector<int> nextLevel;

    nextLevel.reserve(potNextLevel.size());
    for (std::vector<int>::const_iterator ic_b = potNextLevel.begin(); ic_b < potNextLevel.end(); ic_b++) {
      int numUpdated = 0, numChildren = 0;

      for (int const *pc_c = this->GetBVH().GetBV(*ic_b).ChildIndices; pc_c < this->GetBVH().GetBV(*ic_b).ChildIndices + BVH::BVHOrder; pc_c++) {
	if (*pc_c >= 0) {
	  numUpdated += isUpdated[*pc_c];
	  numChildren += 1;
	}
      }

      if (numUpdated == numChildren) {
	nextLevel.push_back(*ic_b);
      }
    }

    return nextLevel;
}

template <class TBVH>
void tledBottomUpBVHUpdaterGPU<TBVH>::Init() {
  std::vector<int> levels(this->GetStartIndices(), this->GetStartIndices() + this->GetNumberOfStartIndices()), nextLevel;
  std::vector<bool> isUpdated(this->GetBVH().GetNumberOfBVs(), false);
  int2 prevLevelBounds;
  
  m_NumLevels = 1;
  m_LevelBounds.clear();

  prevLevelBounds.x = 0;
  prevLevelBounds.y = levels.size();
  m_LevelBounds.push_back(prevLevelBounds);
  nextLevel.reserve(int(this->GetNumberOfStartIndices()/TBVH::BVHOrder*1.25f));
  _ProgressLevel(isUpdated, nextLevel, levels, prevLevelBounds);
  while (nextLevel.size() > 1 || nextLevel[0] != this->GetSubtreeRootParent()) {
    nextLevel = _FilterNonUpdateable(nextLevel, isUpdated);
    prevLevelBounds.x = levels.size();
    levels.insert(levels.end(), nextLevel.begin(), nextLevel.end());
    prevLevelBounds.y = levels.size();
    m_LevelBounds.push_back(prevLevelBounds);
    _ProgressLevel(isUpdated, nextLevel, levels, prevLevelBounds);
    m_NumLevels += 1;

    assert(nextLevel.size() >= 1);    
  }
  assert(m_NumLevels == int(m_LevelBounds.size()));
  assert(m_LevelBounds[0].y - m_LevelBounds[0].x == this->GetNumberOfStartIndices());

  tledCUDAHelpers::AllocateDeviceMemory(mdp_Levels, levels.size());
  tledCUDAHelpers::AllocateDeviceMemory(mdp_LevelBounds, m_NumLevels);  

  tledCUDAHelpers::CopyToDevice(mdp_Levels, levels);
  tledCUDAHelpers::CopyToDevice(mdp_LevelBounds, m_LevelBounds);
}

template <class TBVH>
tledBottomUpBVHUpdaterGPU<TBVH>::~tledBottomUpBVHUpdaterGPU() {
  tledCheckCUDAErrors(cudaFree(mdp_Levels));
  tledCheckCUDAErrors(cudaFree(mdp_LevelBounds));
}

template <class TBVH>
void tledBottomUpBVHUpdaterGPU<TBVH>::SetStartIndices(const int *startIndices, const int numIndices) { 
  mpc_StartIndices = startIndices; 
  m_NumStartIndices = numIndices;
}

template <class TBVH>
void tledBottomUpBVHUpdaterGPU<TBVH>::_RefitInteriorLevel(const int levelInd) {
  typedef typename ContactMesh::GPUSurface __GPUMesh;
  typedef typename BoundingVolume::GPUBV __GPUBV;

  const int blockSize = 256;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(this->GetLevelBounds(levelInd).y - this->GetLevelBounds(levelInd).x, blockSize);

  tledBottomUpBVHUpdaterGPU_kernels::RefitInteriorLevelKernel<__GPUBV> <<<numBlocks, blockSize>>> (this->GetBVH().GetOnDeviceBVs(), this->GetOnDeviceLevels(), this->GetOnDeviceLevelIndexBounds(), levelInd);
  tledDeviceSyncDebug;
}

template <class TBVH>
void tledBottomUpBVHUpdaterGPU<TBVH>::_RefitLeafLevel() {
  typedef typename ContactMesh::GPUSurface __GPUMesh;
  typedef typename BoundingVolume::GPUBV __GPUBV;

  const int blockSize = 256;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(this->GetNumberOfStartIndices(), blockSize);
  const __GPUMesh *dpc_mesh = static_cast<const __GPUMesh*>(this->GetMesh().GetDeviceGPUSurface());

  tledBottomUpBVHUpdaterGPU_kernels::RefitLeafLevelKernel<__GPUBV, __GPUMesh> <<<numBlocks, blockSize>>> (this->GetBVH().GetOnDeviceBVs(), this->GetOnDeviceLevels(), this->GetOnDeviceLevelIndexBounds(), dpc_mesh, 2*this->GetBVH().GetMargin());
  tledDeviceSyncDebug;
}

template <class TBVH>
void tledBottomUpBVHUpdaterGPU<TBVH>::UpdateBVH() {
  _RefitLeafLevel();
  for (int l = 1; l < m_NumLevels; l++) {
    _RefitInteriorLevel(l);
  }  
}

#endif
