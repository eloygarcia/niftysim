// =========================================================================
// File:       tledParallelSelfCollisionBVHCreator.tpp
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

template <class TBVH>
void tledParallelSelfCollisionBVHCreator<TBVH>::BuilderWorker(tledParallelSelfCollisionBVHCreator &r_builder, std::vector<int> &r_primitiveIndices, const int threadSetTargetSize) {
  r_builder.InitBVHStruct();
  if ((int)r_primitiveIndices.size() <= threadSetTargetSize) {
    r_builder.InitialiseTopDownRecursive(0, r_primitiveIndices.begin(), r_primitiveIndices.end());
  } else r_builder.SplitSetInitialiseTopDownRecursive(r_primitiveIndices, threadSetTargetSize);
}

template <class TBVH>
void tledParallelSelfCollisionBVHCreator<TBVH>::_InitThread(BVH &r_bvh, tledParallelSelfCollisionBVHCreator &r_builder) {
  BoundingVolume rootBV;

  rootBV.PrimitiveIndex = rootBV.ParentIndex = -1;
  for (int *p_c = rootBV.ChildIndices; p_c < rootBV.ChildIndices + BoundingVolume::NumberOfChildBVs; p_c++) *p_c = -1;

  r_bvh.GetBVs().push_back(rootBV);
  r_bvh.SetMargin(this->GetOutput().GetMargin());

  r_builder.SetBVHBuffer(r_bvh);
  r_builder.GetBVPrimitiveIndices() = this->GetBVPrimitiveIndices();
  r_builder.GetPrimitiveCentroids() = this->GetPrimitiveCentroids();
  r_builder.m_PrimitiveBVs = this->m_PrimitiveBVs;
}

template <class TBVH>
void tledParallelSelfCollisionBVHCreator<TBVH>::SplitSetInitialiseTopDownRecursive(std::vector<int> &r_primitiveIndices, const int threadSetTargetSize) {
  std::vector<int>::iterator i_childpInds[TBVH::BVHOrder][2];
  std::vector<tledParallelSelfCollisionBVHCreator*> vp_subCreators;
  std::vector<BVH*> vp_subBVHs;
  std::vector<boost::thread*> vp_threads;
  int newChildInd;

  assert(this->GetOutput().GetNumberOfBVs() == 1);
  assert((int)r_primitiveIndices.size() >= TBVH::BVHOrder);

  this->GetOutput().GetBV(0).PrimitiveIndex = -1;

  this->GetOutput().ComputeBoundsFromPrimitives(this->GetOutput().GetBV(0), &r_primitiveIndices.front(), &r_primitiveIndices.back());
  this->SplitBV(i_childpInds, 0, r_primitiveIndices.begin(), r_primitiveIndices.end());    

  for (int cInd = 1; cInd < TBVH::BVHOrder; cInd++) {
    vp_subBVHs.push_back(new BVH(this->GetMesh()));
    vp_subCreators.push_back(new tledParallelSelfCollisionBVHCreator());

    _InitThread(*vp_subBVHs.back(), *vp_subCreators.back());
    vp_threads.push_back(new boost::thread(BuilderWorker, *vp_subCreators.back(), std::vector<int>(i_childpInds[cInd][0], i_childpInds[cInd][1]), threadSetTargetSize));          
  }

  newChildInd = this->AddBV(0);
  this->InitialiseTopDownRecursive(this->GetOutput().GetBV(0).ChildIndices[0] = newChildInd, i_childpInds[0][0], i_childpInds[0][1]);

  for (int cInd = 1; cInd < TBVH::BVHOrder; cInd++) {
    vp_threads[cInd-1]->join();
    delete vp_subCreators[cInd-1];
    delete vp_threads[cInd-1];

    newChildInd = this->AddBV(0);     
    _CopyAttributes(*vp_subBVHs[cInd-1]);
    _InsertSubBVH(*vp_subBVHs[cInd-1], 0, this->GetOutput().GetBV(0).ChildIndices[cInd] = newChildInd);	
    delete vp_subBVHs[cInd-1];	        
  }

  this->GetOutput().RefitInteriorBV(0);
}

template <class TBVH>
void tledParallelSelfCollisionBVHCreator<TBVH>::_InsertSubBVH(const BVH &subBVH, const int rootBVIndex, const int parentIndex) {
  {
    const int oldParent = this->GetOutput().GetBV(parentIndex).ParentIndex;

    this->GetOutput().GetBV(parentIndex) = subBVH.GetBV(rootBVIndex);
    this->GetOutput().GetBV(parentIndex).ParentIndex = oldParent;
  }

  assert(subBVH.GetBVs()[rootBVIndex].PrimitiveIndex < 0);
  for (int c = 0; c < BoundingVolume::NumberOfChildBVs; c++) {
    const int cInd = subBVH.GetBVs()[rootBVIndex].ChildIndices[c];
    const BoundingVolume &child = subBVH.GetBVs()[cInd];
    const int dstInd = this->GetOutput().GetBVs().size();

    this->GetOutput().GetBVs().push_back(child);
    this->GetOutput().GetBVs().back().ParentIndex = parentIndex;
    this->GetOutput().GetBVs()[parentIndex].ChildIndices[c] = dstInd;
    if (this->GetOutput().GetBVs().back().VolinoAngle >= BoundingVolume::GetVolinoThresholdAngle()) this->GetOutput().AddSelfCollisionCandidate(dstInd);
    assert(this->GetOutput().GetBVs().back().HasDisconnectedGeometry() || (this->GetOutput().GetBVs().back().VolinoAngle >= 0 && this->GetOutput().GetBVs().back().VolinoAngle <= 2*tledPi));

    if (child.PrimitiveIndex < 0) _InsertSubBVH(subBVH, cInd, dstInd);
  }

#ifndef NDEBUG
  assert(parentIndex + 1 == this->GetOutput().GetBV(parentIndex).ChildIndices[0] || this->GetOutput().GetBV(parentIndex).ParentIndex < 0 || this->GetOutput().GetBV(this->GetOutput().GetBV(parentIndex).ParentIndex).HasDisconnectedGeometry());
  for (int const *pc_c = this->GetOutput().GetBV(parentIndex).ChildIndices; pc_c < this->GetOutput().GetBV(parentIndex).ChildIndices + TBVH::BVHOrder; pc_c++) {
    assert(this->GetOutput().GetBV(*pc_c).ParentIndex  == parentIndex);
  }
#endif
}

template <class TBVH>
void tledParallelSelfCollisionBVHCreator<TBVH>::_CopyAttributes(const BVH &subBVH) {
  const int offBVInd = this->GetOutput().GetNumberOfBVs() - 1;

  this->GetOutput().GetSelfCollisionCandidates().reserve(this->GetOutput().GetSelfCollisionCandidates().size() + subBVH.GetSelfCollisionCandidates().size());
  for (std::vector<int>::const_iterator ic_c = subBVH.GetSelfCollisionCandidates().begin(); ic_c < subBVH.GetSelfCollisionCandidates().end(); ic_c++) {
    this->GetOutput().GetSelfCollisionCandidates().push_back(*ic_c + offBVInd);
  }  

  this->GetOutput().GetLeafBVIndices().reserve(this->GetOutput().GetLeafBVIndices().size() + subBVH.GetLeafBVIndices().size());
  for (std::vector<int>::const_iterator ic_c = subBVH.GetLeafBVIndices().begin(); ic_c < subBVH.GetLeafBVIndices().end(); ic_c++) {
    this->GetOutput().GetLeafBVIndices().push_back(*ic_c + offBVInd);
    this->GetOutput().GetAllPrimitiveBVIndices()[subBVH.GetBV(*ic_c).PrimitiveIndex] = *ic_c + offBVInd;
  }
}

template <class TBVH>
void tledParallelSelfCollisionBVHCreator<TBVH>::GenerateMain() {
  if (this->GetMesh().GetNumberOfFacets() < 16*(int)boost::thread::hardware_concurrency()) Superclass::GenerateMain();
  else {
    std::vector<std::vector<int> > primitiveSets;

    this->InitialiseClusters();

    this->GetOutput().GetGeometryClusterSubtreeRootIndices() = this->GetOutput().GetLeafBVIndices();
    this->GetOutput().GetLeafBVIndices().clear();

    for (std::vector<int>::const_iterator ic_clusterBVInd = this->GetOutput().GetGeometryClusterSubtreeRootIndices().begin(); ic_clusterBVInd < this->GetOutput().GetGeometryClusterSubtreeRootIndices().end(); ic_clusterBVInd++) {
      const BoundingVolume &bv = this->GetOutput().GetBVs()[*ic_clusterBVInd];

      if (bv.PrimitiveIndex + 1 >= (int)this->m_ClusterStartInds.size()) {
	primitiveSets.push_back(std::vector<int>(this->GetBVPrimitiveIndices().begin() + this->m_ClusterStartInds[bv.PrimitiveIndex], this->GetBVPrimitiveIndices().end()));
      } else {
	primitiveSets.push_back(std::vector<int>(this->GetBVPrimitiveIndices().begin() + this->m_ClusterStartInds[bv.PrimitiveIndex], this->GetBVPrimitiveIndices().begin() + this->m_ClusterStartInds[bv.PrimitiveIndex+1]));
      }
    }

    for (int numProcessed = 0; numProcessed < (int)primitiveSets.size();) {
      std::vector<BVH*> vp_subBVHs;
      std::vector<tledParallelSelfCollisionBVHCreator*> vp_subCreators;
      std::vector<boost::thread*> vp_threads;
      int targetSize = 0, numAvailable = boost::thread::hardware_concurrency();

      for (int i = 0; i < numAvailable && numProcessed + i < (int)primitiveSets.size(); i++) targetSize += primitiveSets[numProcessed+i].size();
      targetSize /= numAvailable;

      for (int i = 0; i < numAvailable && numProcessed + i < (int)primitiveSets.size(); i++) {
	vp_subBVHs.push_back(new BVH(GetMesh()));
	vp_subCreators.push_back(new tledParallelSelfCollisionBVHCreator());
	_InitThread(*vp_subBVHs.back(), *vp_subCreators.back());     
	vp_threads.push_back(new boost::thread(BuilderWorker, *vp_subCreators.back(), primitiveSets[numProcessed+i], targetSize));      
	numAvailable -= primitiveSets[numProcessed+i].size()/boost::thread::hardware_concurrency() + (primitiveSets[numProcessed+i].size()%boost::thread::hardware_concurrency() > 0) - 1;
      }

      for (size_t i = 0; i < vp_threads.size(); i++) {
	vp_threads[i]->join();
	delete vp_subCreators[i];
	delete vp_threads[i];

	_CopyAttributes(*vp_subBVHs[i]);
	_InsertSubBVH(*vp_subBVHs[i], 0, this->GetOutput().GetGeometryClusterSubtreeRootIndices()[i+numProcessed]);	
	delete vp_subBVHs[i];	
      }

      numProcessed += vp_threads.size();      
    }

#ifndef NDEBUG
    assert((int)this->GetOutput().GetLeafBVIndices().size() == this->GetMesh().GetNumberOfFacets());
    for (int f = 0; f < this->GetMesh().GetNumberOfFacets(); f++) {
      assert(this->GetOutput().GetPrimitiveBVIndex(f) >= 0 && this->GetOutput().GetPrimitiveBVIndex(f) < this->GetOutput().GetNumberOfBVs());
    }

    for (std::vector<int>::const_iterator ic_l = this->GetOutput().GetLeafBVIndices().begin(); ic_l < this->GetOutput().GetLeafBVIndices().end(); ic_l++) {
      assert(this->GetOutput().GetBV(*ic_l).PrimitiveIndex >= 0);
      assert(this->GetOutput().IsLeaf(*ic_l));
      assert(this->GetOutput().GetPrimitiveBVIndex(this->GetOutput().GetBV(*ic_l).PrimitiveIndex) == *ic_l);
    }

    for (std::vector<int>::const_iterator ic_b = this->GetOutput().GetSelfCollisionCandidates().begin(); ic_b < this->GetOutput().GetSelfCollisionCandidates().end(); ic_b++) {
      assert(this->GetOutput().GetBV(*ic_b).VolinoAngle >= BoundingVolume::GetVolinoThresholdAngle() || this->GetOutput().GetBV(*ic_b).HasDisconnectedGeometry());      
    }
#endif
  } /* if mesh big enough */  
}
