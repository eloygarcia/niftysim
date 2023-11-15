// =========================================================================
// File:       tledTrackingDeformableDeformableBVHTraverserGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    December 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledTrackingDeformableDeformableBVHTraverserGPU_CU
#define tledTrackingDeformableDeformableBVHTraverserGPU_CU

#include "tledTrackingDeformableDeformableBVHTraverserGPU.h"
#include "tledDeformableContactSurfaceGPU.h"

#include "tledTrackingBVHTraverserGPU.cu"
#include "tledDeformableDeformableBVHTraverserGPU.cu"

tledDeformableDeformableBVHTraverserGPU*  tledTrackingDeformableDeformableBVHTraverserGPU::CreateTraverser(tledSelfCollisionBVH &r_bvh) {
  tledDeformableDeformableBVHTraverserGPU *p_traverser = NULL;

  if (dynamic_cast<const tledDeformableContactSurfaceT3GPU*>(&r_bvh.GetUnspecifiedMesh()) != NULL) {
    typedef tledDeformableContactSurfaceT3GPU __Mesh;
    typedef tledNarrowConeSelfCollisionBVHImplGPU<__Mesh, tledAABB<2> > __BVH;
    typedef tledDeformableDeformableBVHTraverserImplGPU<__BVH, tledTrackingDeformableDeformableBVHTraverserGPU> __BaseTraverser;

    p_traverser = new tledTrackingDeformableDeformableBVHTraverserImplGPU<__BaseTraverser>(static_cast<__BVH&>(r_bvh));
  } else {
    tledFatalNotYetImplementedError;
  }

  return p_traverser;
}

template <class TBaseTraverser>
tledTrackingDeformableDeformableBVHTraverserImplGPU<TBaseTraverser>::~tledTrackingDeformableDeformableBVHTraverserImplGPU() {
  this->DeallocateSlaveNeighbourData();
}

template <class TBaseTraverser>
void tledTrackingDeformableDeformableBVHTraverserImplGPU<TBaseTraverser>::Init(tledUnstructuredContactManager &r_manager) {
  std::vector<int> list;
  std::vector<std::pair<int, int> > ranges;
  int *dp_list;
  int2 *dp_ranges;
  int maxRange;

  tledLogDebugStream(tledHelper::Info() << "Initialising traverser...");
  Superclass::Init(r_manager);

  this->ExtractNodeFacetNeighbourHood(list, ranges, maxRange, this->GetSlaveMesh());
  tledLogDebugStream(tledHelper::Info() << "Installing node-facet adjacency data list size = " << int(list.size()) << ", range list size = " << int(ranges.size()));
  this->ConvertHostToDeviceNeighbourLists(dp_list, dp_ranges, list, ranges);
  this->SetSlaveNodeFacetNeighbours(dp_list, dp_ranges, maxRange);
  this->SetMasterNodeFacetNeighbours(dp_list, dp_ranges, maxRange);

  this->ExtractNodeEdgeNeighbourHood(list, ranges, maxRange, this->GetSlaveMesh());
  tledLogDebugStream(tledHelper::Info() << "Installing node-edge adjacency data list size = " << int(list.size()) << ", range list size = " << int(ranges.size()));
  this->ConvertHostToDeviceNeighbourLists(dp_list, dp_ranges, list, ranges);
  this->SetSlaveNodeEdgeNeighbours(dp_list, dp_ranges, maxRange);
  this->SetMasterNodeEdgeNeighbours(dp_list, dp_ranges, maxRange);

  tledLogDebugStream(tledHelper::Info() << "Initialisation traverser complete.");
}

#endif
