// =========================================================================
// File:       tledParallelGeometricSelfCollisionBVHUpdaterCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBVH>
void tledParallelGeometricSelfCollisionBVHUpdaterCPU<TBVH>::UpdateUpdateStatusWorker(tledParallelGeometricSelfCollisionBVHUpdaterCPU *p_updater, UpdateNodeInfo *p_updateNode) {
  p_updater->UpdateSubtreeUpdateStatus(*p_updateNode);
}

template <class TBVH>
void tledParallelGeometricSelfCollisionBVHUpdaterCPU<TBVH>::UpdateBVH() {
  std::vector<boost::thread*> vp_threads;

  assert(this->GetNumberOfThreads() > 0);
  for (typename std::vector<UpdateNodeInfo>::iterator i_updateNode = this->GetUpdateNodes().begin(); i_updateNode < this->GetUpdateNodes().end();) {
    vp_threads.clear();

    while ((int)vp_threads.size() < this->GetNumberOfThreads() && i_updateNode < this->GetUpdateNodes().end()) {
      vp_threads.push_back(new boost::thread(UpdateUpdateStatusWorker, this, &(*i_updateNode)));
      i_updateNode++;
    }

    for (std::vector<boost::thread*>::iterator ip_t = vp_threads.begin(); ip_t < vp_threads.end(); ip_t++) {
      (*ip_t)->join();
      delete *ip_t;
    }
  }

  for (typename std::vector<UpdateNodeInfo>::iterator i_updateNode = this->GetUpdateNodes().begin(); i_updateNode < this->GetUpdateNodes().end(); i_updateNode++) {
    this->PerformSubtreeUpdate(*i_updateNode);
  }
}
