// =========================================================================
// File:       tledRigidMotionBVHUpdater.tpp
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

template <class TBVH>
void tledRigidMotionBVHUpdater<TBVH>::Init() {
  std::fill(m_RootTranslation, m_RootTranslation + 3, 0.0f);
}

