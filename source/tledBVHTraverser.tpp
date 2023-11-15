// =========================================================================
// File:       tledBVHTraverser.tpp
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

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImpl<TMasterBVH, TSlaveBVH, TAPI>::FindCollisions() {
  this->RunBroadPhase();
  this->RunNarrowPhase();
}
