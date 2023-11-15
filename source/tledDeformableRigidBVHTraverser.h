// =========================================================================
// File:       tledDeformableRigidBVHTraverser.h
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
#ifndef tledDeformableRigidBVHTraverser_H
#define tledDeformableRigidBVHTraverser_H

#include "tledBVHTraverser.h"

class tledDeformableRigidBVHTraverser : public tledBVHTraverser {
  /**
   * \name Slave Node Masks and Lists
   * @{
   */
private:
  

public:
  bool IsSlaveNode(const int nodeIndex) const { 
  /** @} */
};

#endif
