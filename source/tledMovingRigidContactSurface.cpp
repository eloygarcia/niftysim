// =========================================================================
// File:       tledMovingRigidContactSurface.cpp
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
#include "tledMovingRigidContactSurface.h"
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledHelper.h"
#ifdef GPU_GP_CONTACT
#include "tledMovingRigidContactSurfaceGPU.h"
#endif

#include <algorithm>

tledMovingRigidContactSurface::tledMovingRigidContactSurface(void) : m_TotalNumSteps(0), m_CurrentStep(0), m_HistoryLength(10) {
  std::fill(m_RotationAngles, m_RotationAngles + 3, std::numeric_limits<float>::quiet_NaN());
  std::fill(m_Translation, m_Translation + 3, 0.0f);
}

tledRigidContactSurface* tledMovingRigidContactSurface::CreateSurface(const std::string &type, const bool useGPU) {
  if (useGPU) {
#ifdef GPU_GP_CONTACT
    return tledMovingRigidContactSurfaceGPU::CreateSurface(type);
#else
    tledFatalFeatureNotEnabledError;
    return NULL;
#endif
  } else {
    return tledMovingRigidContactSurfaceCPU::CreateSurface(type);
  }
}
