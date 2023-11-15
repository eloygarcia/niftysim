// =========================================================================
// File:       tledContactSurfaceExporter.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactSurfaceExporter_H
#define tledContactSurfaceExporter_H

#include "tledContactSurface.h"
#include "tledMatrixFunctions.h"
#include "tledDeformableContactSurface.h"

#include <algorithm>
#include <vector>
#include <sstream>

class tledContactSurfaceExporter {
  /**
   * \name Contact Forces
   * @{
   */
private:
  std::vector<float> m_ContactForceAccumulator;
  int m_NumContactForceSamples;

public:
  /** Resets the contact force accumulator */
  void ResetContactForces(void);

  /** Adds a new sample to the contact force accumulator */
  void AddContactForceSample(const std::vector<float> &contactForces);
  /** @} */

  /**
   * \name Surface
   * @{
   */
private:
  const tledContactSurface *mpc_Surface;

public:
  void SetSurface(const tledContactSurface &surface) { mpc_Surface = &surface; }  
  /** @} */

  /**
   * \name I/O
   * @{
   */
private:
  std::string m_FilePrefix;
  
public:
  /** Prefix for output files */
  void SetFilePrefix(const std::string &pref) { m_FilePrefix = pref; }  

  /** Writes the current surface and contact forces (average of forces since last reset) to a file located at <FILE PREFIX>_t=<TIME>.vtk */
  void Write(const float t);
  /** @} */
};

#endif
