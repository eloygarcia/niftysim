// =========================================================================
// File:       tledMovingRigidContactSurface.h
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
#ifndef tledMovingRigidContactSurface_H
#define tledMovingRigidContactSurface_H

#include "tledRigidContactSurface.h"
#include "tledVectorArithmetic.h"
#include "tledXMLImporter.h"

#include <string>

/**
 * \brief Transformable rigid contact surface class.
 * \ingroup contact
 */
class tledMovingRigidContactSurface {
  /**
   * \name Transform
   * @{
   */
private:
  float m_RotationAngles[3], m_RotationCentre[3];
  float m_Translation[3];
  int m_TotalNumSteps, m_CurrentStep;  
  int m_HistoryLength;

public:
  /** Rotation applied over length of simulation. One rotation angle for each spatial axis. */
  void SetRotations(const float angles[], const float cor[]) { std::copy(angles, angles + 3, m_RotationAngles), std::copy(cor, cor + 3, m_RotationCentre); }

  /** Translation applied over length of simulation */
  void SetTranslation(const float t[]) { std::copy(t, t + 3, m_Translation); }

  /** Translation applied over length of simulation */
  const float* GetTotalTranslation(void) const { return m_Translation; }
  /** Translation at given time step */
  float* GetTranslation(float *t, const int timeStep) const { return tledVectorArithmetic::ScalarMul(t, this->GetTotalTranslation(), ((float)timeStep)/this->GetTotalNumberOfSteps()); }
  /** Translation at current time step */
  float* GetTranslation(float *t) const { return this->GetTranslation(t, this->GetCurrentStep()); }
  
  /** Centre of rotation */
  const float* GetRotationCentre(void) const { return m_RotationCentre; }
  
  /** Total angle of rotation for given axis */
  float GetTotalRotationAngle(const int axisIndex) const { return m_RotationAngles[axisIndex]; }
  const float* GetAllTotalRotationAngles(void) const { return m_RotationAngles; }
  bool HasRotation(void) const { return tledVectorArithmetic::Norm(this->GetAllTotalRotationAngles()) > 0; }
  
  /** Rotation angle at specified time step */
  float GetRotationAngle(const int axisIndex, const int timeStep) const { return (this->GetTotalRotationAngle(axisIndex)*timeStep)/this->GetTotalNumberOfSteps(); }
  /** Current rotation angle */
  float GetRotationAngle(const int axisIndex) const { return this->GetRotationAngle(axisIndex, this->GetCurrentStep()); }

  /** Total number of transform increments */
  void SetTotalNumberOfSteps(const int numSteps) { m_TotalNumSteps = numSteps; }
  int GetTotalNumberOfSteps(void) const { return m_TotalNumSteps; }

  /** Current transform increment */
  int GetCurrentStep(void) const { return m_CurrentStep; }  

  /** Apply current time step transform */
  virtual void Update(void) { m_CurrentStep += 1; }

  /** Number of time steps to look back for "old" node positions, normals, etc. (default: 10) */
  void SetHistoryLength(const int hLen) { m_HistoryLength = hLen; }
  
  int GetHistoryLength(void) const { return m_HistoryLength; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledMovingRigidContactSurface(void);

public:
  tledRigidContactSurface* CreateSurface(const std::string &type, const bool useGPU);

  virtual ~tledMovingRigidContactSurface(void) {}
  /** @} */
};

template <const int t_numFacetVertices, class TAPI = tledMovingRigidContactSurface> 
class tledMovingRigidContactSurfaceImpl : public tledRigidContactSurfaceImpl<t_numFacetVertices, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledRigidContactSurfaceImpl<t_numFacetVertices, TAPI> Superclass;
  /** @} */

  /**
   * \name Surface Type Queries
   * @{
   */
public:
  virtual bool IsMoving(void) const { return true; }
  /** @} */

  /**
   * \name Output Untemplated Interface
   * @{
   */
public:
  float* GetFinalDisplacement(float *p_dst) const;
  /** @} */

  /**
   * \name Construction, Destruction
   */
public:
  class MovingSurfaceXMLImporter;

public:
  /** Should be instantiated through tledRigidContactSurface::CreateSurface */
  tledMovingRigidContactSurfaceImpl(void) {}

  virtual ~tledMovingRigidContactSurfaceImpl(void) {}
  /** @} */
};

template <const int t_numFacetVertices, class TAPI> 
class tledMovingRigidContactSurfaceImpl<t_numFacetVertices, TAPI>::MovingSurfaceXMLImporter : public Superclass::BasicSurfaceXMLImporter {
protected:
  virtual void ReadBasicMeshSpec(void);
  virtual void ReadMotion(void);

public:
  virtual void Import(void);

public:
  virtual ~MovingSurfaceXMLImporter(void) {}
};

#include "tledMovingRigidContactSurface.tpp"
#endif
