// =========================================================================
// File:       tledBasicMeshFileReader.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBasicMeshFileReader_H
#define tledBasicMeshFileReader_H

#include "tledBasicMeshCreator.h"
#include "tledMatrixFunctions.h"
#include "tledVectorArithmetic.h"

#include <algorithm>
#include <string>
#include <limits>

/**
 * \defgroup fileio File I/O
 */

/**
 * \defgroup import Mesh-File Readers
 * \ingroup fileio
 */

/**
 * \brief Base class of all file readers
 * \ingroup import
 *
 * File readers are not intended to be re-usable, a new instance ought to be created for every read.
 */
template <class TMesh>
class tledBasicMeshFileReader : public tledBasicMeshCreator<TMesh> {
  /**
   * \name Transforms
   * @{
   */
private:
  float m_ScaleFactor, m_Translation[3], m_CentreOfRotation[3], m_RotationAngles[3];

protected:
  virtual void ApplyTransforms(void) = 0;
  void AssembleRotationMatrix(float *p_dst) const;

public:
  /** Sets the translation vector, this is the last transform to be applied */
  void SetTranslation(const float trans[]) { std::copy(trans, trans + 3, m_Translation); }
  const float* GetTranslation(void) const { return m_Translation; }

  /** Mesh scaling: expansion/contraction about origin */
  void SetScaleFactor(const float s) { m_ScaleFactor = s; }
  float GetScaleFactor(void) const { return m_ScaleFactor; }
  
  /** @{ */
  /** 
   * Can be used to specify a rotation about \"cor\" (centre of rotation). Angles are defined with respect to coordinate system axes: angleX = angle of rotation about x-axis, etc.
   * The order of application is \f$x' = R(\alpha_{X})R(\alpha_{Y})R(\alpha_{Z})x\f$. Rotation, if requested, is the first transform to be applied to the mesh (before scaling, translation).
   * NaNs are used to indicate that no rotation was requested.
   */
  void SetRotations(const float cor[], const float angleX, const float angleY, const float angleZ);
  const float* GetCentreOfRotation(void) const { return m_CentreOfRotation; }
  float GetRotationX(void) const { return m_RotationAngles[0]; }
  float GetRotationY(void) const { return m_RotationAngles[1]; }
  float GetRotationZ(void) const { return m_RotationAngles[2]; }
  /** @} */
  /** @} */

  /**
   * \name FS I/O
   * @{
   */
private:
  std::string m_Filename;

public:
  void SetFilename(const std::string &fn) { m_Filename = fn; }
  const std::string& GetFilename(void) const { return m_Filename; }

  virtual void Read(void) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledBasicMeshFileReader(void);
  virtual ~tledBasicMeshFileReader(void) {}
  /** @} */
};

#include "tledBasicMeshFileReader.tpp"
#endif
