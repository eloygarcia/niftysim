// =========================================================================
// File:       tledBasicMeshFileReader.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TMesh>
void tledBasicMeshFileReader<TMesh>::SetRotations(const float cor[], const float angleX, const float angleY, const float angleZ) {
  std::copy(cor, cor + 3, m_CentreOfRotation);
  m_RotationAngles[0] = tledPi*angleX/180.f;
  m_RotationAngles[1] = tledPi*angleY/180.f;
  m_RotationAngles[2] = tledPi*angleZ/180.f;
}

template <class TMesh>
void tledBasicMeshFileReader<TMesh>::AssembleRotationMatrix(float *p_dst) const {
  using namespace tledVectorArithmetic;

  const float a = std::cos(this->GetRotationX()), b = std::sin(this->GetRotationX());
  const float c = std::cos(this->GetRotationY()), d = std::sin(this->GetRotationY());
  const float e = std::cos(this->GetRotationZ()), f = std::sin(this->GetRotationZ());
  
  p_dst[0] = c*e, p_dst[1] = -c*f, p_dst[2] = d;
  p_dst[3] = b*d*e + a*f, p_dst[4] = -b*d*f + a*e, p_dst[5] = -b*c;
  p_dst[6] = -a*d*e + b*f, p_dst[7] = a*d*f + b*e, p_dst[8] = a*c;

  MatTranspose(p_dst, 3, 3);
  ScalarDiv(p_dst, Norm(p_dst));
  for (int c0 = 1; c0 < 3; c0++) {    
    for (int c1 = c0 - 1; c1 >= 0; c1--) {
      float proj[3];

      Sub(p_dst + 3*c0, p_dst + 3*c0, ScalarMul(proj, p_dst + 3*c1, Dot(p_dst + 3*c0, p_dst + 3*c1)));
    }
    ScalarDiv(p_dst + 3*c0, Norm(p_dst + 3*c0));
  }
  MatTranspose(p_dst, 3, 3);
}

template <class TMesh>
tledBasicMeshFileReader<TMesh>::tledBasicMeshFileReader(void) : m_ScaleFactor(1) { 
  std::fill(m_Translation, m_Translation + 3, 0.0f); 
  std::fill(m_RotationAngles, m_RotationAngles + 3, std::numeric_limits<float>::quiet_NaN());
  std::fill(m_CentreOfRotation, m_CentreOfRotation + 3, std::numeric_limits<float>::quiet_NaN());
}
