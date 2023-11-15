// =========================================================================
// File:       tledOBB.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::ComputeFromNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]) {
  using namespace tledVectorArithmetic;

  std::fill(this->Centroid, this->Centroid + 3, 0.f);
  for (int const *pc_n = nodeListStart; pc_n < nodeListEnd; pc_n++) Add(this->Centroid, this->Centroid, nodes + (*pc_n)*3);
  ScalarDiv(this->Centroid, float(nodeListEnd - nodeListStart));

  if (nodeListEnd - nodeListStart > 3) {
    {
      float C[3][3], EVals[3];

      std::fill(&C[0][0], &C[0][0] + 3*3, 0.f);
      for (int const *pc_n = nodeListStart; pc_n < nodeListEnd; pc_n++) {
	float d[3];

	Sub(d, nodes + (*pc_n)*3, this->Centroid);
	for (int r = 0; r < 3; r++) for (int c = r; c < 3; c++) C[r][c] += d[r]*d[c];
      }
      for (int r = 0; r < 3; r++) for (int c = 0; c < r; c++) C[r][c] = C[c][r];
  
      SymmetricEVD3x3(&this->Axes[0][0], EVals, &C[0][0]);
      for (int r = 0; r < 3; r++) for (int c = r; c < 3; c++) std::iter_swap(&this->Axes[r][c], &this->Axes[c][r]);
    }

    {
      float maxs[3], mins[3], tmp[3];

      for (int c = 0; c < 3; c++) maxs[c] = -(mins[c] = std::numeric_limits<float>::max());
      for (int const *pc_n = nodeListStart; pc_n < nodeListEnd; pc_n++) {
	for (int c = 0; c < 3; c++) {
	  const float proj = Dot(this->Axes[c], nodes + 3*(*pc_n));

	  maxs[c] = std::max(maxs[c], proj);
	  mins[c] = std::min(mins[c], proj);
	}
      }

      std::fill(this->Centroid, this->Centroid + 3, 0.f);
      for (int c = 0; c < 3; c++) Add(this->Centroid, this->Centroid, ScalarMul(tmp, this->Axes[c], (maxs[c] + mins[c])/2));
    }
  } else if (nodeListEnd - nodeListStart == 3) {
    float maxEdgeMag, normA1;
    float const *pc_v0;

    pc_v0 = nodes + 3*(*nodeListStart);
    maxEdgeMag = Norm(Sub(this->Axes[0], pc_v0, nodes + 3*(*(nodeListEnd - 1))));
    ScalarDiv(this->Axes[0], maxEdgeMag);
    for (int const *pc_n = nodeListStart; pc_n < nodeListEnd - 1; pc_n++) {
      float edge[3], mag;
      
      if ((mag = Norm(Sub(edge, nodes + 3*(*pc_n), nodes + 3*(*(pc_n + 1))))) > maxEdgeMag) {
	maxEdgeMag = mag;
	pc_v0 = nodes + 3*(*pc_n);
	ScalarDiv(this->Axes[0], edge, mag);
      }
    }

    OrthogonaliseNormalised(this->Axes[1], this->Axes[0], Sub(this->Axes[1], pc_v0, this->Centroid), false);
    if ((normA1 = Norm(this->Axes[1])) < 1e-5f) {
      int minAxis = 0;
      float minAxisAbsV = std::fabs(this->Axes[1][0]);

      /* points lie on a line */
      for (int a = 1; a < 3; a++) {
	if (std::fabs(this->Axes[1][a]) < minAxisAbsV) {
	  minAxisAbsV = std::fabs(this->Axes[1][a]);
	  minAxis = a;
	}
      }

      std::fill(this->Axes[1], this->Axes[1] + 3, 0.f);
      this->Axes[1][minAxis] = 1.f;
      OrthogonaliseNormalised(this->Axes[1], this->Axes[0], this->Axes[1], true);      
    } else {
      ScalarDiv(this->Axes[1], normA1);
    }
    Cross(this->Axes[2], this->Axes[0], this->Axes[1]);
  } else {
    tledFatalError("Need at least 3 points for fit.");
  }

  _ResetExtent();
  for (int const *pc_n = nodeListStart; pc_n < nodeListEnd; pc_n++) {
    _UpdateExtent(nodes + 3*(*pc_n));
  }

  _SortAxes();
  _CorrectAxisInversion();

#if !defined NDEBUG && defined __OBB_QUATERNION_MERGE
  for (int r = 0; r < 3; r++) {
    assert(this->Axes[r][r] + 1 >= 1e-4f);
  }
#endif
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::_SwapAxes(const int a0Ind, const int a1Ind) {
  for (int c = 0; c < 3; c++) {
    std::iter_swap(this->Axes[a0Ind] + c, this->Axes[a1Ind] + c);
  }
  std::iter_swap(this->Extents + a0Ind, this->Extents + a1Ind);
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::_SortAxes() {
  assert(this->Extents[0] >= 0 && this->Extents[1] >= 0 && this->Extents[2] >= 0);

  if (this->Extents[1] > this->Extents[0] && this->Extents[1] > this->Extents[2]) _SwapAxes(0, 1);
  else if (this->Extents[2] > this->Extents[0]) _SwapAxes(0, 2);
  assert(this->Extents[0] >= this->Extents[1] && this->Extents[0] >= this->Extents[2]);

  if (this->Extents[1] < this->Extents[2]) _SwapAxes(1, 2);
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::_CorrectAxisInversion() {
#ifdef __OBB_QUATERNION_MERGE
  for (int c = 0; c < 3; c++) {
    float maxAbsComp = 0.f, maxAbsVal = 0.f;

    for (int r = 0; r < 3; r++) if (std::fabs(this->Axes[c][r]) > maxAbsVal) {
	maxAbsVal = std::fabs(maxAbsComp = this->Axes[c][r]);
      }
    if (maxAbsComp < 0.f) tledVectorArithmetic::ScalarMul(this->Axes[c], -1);
  }
#endif
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::_ResetExtent() {
  std::fill(this->Extents, this->Extents + 3, 0.f);
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::_UpdateExtent(const float x[]) {
  using namespace tledVectorArithmetic;

  float lx[3];

  Sub(lx, x, this->Centroid);
  for (int c = 0; c < 3; c++) {
    this->Extents[c] = std::max(this->Extents[c], std::fabs(tledVectorArithmetic::Dot(lx, this->Axes[c])));
  }
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::ExpandWithNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]) {
  for (int const *pc_n = nodeListStart; pc_n < nodeListEnd; pc_n++) {
    _UpdateExtent(nodes + *pc_n*3);
  }
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::_ComputeVertices(float (*p_vDst)[3]) const {
  using namespace tledVectorArithmetic;

  float halfEdges[3][3];

  for (int a = 0; a < 3; a++) {
    ScalarMul(halfEdges[a], this->Axes[a], this->Extents[a]);
  }
  
  for (int x = -1, v = 0; x <= 1; x += 2) {
    float bx[3], by[3];

    Add(bx, this->Centroid, ScalarMul(bx, halfEdges[0], float(x)));
    for (int y = -1; y <= 1; y += 2) {
      ScalarMul(by, halfEdges[1], float(y));
      Add(by, by, bx);
      for (int z = -1; z <= 1; z += 2, v++) {
	assert(v < 8);
	Add(p_vDst[v], ScalarMul(p_vDst[v], halfEdges[2], float(z)), by);
      }
    }
  }
}

#if defined __OBB_QUATERNION_MERGE
template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::Merge(tledOBB &r_bv0, const tledOBB &bv1) {
  if (r_bv0.ComputeVolume() < bv1.ComputeVolume()) {
    tledOBB tmpOBB;

    tledOBB::CopyBoundsFromBV(tmpOBB, bv1);
    tledOBB::Merge(tmpOBB, r_bv0);
    tledOBB::CopyBoundsFromBV(r_bv0, tmpOBB);
  } else {
    using namespace tledVectorArithmetic;

    float vertices[2][8][3];
    bool isContained = true;

    bv1._ComputeVertices(vertices[1]);    
    for (int v = 0; v < 8 && isContained; v++) isContained &= r_bv0.IsInside(vertices[1][v]);
    if (!isContained) {
      r_bv0._ComputeVertices(vertices[0]);

      {
	float q0[4], q1[4], qm[4], qMag;

	ComputeQuaternionFromRotation(q0, &r_bv0.Axes[0][0]);
	ComputeQuaternionFromRotation(q1, &bv1.Axes[0][0]);
	MatAdd(q0, q1, 1, 4, qm);
	if ((qMag = Norm(qm, 4)) < 1e-5f) {
	  if (r_bv0.ComputeVolume() < bv1.ComputeVolume()) std::copy(&bv1.Axes[0][0], &bv1.Axes[0][0] + 3*3, &r_bv0.Axes[0][0]);
	} else {
	  MatMultScalar(qm, 1, 4, 1.f/qMag, qm);    
	  ComputeRotationFromQuaternion(&r_bv0.Axes[0][0], qm);
	}
      }

      std::fill(r_bv0.Centroid, r_bv0.Centroid + 3, 0.f);
      for (int b = 0; b < 2; b++) {
	for (float const (*ppc_v)[3] = vertices[b]; ppc_v < vertices[b] + 8; ppc_v++) Add(r_bv0.Centroid, r_bv0.Centroid, *ppc_v);
      }
      ScalarDiv(r_bv0.Centroid, 2.f*8.f);

      r_bv0._ResetExtent();
      for (int b = 0; b < 2; b++) {
	for (float const (*ppc_v)[3] = vertices[b]; ppc_v < vertices[b] + 8; ppc_v++) r_bv0._UpdateExtent(*ppc_v);
      }
      r_bv0._SortAxes();
      r_bv0._CorrectAxisInversion();
    } /* if is contained */
  } /* if r_bv0 is smaller else .. */

#ifndef NDEBUG
  for (int r = 0; r < 3; r++) {
    assert(r_bv0.Axes[r][r] + 1 >= 1e-4f);
  }
#endif
}
#elif defined __OBB_AXIS_MERGE
template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::Merge(tledOBB &r_bv0, const tledOBB &bv1) {
  if (r_bv0.ComputeVolume() < bv1.ComputeVolume()) {
    tledOBB tmpOBB;

    tledOBB::CopyBoundsFromBV(tmpOBB, bv1);
    tledOBB::Merge(tmpOBB, r_bv0);
    tledOBB::CopyBoundsFromBV(r_bv0, tmpOBB);
  } else {
    using namespace tledVectorArithmetic;

    float vertices[2][8][3];
    bool isContained = true;

    bv1._ComputeVertices(vertices[1]);    
    for (int v = 0; v < 8 && isContained; v++) isContained &= r_bv0.IsInside(vertices[1][v]);
    if (!isContained) {
      bool isUsed[3];

      r_bv0._ComputeVertices(vertices[0]);
      std::fill(isUsed, isUsed + 3, false);
      for (int c0 = 0; c0 < 2; c0++) {
	int maxAxis = -1, aSign;
	float maxOlap = -1.f, olap, tmpAxis[3];

	for (int c1 = 0; c1 < 3; c1++) if (!isUsed[c1]) {
	    olap = Dot(r_bv0.Axes[c0], bv1.Axes[c1]);	    
	    if (std::fabs(olap) >= maxOlap) {
	      aSign = olap < 0? -1 : 1;
	      maxOlap = aSign*olap;
	      maxAxis = c1;
	    }
	  }
	assert(maxAxis >= 0 && maxAxis < 3);
	assert(aSign == -1 || aSign == 1);

	isUsed[maxAxis] = true;
	Add(r_bv0.Axes[c0], r_bv0.Axes[c0], ScalarMul(tmpAxis, bv1.Axes[maxAxis], aSign));
	ScalarDiv(r_bv0.Axes[c0], Norm(r_bv0.Axes[c0]));
      }      
      OrthogonaliseNormalised(r_bv0.Axes[1], r_bv0.Axes[0], r_bv0.Axes[1], true);
      Cross(r_bv0.Axes[2], r_bv0.Axes[0], r_bv0.Axes[1]);

      std::fill(r_bv0.Centroid, r_bv0.Centroid + 3, 0.f);
      for (int b = 0; b < 2; b++) {
	for (float const (*ppc_v)[3] = vertices[b]; ppc_v < vertices[b] + 8; ppc_v++) Add(r_bv0.Centroid, r_bv0.Centroid, *ppc_v);
      }
      ScalarDiv(r_bv0.Centroid, 2.f*8.f);

      r_bv0._ResetExtent();
      for (int b = 0; b < 2; b++) {
	for (float const (*ppc_v)[3] = vertices[b]; ppc_v < vertices[b] + 8; ppc_v++) r_bv0._UpdateExtent(*ppc_v);
      }
      r_bv0._SortAxes();
    } /* if is contained */
  } /* if r_bv0 is smaller else .. */
}
#else
template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::Merge(tledOBB &r_bv0, const tledOBB &bv1) {
  const float vol0 = r_bv0.ComputeVolume(), vol1 = bv1.ComputeVolume();

  if (vol0 < vol1) {
    tledOBB tmpOBB;

    tledOBB::CopyBoundsFromBV(tmpOBB, bv1);
    tledOBB::Merge(tmpOBB, r_bv0);
    tledOBB::CopyBoundsFromBV(r_bv0, tmpOBB);
  } else {
    using namespace tledVectorArithmetic;

    float vertices[2][8][3];
    bool isContained = true;

    bv1._ComputeVertices(vertices[1]);    
    for (int v = 0; v < 8 && isContained; v++) isContained &= r_bv0.IsInside(vertices[1][v]);
    if (!isContained) {
      bool isUsed[3];

      r_bv0._ComputeVertices(vertices[0]);
      std::fill(isUsed, isUsed + 3, false);
      for (int c0 = 0; c0 < 2; c0++) {
	int maxAxis = -1, aSign;
	float maxOlap = -1.f, olap, axis0[3], axis1[3];

	for (int c1 = 0; c1 < 3; c1++) if (!isUsed[c1]) {
	    olap = Dot(r_bv0.Axes[c0], bv1.Axes[c1]);	    
	    if (std::fabs(olap) >= maxOlap) {
	      aSign = olap < 0? -1 : 1;
	      maxOlap = aSign*olap;
	      maxAxis = c1;
	    }
	  }
	assert(maxAxis >= 0 && maxAxis < 3);
	assert(aSign == -1 || aSign == 1);

	isUsed[maxAxis] = true;
	Add(r_bv0.Axes[c0], ScalarMul(axis0, r_bv0.Axes[c0], vol0), ScalarMul(axis1, bv1.Axes[maxAxis], aSign*vol1));
	ScalarDiv(r_bv0.Axes[c0], Norm(r_bv0.Axes[c0]));
      }      
      OrthogonaliseNormalised(r_bv0.Axes[1], r_bv0.Axes[0], r_bv0.Axes[1], true);
      Cross(r_bv0.Axes[2], r_bv0.Axes[0], r_bv0.Axes[1]);

      std::fill(r_bv0.Centroid, r_bv0.Centroid + 3, 0.f);
      for (int b = 0; b < 2; b++) {
	for (float const (*ppc_v)[3] = vertices[b]; ppc_v < vertices[b] + 8; ppc_v++) Add(r_bv0.Centroid, r_bv0.Centroid, *ppc_v);
      }
      ScalarDiv(r_bv0.Centroid, 2.f*8.f);

      r_bv0._ResetExtent();
      for (int b = 0; b < 2; b++) {
	for (float const (*ppc_v)[3] = vertices[b]; ppc_v < vertices[b] + 8; ppc_v++) r_bv0._UpdateExtent(*ppc_v);
      }
      r_bv0._SortAxes();
    } /* if is contained */
  } /* if r_bv0 is smaller else .. */
}
#endif

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::AddMargin(const float margin) {
  for (int c = 0; c < 3; c++) this->Extents[c] += margin;
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::Translate(const float t[]) {
  tledVectorArithmetic::Add(this->Centroid, this->Centroid, t);
}

template <const int t_numChildBVs>
bool tledOBB<t_numChildBVs>::IsInside(const float pt[]) const {
  using namespace tledVectorArithmetic;

  float dCent[3];
  int c;

  /* Return false on NaN bounds, do not early exit on proj > extent! */
  Sub(dCent, pt, this->Centroid);
  for (c = 0; c < 3 && std::fabs(Dot(this->Axes[c], dCent)) <= this->Extents[c]; c++);

  return c == 3;
}

template <const int t_numChildBVs>
bool tledOBB<t_numChildBVs>::DoesIntersect(const float a[], const float b[]) const {
  tledFatalNotYetImplementedError;

  return false;
}

template <const int t_numChildBVs>
const float* tledOBB<t_numChildBVs>::_ComputeSATAxis(float* computeBuffer, bool &r_isSameAxis, const int axisIndex, const tledOBB &obb) const {
  r_isSameAxis = false;

  if (axisIndex < 3) return this->Axes[axisIndex];
  else if (axisIndex < 6) return obb.Axes[axisIndex-3];
  else {
    using namespace tledVectorArithmetic;

    const int crossAxisIndex = axisIndex - 6;
    
    float nAxis;

    assert(crossAxisIndex >= 0 && crossAxisIndex < 9);
    nAxis = Norm(Cross(computeBuffer, this->Axes[crossAxisIndex/3], obb.Axes[crossAxisIndex%3]));

    if (nAxis > 1e-6f) return ScalarDiv(computeBuffer, nAxis);
    else {
      r_isSameAxis = true;

      return computeBuffer;
    }
  }
}

template <const int t_numChildBVs>
float tledOBB<t_numChildBVs>::_ComputeBoundsOnAxis(const float axis[]) const {
  float r = 0.f;

  for (int c = 0; c < 3; c++) {
    r = std::max(r, std::fabs(tledVectorArithmetic::Dot(axis, this->Axes[c]))*this->Extents[c]);
  }

  return r;
}

template <const int t_numChildBVs>
bool tledOBB<t_numChildBVs>::_DoOverlapOnAxis(const float axis[], const tledOBB &obb) const {
  using namespace tledVectorArithmetic;

  float dCent[3];

  Sub(dCent, this->Centroid, obb.Centroid);

  return std::fabs(Dot(dCent, axis)) <= _ComputeBoundsOnAxis(axis) + obb._ComputeBoundsOnAxis(axis);
}

template <const int t_numChildBVs>
bool tledOBB<t_numChildBVs>::DoesIntersect(const tledBV<t_numChildBVs> &bv) const {
  const tledOBB &obb = static_cast<const tledOBB&>(bv);

  assert(dynamic_cast<const tledOBB*>(&bv) != NULL);
  for (int a = 0; a < 15; a++) {
    float const *pc_axis;
    float axisBuffer[3];
    bool isIdentAxis;

    pc_axis = _ComputeSATAxis(axisBuffer, isIdentAxis, a, obb);
  
    if (!isIdentAxis && !_DoOverlapOnAxis(pc_axis, obb)) return false;
  }

  return true;
}

template <const int t_numChildBVs>
void tledOBB<t_numChildBVs>::CopyBoundsFromBV(tledOBB &r_dstBV, const tledOBB &srcBV) {
  std::copy(srcBV.Centroid, srcBV.Centroid + 3, r_dstBV.Centroid);
  std::copy(srcBV.Extents, srcBV.Extents + 3, r_dstBV.Extents);
  std::copy(&srcBV.Axes[0][0], &srcBV.Axes[0][0] + 3*3, &r_dstBV.Axes[0][0]);
}

#if defined _GPU_ && defined __CUDACC__
template <>
void tledOBB<2>::InitGPU(Superclass::GPUBV &r_dst);

template <>
void tledOBB<4>::InitGPU(Superclass::GPUBV &r_dst);
#endif
