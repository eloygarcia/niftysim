template <const int t_NumChildBVs>
void tledAABB<t_NumChildBVs>::ComputeFromNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]) {
  for (int cInd = 0; cInd < 3; cInd++) this->Bounds[cInd][1] = this->Bounds[cInd][0] = nodes[3*nodeListStart[0]+cInd];
  this->ExpandWithNodeList(nodeListStart + 1, nodeListEnd, nodes);
}

template <const int t_NumChildBVs>
void tledAABB<t_NumChildBVs>::ExpandWithNodeList(const int *nodeListStart, const int *nodeListEnd, const float nodes[]) {
  for (int const *pc_nodeInd = nodeListStart; pc_nodeInd < nodeListEnd; pc_nodeInd++) {
    const float *x = nodes + 3*(*pc_nodeInd);

    for (int cInd = 0; cInd < 3; cInd++) {
      this->Bounds[cInd][0] = std::min(this->Bounds[cInd][0], x[cInd]);
      this->Bounds[cInd][1] = std::max(this->Bounds[cInd][1], x[cInd]);
    }
  }
}

template <const int t_NumChildBVs>
void tledAABB<t_NumChildBVs>::Translate(const float t[]) {
  for (int cInd = 0; cInd < 3; cInd++) this->Bounds[cInd][0] += t[cInd], this->Bounds[cInd][1] += t[cInd];
}

template <const int t_numChildBVs>
void tledAABB<t_numChildBVs>::CopyBoundsFromBV(tledAABB &r_dstBV, const tledAABB &src)  { 
  std::copy(&src.Bounds[0][0], &src.Bounds[0][0] + 2*3, &r_dstBV.Bounds[0][0]);
}

template <const int t_NumChildBVs>
void tledAABB<t_NumChildBVs>::Merge(tledAABB &r_bv0, const tledAABB &bv1) {
  for (int compInd = 0; compInd < 3; compInd++) {
    r_bv0.Bounds[compInd][0] = std::min(r_bv0.Bounds[compInd][0], bv1.Bounds[compInd][0]);
    r_bv0.Bounds[compInd][1] = std::max(r_bv0.Bounds[compInd][1], bv1.Bounds[compInd][1]);
  }
}

template <const int t_numChildBVs>
bool tledAABB<t_numChildBVs>::DoesIntersect(const tledBV<t_numChildBVs> &bv) const {
  const tledAABB &aabb1 = static_cast<const tledAABB&>(bv);

  int c;

  /* Check for overlap in all dimensions. Use positive statement to avoid triggering contact handling when simulation is diverging, bounds are NaN. */
  for (c = 0; c < 3 && (this->Bounds[c][1] >= aabb1.Bounds[c][0] && aabb1.Bounds[c][1] >= this->Bounds[c][0]); c++);

  return c == 3;
}

template <const int t_NumChildBVs>
bool tledAABB<t_NumChildBVs>::DoesIntersect(const float edge0[], const float edge1[]) const {
  using namespace tledVectorArithmetic;

  /*
   * Does intersect if:
   * - either endpoint is inside the box 
   * - ray intersects with any face
   */
  for (int s = 0; s < 3; s++) {
    if ((edge1[s] < this->Bounds[s][0] && edge0[s] < this->Bounds[s][0]) || (edge1[s] > this->Bounds[s][1] && edge0[s] > this->Bounds[s][1])) return false;
  }

  for (int s = 0; s < 2; s++) {
    float xcoll[3], tedge0[3], tedge1[3];
    float t;

    t = (this->Bounds[0][s] - edge0[0])/(edge1[0] - edge0[0]);
    if (t >= 0 && t <= 1) {
      Add(xcoll, ScalarMul(tedge0, edge0, (1 - t)), ScalarMul(tedge1, edge1, t));
      if (xcoll[1] > this->Bounds[1][0] && xcoll[1] < this->Bounds[1][1] && xcoll[2] > this->Bounds[2][0] && xcoll[2] < this->Bounds[2][1]) return true;
    }

    t = (this->Bounds[1][s] - edge0[1])/(edge1[1] - edge0[1]);
    if (t >= 0 && t <= 1) {
      Add(xcoll, ScalarMul(tedge0, edge0, (1 - t)), ScalarMul(tedge1, edge1, t));
      if (xcoll[0] > this->Bounds[0][0] && xcoll[0] < this->Bounds[0][1] && xcoll[2] > this->Bounds[2][0] && xcoll[2] < this->Bounds[2][1]) return true;
    }

    t = (this->Bounds[2][s] - edge0[2])/(edge1[2] - edge0[2]);
    if (t >= 0 && t <= 1) {
      Add(xcoll, ScalarMul(tedge0, edge0, (1 - t)), ScalarMul(tedge1, edge1, t));
      if (xcoll[0] > this->Bounds[0][0] && xcoll[0] < this->Bounds[0][1] && xcoll[1] > this->Bounds[1][0] && xcoll[1] < this->Bounds[1][1]) return true;
    }
  }

  return (this->IsInside(edge0) || this->IsInside(edge1));    
} /* DoesIntersect */

template <const int t_NumChildBVs>
void tledAABB<t_NumChildBVs>::AddMargin(const float margin) {
  assert(margin >= 0);
  for (int cInd = 0; cInd < 3; cInd++) {
    this->Bounds[cInd][0] -= margin;
    this->Bounds[cInd][1] += margin;
  }
}

#ifdef _GPU_
template <>
void tledAABB<2>::InitGPU(Superclass::GPUBV &r_dst);

template <>
void tledAABB<4>::InitGPU(Superclass::GPUBV &r_dst);
#endif
