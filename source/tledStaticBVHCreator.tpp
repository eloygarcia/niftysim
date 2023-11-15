// =========================================================================
// File:       tledStaticBVHCreator.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2015
// 
// Copyright (c) 2015, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

/*
 * Finds the axis on which the centroids have the widest spread
 */
template <class TBVH>
int tledStaticBVHCreator<TBVH>::FindCentroidMaxAxis(const std::vector<int>::const_iterator ic_pIndsBegin, const std::vector<int>::const_iterator ic_pIndsEnd) const {
  float cbnds[3][2];

  for (int cInd = 0; cInd < 3; cInd++) {
    cbnds[cInd][0] = cbnds[cInd][1] = this->GetPrimitiveCentroids()[3*(*ic_pIndsBegin)+cInd];
  }

  for (std::vector<int>::const_iterator ic_pidx = ic_pIndsBegin; ic_pidx < ic_pIndsEnd; ic_pidx++) {
    const float *pcentroid = &this->GetPrimitiveCentroids()[3*(*ic_pidx)];

    for (int cInd = 0; cInd < 3; cInd++) {
      if (cbnds[cInd][0] > pcentroid[cInd]) cbnds[cInd][0] = pcentroid[cInd];
      else if (cbnds[cInd][1] < pcentroid[cInd]) cbnds[cInd][1] = pcentroid[cInd];
    }
  }

  {
    float max_dim;
    int max_cInd;

    max_dim = cbnds[0][1] - cbnds[0][0], max_cInd = 0;
    assert(max_dim >= 0);

    for (int cInd = 1; cInd < 3; cInd++) if (max_dim < cbnds[cInd][1] - cbnds[cInd][0]) {
      max_dim = cbnds[cInd][1] - cbnds[cInd][0], max_cInd = cInd;
      assert(cbnds[cInd][1] > cbnds[cInd][0]);
    }

    return max_cInd;
  }  
}

/*
 * Finds the average of the centroid values on one axis
 */
template <class TBVH>
float tledStaticBVHCreator<TBVH>::FindCentroidAxisAvg(const int cInd, const std::vector<int>::const_iterator ic_pIndsBegin, const std::vector<int>::const_iterator ic_pIndsEnd) const {
  std::vector<int>::const_iterator ic_pidx;
  float avg;

  avg = 0;
  for (ic_pidx = ic_pIndsBegin; ic_pidx < ic_pIndsEnd; ic_pidx++) avg += this->GetPrimitiveCentroids()[3*(*ic_pidx)+cInd];

  return avg/(ic_pIndsEnd - ic_pIndsBegin);
} 

template <class TBVH>
void tledStaticBVHCreator<TBVH>::SplitBV(std::vector<int>::iterator ppi_childIndexBounds[][2], const int bvIndex, const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd) {
  if (TBVH::BVHOrder == 2) {
    _SplitBinary(ppi_childIndexBounds, pIndsBegin, pIndsEnd);
  } else {
    std::vector<int>::iterator tmpChildIndexBoundsLvl1[2][2];

    _SplitBinary(tmpChildIndexBoundsLvl1, pIndsBegin, pIndsEnd);
    if (TBVH::BVHOrder == 4) {
      int cCtr = 0;
      
      for (int tc = 0; tc < 2; tc++) {
	if (tmpChildIndexBoundsLvl1[tc][1] - tmpChildIndexBoundsLvl1[tc][0] > 1) {
	  _SplitBinary(ppi_childIndexBounds + cCtr, tmpChildIndexBoundsLvl1[tc][0], tmpChildIndexBoundsLvl1[tc][1]);
	  cCtr += 2;
	} else {
	  ppi_childIndexBounds[cCtr][0] = tmpChildIndexBoundsLvl1[tc][0];
	  ppi_childIndexBounds[cCtr][1] = tmpChildIndexBoundsLvl1[tc][1];
	  cCtr += 1;
	}
      }
      for (; cCtr < TBVH::BVHOrder; cCtr++) ppi_childIndexBounds[cCtr][0] = ppi_childIndexBounds[cCtr][1] = pIndsEnd;
    } else {
      tledFatalNotYetImplementedError;
    }
  }
} /* SplitBV */

template <class TBVH>
void tledStaticBVHCreator<TBVH>::_SplitBinary(std::vector<int>::iterator p_childIndexBounds[][2], const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd) {
  const int div_cInd = this->FindCentroidMaxAxis(pIndsBegin, pIndsEnd);
  const float div_val = this->FindCentroidAxisAvg(div_cInd, pIndsBegin, pIndsEnd);

  std::vector<int>::iterator i_fst, i_startSnd;

  /*
   * Traverse pInds array from two sides
   */
  i_fst = pIndsBegin, i_startSnd = pIndsEnd - 1;
  while (i_fst < i_startSnd) {
    if (this->GetPrimitiveCentroids()[(*i_fst)*3+div_cInd] < div_val) i_fst += 1;
    else {
      std::iter_swap(i_startSnd, i_fst);

      assert(this->GetPrimitiveCentroids()[(*i_startSnd)*3+div_cInd] >= div_val);
      i_startSnd -= 1;
    }
    assert(i_startSnd < pIndsEnd);
    assert(i_fst >= pIndsBegin);

    assert(i_fst == pIndsBegin || this->GetPrimitiveCentroids()[(*(i_fst-1))*3+div_cInd] < div_val);
  }

  if (this->GetPrimitiveCentroids()[(*i_fst)*3+div_cInd] < div_val) i_fst += 1;

  assert(i_fst - pIndsBegin > 0 && pIndsEnd - i_fst > 0);
  assert(i_fst < pIndsEnd && this->GetPrimitiveCentroids()[(*(i_fst))*3+div_cInd] >= div_val);

#ifndef NDEBUG
  {
    std::vector<int> test(pIndsBegin, pIndsEnd);
    std::vector<int>::const_iterator ic_pInd;
      
    for (ic_pInd = pIndsBegin; ic_pInd < i_fst; ic_pInd++) {
      assert(this->GetPrimitiveCentroids()[(*ic_pInd)*3+div_cInd] < div_val);
      assert(*ic_pInd < this->GetMesh().GetNumberOfFacets());
    }
    for (; ic_pInd < pIndsEnd; ic_pInd++) {
      assert(this->GetPrimitiveCentroids()[(*ic_pInd)*3+div_cInd] >= div_val);
      assert(*ic_pInd < this->GetMesh().GetNumberOfFacets());
    }
  }
#endif

  p_childIndexBounds[0][0] = pIndsBegin;
  p_childIndexBounds[0][1] = i_fst;
  p_childIndexBounds[1][0] = i_fst;
  p_childIndexBounds[1][1] = pIndsEnd;    
} /* SplitBV */

template <class TBVH>
void tledStaticBVHCreator<TBVH>::GenerateMain() {
  typename TBVH::BoundingVolume root;

  root.ParentIndex = -1;
  this->GetOutput().GetBVs().push_back(root);

  this->GetOutput().GetBVs().reserve(2*this->GetMesh().GetNumberOfFacets() + 1);
  this->InitialiseTopDownRecursive(0, this->GetBVPrimitiveIndices().begin(), this->GetBVPrimitiveIndices().end());
}
