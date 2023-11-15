// =========================================================================
// File:       tledMembraneSelfCollisionBVHCreator.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBVH>
int tledMembraneSelfCollisionBVHCreator<TBVH>::_FindMembraneClusterRoot(const int primitiveInd) const {
  int bvInd = this->GetOutput().GetPrimitiveBVIndex(primitiveInd);

  do {
    bvInd = this->GetOutput().GetBV(bvInd).ParentIndex;
    assert(bvInd >= 0);
  } while (!this->GetOutput().GetBV(bvInd).HasDisconnectedGeometry());

  return bvInd;
}

template <class TBVH>
void tledMembraneSelfCollisionBVHCreator<TBVH>::PreBottomUpHook(std::vector<int> &r_activeBVs) {
  const ContactMesh &surface = this->GetMesh();

  std::vector<std::pair<int, int> > membraneClusterPairs;
  std::vector<int> deleteBVs;

  if (this->GetOutput().GetBVHOrder() > 2) {
    tledFatalNotYetImplementedError;
  }

  for (std::vector<int>::const_iterator ic_cs = this->GetClusterStartIndices().begin(); ic_cs < this->GetClusterStartIndices().end(); ic_cs++) {
    /* Based on assumption that membrane clusters only comprise membrane primitives, hence checking first element of each cluster is sufficient */
    if (*(this->GetBVPrimitiveIndices().begin() + *ic_cs) >= surface.GetMembraneFacetBaseIndex() && *(this->GetBVPrimitiveIndices().begin() + *ic_cs) < surface.GetMembraneFacetBaseIndex() + surface.GetNumberOfMembraneElements()) {
      const int basePInd = *(this->GetBVPrimitiveIndices().begin() + *ic_cs);

      std::pair<int, int> inds;

      inds.first = ic_cs - this->GetClusterStartIndices().begin();
#ifndef NDEBUG
      inds.second = -1;
#endif
      m_MembraneSeedFacetInds.push_back(basePInd);

      for (std::vector<int>::const_iterator ic_bcs = this->GetClusterStartIndices().begin(); ic_bcs < this->GetClusterStartIndices().end(); ic_bcs++) if (*(this->GetBVPrimitiveIndices().begin() + *ic_bcs) >= surface.GetMembraneFacetBaseIndex() + surface.GetNumberOfMembraneElements()) {
	  for (std::vector<int>::const_iterator ic_pi = this->GetBVPrimitiveIndices().begin() + *ic_bcs; ic_pi < this->GetBVPrimitiveIndices().end() && (ic_bcs == this->GetClusterStartIndices().end() - 1 || ic_pi < this->GetBVPrimitiveIndices().begin() + *(ic_bcs + 1)); ic_pi++) {
	    if (*ic_pi == basePInd + surface.GetNumberOfMembraneElements()) {
	      inds.second = ic_bcs - this->GetClusterStartIndices().begin();
	      break;
	    }
	  }
	} /* for membrane backside clusters */

      assert(inds.second >= 0 && inds.first >= 0);

      r_activeBVs.push_back(this->InsertParentBinary(inds.first, inds.second));

      deleteBVs.push_back(inds.first);
      deleteBVs.push_back(inds.second);
    } /* if membrane front side cluster */
  }  

  for (std::vector<int>::const_iterator ic_bv = deleteBVs.begin(); ic_bv < deleteBVs.end(); ic_bv++) r_activeBVs.erase(std::find(r_activeBVs.begin(), r_activeBVs.end(), *ic_bv));
}


template <class TBVH>
void tledMembraneSelfCollisionBVHCreator<TBVH>::GenerateMain() {
  Superclass::GenerateMain();

  for (std::vector<int>::const_iterator ic_s = this->m_MembraneSeedFacetInds.begin(); ic_s < this->m_MembraneSeedFacetInds.end(); ic_s++) {
    std::vector<int> &r_nonAdjCheckBVs = this->GetOutput().GetNonAdjacentGeometryNodes();

    assert(_FindMembraneClusterRoot(*ic_s) == _FindMembraneClusterRoot(*ic_s + this->GetMesh().GetNumberOfMembraneElements()));
    assert(!this->GetOutput().GetBV(this->GetOutput().GetBV(_FindMembraneClusterRoot(*ic_s)).ChildIndices[0]).HasDisconnectedGeometry() && !this->GetOutput().GetBV(this->GetOutput().GetBV(_FindMembraneClusterRoot(*ic_s)).ChildIndices[0]).HasDisconnectedGeometry());
    r_nonAdjCheckBVs.erase(std::find(r_nonAdjCheckBVs.begin(), r_nonAdjCheckBVs.end(), _FindMembraneClusterRoot(*ic_s)));
  }
}
