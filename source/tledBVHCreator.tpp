#ifndef NDEBUG
template <class TMesh>
bool CheckConnectivity(const std::vector<int> &primIndList, const TMesh &mesh) {
  typedef typename TMesh::Facet __Facet;

  std::vector<int> uniqueInds;
  std::vector<int>::const_iterator ic_fInd;

  uniqueInds = tledHelper::MakeSortedUnique(primIndList);
  if (uniqueInds.size() == 1) return true;
  for (ic_fInd = uniqueInds.begin(); ic_fInd < uniqueInds.end(); ic_fInd++) {
    const int *tVtxInds = mesh.GetFacet(*ic_fInd).NodeIndices;

    std::vector<int> otherPrimInds;
    std::vector<int> otherNodeInds;
    int vInd;

    otherPrimInds.reserve(uniqueInds.size() - 1);
    std::copy(std::vector<int>::const_iterator(uniqueInds.begin()), ic_fInd, std::back_inserter(otherPrimInds));
    std::copy(ic_fInd + 1, std::vector<int>::const_iterator(uniqueInds.end()), std::back_inserter(otherPrimInds));
    assert(otherPrimInds.size() == uniqueInds.size() - 1);

    otherNodeInds = mesh.CompileNodeListFromFacetList(otherPrimInds.begin(), otherPrimInds.end());
    for (vInd = 0; vInd < __Facet::NumberOfVertices && std::find(otherNodeInds.begin(), otherNodeInds.end(), tVtxInds[vInd]) == otherNodeInds.end(); vInd++);
    if (vInd == __Facet::NumberOfVertices) {
      std::cerr << "Didn't find any of (", std::copy(tVtxInds, tVtxInds + 3, std::ostream_iterator<int>(std::cerr, " ")), std::cerr << ") in:\n";
      std::copy(otherNodeInds.begin(), otherNodeInds.end(), std::ostream_iterator<int>(std::cerr, " ")), std::cerr << std::endl;
      break;
    }
  }

  return ic_fInd == uniqueInds.end();
}
#endif

template <class TBVH, class TMesh>
int tledBVHCreatorImpl<TBVH, TMesh>::AddBV(const int parentIndex) {
  BoundingVolume newBV;

  newBV.ParentIndex = parentIndex;
#ifndef NDEBUG
  std::fill(newBV.ChildIndices, newBV.ChildIndices + BoundingVolume::NumberOfChildBVs, -1);
  newBV.PrimitiveIndex = -1;
#endif
  this->GetOutput().GetBVs().push_back(newBV);

  return this->GetOutput().GetBVs().size() - 1;
}

template <class TBVH, class TMesh>
void tledBVHCreatorImpl<TBVH, TMesh>::InitialiseTopDownRecursive(const int bvInd, const std::vector<int>::iterator pIndsBegin, const std::vector<int>::iterator pIndsEnd) {
  std::vector<int>::iterator i_childpInds[TBVH::BVHOrder][2];

  assert(bvInd < (int)this->GetOutput().GetBVs().size() && bvInd >= 0);
  
  if (pIndsEnd - pIndsBegin == 1) {
    this->GetOutput().GetBV(bvInd).PrimitiveIndex = *pIndsBegin;
    std::fill(this->GetOutput().GetBV(bvInd).ChildIndices, this->GetOutput().GetBV(bvInd).ChildIndices + TBVH::BVHOrder, -1);
    this->GetOutput().RefitLeafBV(bvInd);
    this->GetOutput().GetLeafBVIndices().push_back(bvInd);
    this->GetOutput().GetAllPrimitiveBVIndices()[this->GetOutput().GetBV(bvInd).PrimitiveIndex] = bvInd;

    return;
  } 

  tledLogDebugStream(tledHelper::Info() << "Performing top-down split for primitive set of size " << int(pIndsEnd - pIndsBegin));
  this->GetOutput().GetBV(bvInd).PrimitiveIndex = -1;
  this->GetOutput().ComputeBoundsFromPrimitives(this->GetOutput().GetBV(bvInd), &*pIndsBegin, &*pIndsEnd);
  assert(pIndsEnd - pIndsBegin > 1);

  this->SplitBV(i_childpInds, bvInd, pIndsBegin, pIndsEnd);    

#ifndef NDEBUG
  {
    using namespace tledHelper;

    std::vector<int> primIndCheckSet, primIndRefSet(pIndsBegin, pIndsEnd);

    assert(primIndRefSet.size() > 0);
    for (int cInd = 0; cInd < BoundingVolume::NumberOfChildBVs; cInd++) {
      assert(i_childpInds[cInd][1] <= pIndsEnd);
      assert(i_childpInds[cInd][1] - i_childpInds[cInd][0] >= 0);
    }

    for (int c0Ind = 0; c0Ind < BoundingVolume::NumberOfChildBVs; c0Ind++) {
      std::vector<int> subSet0(i_childpInds[c0Ind][0], i_childpInds[c0Ind][1]);

      subSet0 = MakeSortedUnique(subSet0);
      primIndCheckSet.insert(primIndCheckSet.end(), i_childpInds[c0Ind][0], i_childpInds[c0Ind][1]);
      for (int c1Ind = c0Ind + 1; c1Ind < BoundingVolume::NumberOfChildBVs; c1Ind++) {
        std::vector<int> subSet1(i_childpInds[c1Ind][0], i_childpInds[c1Ind][1]);
	std::vector<int> intersection;

	subSet1 = MakeSortedUnique(subSet1);
        std::set_intersection(subSet0.begin(), subSet0.end(), subSet1.begin(), subSet1.end(), std::back_inserter(intersection));
        assert(intersection.size() == 0);
      }
    }

    assert(primIndCheckSet.size() == MakeSortedUnique(primIndCheckSet).size());
    primIndCheckSet = MakeSortedUnique(primIndCheckSet);

    assert(primIndRefSet.size() == MakeSortedUnique(primIndRefSet).size());
    primIndRefSet = MakeSortedUnique(primIndRefSet);

    assert(primIndRefSet.size() == primIndCheckSet.size());
    assert(std::equal(primIndRefSet.begin(), primIndRefSet.end(), primIndCheckSet.begin()));
  }
#endif

  /*
   * Improve look-up speed by having children next to parent in array (always descend in left tree first!) ->
   * first add children, then subdivide!
   */
  for (int cInd = 0; cInd < TBVH::BVHOrder; cInd++) {
    if (TBVH::BVHOrder == 2 || i_childpInds[cInd][1] - i_childpInds[cInd][0] >= 1) {
      const int newChildInd = this->AddBV(bvInd);

      this->GetOutput().GetBV(bvInd).ChildIndices[cInd] = newChildInd;     
      assert(i_childpInds[cInd][1] - i_childpInds[cInd][0] >= 1);
      this->InitialiseTopDownRecursive(newChildInd, i_childpInds[cInd][0], i_childpInds[cInd][1]);
      assert(this->GetOutput().GetBV(bvInd).ChildIndices[cInd] >= 0 && this->GetOutput().GetBV(bvInd).ChildIndices[cInd] < (int)this->GetOutput().GetBVs().size());
    } else this->GetOutput().GetBV(bvInd).ChildIndices[cInd] = -1;
  }

  this->GetOutput().RefitInteriorBV(bvInd);

  assert(this->GetOutput().GetBV(bvInd).PrimitiveIndex == -1 || (this->GetOutput().GetBV(bvInd).PrimitiveIndex >= 0 && this->GetOutput().GetBV(bvInd).PrimitiveIndex < (int)this->GetMesh().GetNumberOfFacets()));
  assert(this->GetOutput().GetBV(bvInd).ParentIndex == -1 || (this->GetOutput().GetBV(bvInd).ParentIndex >= 0 && this->GetOutput().GetBV(bvInd).ParentIndex < bvInd));
} 

template <class TBVH, class TMesh>
void tledBVHCreatorImpl<TBVH, TMesh>::InitBVHCreationData() {
  m_PrimitiveCentroids.resize(3*this->GetMesh().GetNumberOfFacets());
  this->ComputeCentroids(&m_PrimitiveCentroids.front());

  m_BVPrimitiveIndices = tledSequenceGenerator::MakeSequence(0, this->GetMesh().GetNumberOfFacets());
}

template <class TBVH, class TMesh>
void tledBVHCreatorImpl<TBVH, TMesh>::FreeBVHCreationData() {
  m_PrimitiveCentroids.resize(0);
  m_BVPrimitiveIndices.resize(0);
}

template <class TBVH, class TMesh>
void tledBVHCreatorImpl<TBVH, TMesh>::SetBVHBuffer(TBVH &r_dst) { 
  mp_BVH = &r_dst; 
  this->SetMesh(r_dst.GetMesh());
}

template <class TBVH, class TMesh>
void tledBVHCreatorImpl<TBVH, TMesh>::InitBVHStruct() {
  this->GetOutput().GetAllPrimitiveBVIndices().clear();
  this->GetOutput().GetAllPrimitiveBVIndices().insert(this->GetOutput().GetAllPrimitiveBVIndices().end(), this->GetMesh().GetNumberOfFacets(), -1);
  this->GetOutput().GetBVs().reserve(this->GetMesh().GetNumberOfFacets());
  this->GetOutput().GetLeafBVIndices().reserve(this->GetMesh().GetNumberOfFacets());
#ifndef __CUDACC__
  assert(!std::isnan(this->GetOutput().GetMargin()));
#endif
}

template <class TBVH, class TMesh>
void tledBVHCreatorImpl<TBVH, TMesh>::Generate() {
  this->InitBVHCreationData();
  this->InitBVHStruct();
  this->GenerateMain();
  this->FreeBVHCreationData();
} /* _InitialiseTopDown */

template <class TBVH, class TMesh>
void tledBVHCreatorImpl<TBVH, TMesh>::InitLeafData() {
  assert(this->GetOutput().GetLeafBVIndices().size() > 0);
#ifndef NDEBUG
  std::fill(this->GetOutput().GetAllPrimitiveBVIndices().begin(), this->GetOutput().GetAllPrimitiveBVIndices().end(), -1);
#endif

  for (std::vector<int>::const_iterator ic_leafInd = this->GetOutput().GetLeafBVIndices().begin(); ic_leafInd < this->GetOutput().GetLeafBVIndices().end(); ic_leafInd++) {
    this->GetOutput().GetAllPrimitiveBVIndices()[this->GetOutput().GetBV(*ic_leafInd).PrimitiveIndex] = *ic_leafInd;
  }

  assert(std::find(this->GetOutput().GetAllPrimitiveBVIndices().begin(), this->GetOutput().GetAllPrimitiveBVIndices().end(), -1) == this->GetOutput().GetAllPrimitiveBVIndices().end());
}

template <class TBVH, class TMesh>
void tledBVHCreatorImpl<TBVH, TMesh>::ComputeCentroids(float *p_dst) const {
  int pInd;
  float *p_currCent;

  for (pInd = 0, p_currCent = p_dst; pInd < this->GetMesh().GetNumberOfFacets(); pInd++, p_currCent += 3) {
    this->GetMesh().ComputeCentroid(p_currCent, pInd);
  }
} /* ComputeCentroids */
