// =========================================================================
// File:       tledBVH.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TContactMesh, class TBV, class TAPI>
void tledBVHImpl<TContactMesh, TBV, TAPI>::CompilePrimitiveListRecursive(std::vector<int> &r_dstList, const int subtreeRootAABBInd) const {
  const BoundingVolume &subtreeRoot = this->GetBVs()[subtreeRootAABBInd];

  if (subtreeRoot.PrimitiveIndex >= 0) {
    r_dstList.push_back(subtreeRoot.PrimitiveIndex);
  } else {
    int const *pc_childInd;

    for (pc_childInd = subtreeRoot.ChildIndices; pc_childInd < subtreeRoot.ChildIndices + TBV::NumberOfChildBVs; pc_childInd++) if (BVHOrder == 2 || *pc_childInd >= 0) {
	this->CompilePrimitiveListRecursive(r_dstList, *pc_childInd);
      }
  }
} 

template <class TContactMesh, class TBV, class TAPI>
void tledBVHImpl<TContactMesh, TBV, TAPI>::ComputeBoundsFromNodes(BoundingVolume &r_bv, const int *nodeListStart, const int *nodeListEnd) const {
  r_bv.ComputeFromNodeList(nodeListStart, nodeListEnd, this->GetMesh().GetAllNodeCoordinates());
}

template <class TContactMesh, class TBV, class TAPI>
void tledBVHImpl<TContactMesh, TBV, TAPI>::ComputeBoundsFromChildren(BoundingVolume &r_bv) const {
  assert(r_bv.ChildIndices[0] >= 0 && r_bv.ChildIndices[0] < this->GetNumberOfBVs());
  TBV::CopyBoundsFromBV(r_bv, this->GetBV(r_bv.ChildIndices[0]));
  for (int const *pc_childInd = r_bv.ChildIndices + 1; pc_childInd < r_bv.ChildIndices + BVHOrder; pc_childInd++) if (BVHOrder == 2 || *pc_childInd >= 0) {
      assert(*pc_childInd >= 0 && *pc_childInd < this->GetNumberOfBVs());
      TBV::Merge(r_bv, this->GetBV(*pc_childInd));
    }
}

template <class TContactMesh, class TBV, class TAPI>
void tledBVHImpl<TContactMesh, TBV, TAPI>::ComputeBoundsFromPrimitives(BoundingVolume &r_bv, const int *primitiveListStart, const int *primitiveListEnd) const {
  const std::vector<int> nodeList = this->GetMesh().CompileNodeListFromFacetList(primitiveListStart, primitiveListEnd);

  this->ComputeBoundsFromNodes(r_bv, &nodeList.front(), &nodeList.back() + 1);
}

template <class TContactMesh, class TBV, class TAPI>
void tledBVHImpl<TContactMesh, TBV, TAPI>::RefitLeafBV(const int bvInd) {
  BoundingVolume &r_bv = this->GetBV(bvInd);

  assert(r_bv.PrimitiveIndex >= 0 && r_bv.PrimitiveIndex < GetMesh().GetNumberOfFacets());
  assert(r_bv.ChildIndices[0] < 0);
  this->ComputeBoundsFromNodes(r_bv, this->GetMesh().GetFacet(r_bv.PrimitiveIndex).NodeIndices, this->GetMesh().GetFacet(r_bv.PrimitiveIndex).NodeIndices + Facet::NumberOfVertices);
  r_bv.AddMargin(this->GetMargin());
  this->RefitBVCommon(bvInd);
}

template <class TContactMesh, class TBV, class TAPI>
void tledBVHImpl<TContactMesh, TBV, TAPI>::RefitInteriorBV(const int bvInd) {
  BoundingVolume &r_bv = this->GetBV(bvInd);

  assert(r_bv.PrimitiveIndex < 0);
  assert(r_bv.ChildIndices[0] >= 0);
  this->ComputeBoundsFromChildren(r_bv);
  this->RefitBVCommon(bvInd);
}

template <class TContactMesh, class TBV, class TAPI>
void tledBVHImpl<TContactMesh, TBV, TAPI>::Init(tledBVHCreator &r_bvhBuilder) {
  BVHBuilder &r_specBuilder = static_cast<BVHBuilder&>(r_bvhBuilder);

  r_specBuilder.SetBVHBuffer(*this);
  r_specBuilder.SetMesh(this->GetMesh());
  r_specBuilder.Generate();
}

template <class TContactMesh, class TBV, class TAPI>
void tledBVHImpl<TContactMesh, TBV, TAPI>::GetSubtreeLeafs(std::vector<int> &r_leafIndexBuffer, const int rootIndex) const {
  if (this->IsLeaf(rootIndex)) {
    r_leafIndexBuffer.push_back(rootIndex);
  } else {
    const BoundingVolume &bv = this->GetBV(rootIndex);

    for (int const *pc_cInd = bv.ChildIndices; pc_cInd < bv.ChildIndices + BVHOrder && (BVHOrder == 2 || *pc_cInd >= 0); pc_cInd++) {
      this->GetSubtreeLeafs(r_leafIndexBuffer, *pc_cInd);
    }
  }  
}

template <class TContactMesh, class TBV, class TAPI>
int tledBVHImpl<TContactMesh, TBV, TAPI>::GetSubtreeMaxDepth(const int rootIndex) const {
  if (this->IsLeaf(rootIndex)) return 0;
  else {
    const BoundingVolume &bv = this->GetBV(rootIndex);

    int depth = this->GetSubtreeMaxDepth(*bv.ChildIndices);

    for (int const *pc_c = bv.ChildIndices + 1; pc_c < bv.ChildIndices + BoundingVolume::NumberOfChildBVs; pc_c++) if (BVHOrder == 2 || *pc_c >= 0) {
	depth = std::max(depth, this->GetSubtreeMaxDepth(*pc_c));
      }

    return depth + 1;
  }
}

template <class TContactMesh, class TBV, class TAPI>
XMLNode tledBVHImpl<TContactMesh, TBV, TAPI>::ExportToXML() const {
  tledBVHXMLExporter<tledBVHImpl> exporter;
  XMLNode root;

  exporter.SetInput(*this);
  exporter.Export();
  root = exporter.GetRootNode();

  return root;
}
