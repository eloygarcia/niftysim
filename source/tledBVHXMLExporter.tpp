// =========================================================================
// File:       tledBVHXMLExporter.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBVH>
void tledBVHXMLExporter<TBVH>::WriteBVHBaseStats() {
  std::string bvGeomType = "INV", bvTypeSfx = "";

  if (this->GetInput().GetBVHOrder() > 2) {
    char numStr[] = {char(this->GetInput().GetBVHOrder() + '0'), '\0'};

    bvTypeSfx = std::string(numStr);
  } 

  switch (this->GetInput().GetBVGeometryID()) {
  case tledAABB<2>::BVGeometryID:
    bvGeomType = "AABB";
    break;

  case tledOBB<2>::BVGeometryID:
    bvGeomType = "OBB";
    break;

  default:
    tledLogErrorStream(tledHelper::FatalError() << "BV type ID " << this->GetInput().GetBVGeometryID() << " not recognised.");
  }

  this->CreateTextNode("BVType", bvGeomType + bvTypeSfx);
  this->CreateNumericNode("Margin", this->GetInput().GetMargin());
}

template <class TBVH>
void tledBVHXMLExporter<TBVH>::WriteLeafData() {
  this->CreateNumericListNode("LeafBVIndices", this->GetInput().GetLeafBVIndices());
}

template <class TBVH>
void tledBVHXMLExporter<TBVH>::WriteBVs() {
  XMLNode bvArrayRoot = this->GetRootNode().addChild("BoundingVolumes");
    
  this->AddNumericAttribute(bvArrayRoot, "NumBVs", this->GetInput().GetNumberOfBVs());
  for (int b = 0; b < this->GetInput().GetNumberOfBVs(); b++) {
    tledBVHExporterBVXMLExporterAdapter<BoundingVolume> bvExp;
    
    bvExp.SetInput(this->GetInput().GetBV(b));
    bvExp.Export();
    bvArrayRoot.addChild(bvExp.GetRootNode());
  }
}

template <class TBVH>
void tledBVHXMLExporter<TBVH>::WriteBody() {
  this->WriteBVHBaseStats();
  this->WriteLeafData();
  this->WriteBVs();
}
