// =========================================================================
// File:       tledBVHXMLImporter.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBVH>
void tledBVHXMLImporter<TBVH>::ParseLeafData() {
  this->GetOutput().GetLeafBVIndices() = GetXMLTextAsVector<int>(this->GetUniqueChild("LeafBVIndices", true));
}

template <class TBVH>
void tledBVHXMLImporter<TBVH>::ParseBaseStatsData() {
  this->GetOutput().SetMargin(this->template GetNumericElementValue<float>(this->GetUniqueChild("Margin", true)));
}

template <class TBVH>
void tledBVHXMLImporter<TBVH>::ParseBVData() {
  const XMLNode bvArrayRoot = this->GetUniqueChild("BoundingVolumes", true);

  tledBVHExporterBVXMLExporterAdapter<BoundingVolume> exporter;
    
  this->GetOutput().GetBVs().resize(this->template GetNumericAttribute<int>(bvArrayRoot, "NumBVs"));
  if (bvArrayRoot.nChildNode(exporter.GetRootElementName()) != this->GetOutput().GetNumberOfBVs()) {
    tledLogErrorStream(tledHelper::FatalError() << "Dimension mismatch between BV array and XML node array: " << this->GetOutput().GetNumberOfBVs() << " - " << bvArrayRoot.nChildNode(exporter.GetRootElementName()));
  }

  for (int b = 0; b < this->GetOutput().GetNumberOfBVs(); b++) {
    tledBVHImporterBVXMLImporterAdapter<BoundingVolume> bvImp;
    XMLNode bvNode = bvArrayRoot.getChildNode(exporter.GetRootElementName(), b);
    
    bvImp.SetRootNode(bvNode);
    bvImp.SetOuputObject(this->GetOutput().GetBV(b));
    bvImp.Import();
  }
}

template <class TBVH>
void tledBVHXMLImporter<TBVH>::Import() {
  this->ParseBaseStatsData();
  this->ParseLeafData();
  this->ParseBVData();
  this->GetOutput().LoadFromXMLPostloadHook();
}
