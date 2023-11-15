template <class TBVH>
void tledSelfCollisionBVHXMLImporter<TBVH>::ParseAdjacencyData() {
  XMLNode node;

  node = this->GetUniqueChild("NonAdjacentGeometryNodes", false);
  if (!node.isEmpty()) this->GetOutput().GetNonAdjacentGeometryNodes() = GetXMLTextAsVector<int>(node);
  
  node = this->GetUniqueChild("GeometryClusterSubtreeRootIndices", false);
  if (!node.isEmpty()) this->GetOutput().GetGeometryClusterSubtreeRootIndices() = GetXMLTextAsVector<int>(node);
}

template <class TBVH>
void tledSelfCollisionBVHXMLImporter<TBVH>::ParseCollisionCandidates() {
  this->GetOutput().GetSelfCollisionCandidates() = GetXMLTextAsVector<int>(this->GetUniqueChild("SelfCollisionCandidates", true));
}

template <class TBVH>
void tledSelfCollisionBVHXMLImporter<TBVH>::Import() {
  this->GetOutput().ResetUpdateCounter();

  Superclass::Import();
  this->ParseAdjacencyData();
  this->ParseCollisionCandidates();
}
