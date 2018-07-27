#pragma once

// Plugin implementation class
class Cpaleosea2d
{
public:
  Cpaleosea2d(IfmDocument pDoc);
  ~Cpaleosea2d();
  static Cpaleosea2d* FromHandle(IfmDocument pDoc);

  int editmenu_nEnum;

#pragma region IFM_Definitions
  // Implementation
public:
  void Serialize (IfmDocument pDoc, IfmArchive pArc);
  void OnEditDocument (IfmDocument pDoc, Widget wParent);
  void PreSimulation (IfmDocument pDoc);
  void PostSimulation (IfmDocument pDoc);
  void PreTimeStep (IfmDocument pDoc);
  void PostTimeStep (IfmDocument pDoc);
  void PreFlowSimulation (IfmDocument pDoc);
  void PostFlowSimulation (IfmDocument pDoc);
  void PreMassSimulation (IfmDocument pDoc, int iSpecies);
  void PostMassSimulation (IfmDocument pDoc, int iSpecies);
  IfmBool OnTimeStepConstraint (IfmDocument pDoc, double tNow, double* dtProposed);
#pragma endregion

private:
  IfmDocument m_pDoc;
};
