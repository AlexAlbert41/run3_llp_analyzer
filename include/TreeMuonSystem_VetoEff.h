// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information
// changed by Alex to be better equipped for skimmed/merged nTuples with both cluster and muon level information

#ifndef TreeMuonSystem_VetoEff_H
#define TreeMuonSystem_VetoEff_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define N_MAX_DTRECHITS 20000
#define NTriggersMAX 1201 // Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20
#define N_MAX_GPARTICLES 5000
#define N_MAX_LLP 200
#define N_MAX_MATCHED_JETS 10
#define N_MAX_MATCHED_MUONS 10

#include <iostream>
#include <string>
#include <sys/stat.h>
#include "assert.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TTree.h"
#include "DBSCAN.h"

#include "RazorAnalyzer_TnP.h"

#include "RazorHelper.h"

class TreeMuonSystem_VetoEff
{

public:
  TreeMuonSystem_VetoEff();
  ~TreeMuonSystem_VetoEff();
  // TreeMuonSystemCombination_TnP::TreeMuonSystemCombination_TnP()
  // {
  //   InitVariables();
  // };
  // TreeMuonSystemCombination_TnP::~TreeMuonSystemCombination_TnP()
  // {
  //   if (f_) f_->Close();
  // };
  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum, MC_condition, category;
  UInt_t  npv, npu;
  float rho, weight;
  float pileupWeight;
  float pileupWeightUp; float pileupWeightDown;

  float met, metPhi;
  float puppiMet, puppiMetPhi;
  bool Flag_HBHENoiseFilter, Flag_HBHEIsoNoiseFilter, Flag_BadPFMuonFilter, Flag_BadPFMuonDzFilter, Flag_hfNoisyHitsFilter, Flag_globalSuperTightHalo2016Filter,
  Flag_CSCTightHaloFilter, Flag_BadChargedCandidateFilter, Flag_eeBadScFilter, Flag_goodVertices, Flag_ecalBadCalibFilter, Flag_all;
  int mH, mX, ctau;

  bool Flag2_HBHENoiseFilter, Flag2_HBHEIsoNoiseFilter, Flag2_BadPFMuonFilter, Flag2_globalSuperTightHalo2016Filter,
  Flag2_globalTightHalo2016Filter, Flag2_BadChargedCandidateFilter, Flag2_EcalDeadCellTriggerPrimitiveFilter,
  Flag2_ecalBadCalibFilter, Flag2_eeBadScFilter, Flag2_all, Flag_EcalDeadCellTriggerPrimitiveFilter, jetVeto;

  int numProbeMuons; int numTag;

  float sampledEta, sampledPhi;
  bool sampledCSC, sampledDT;
  int numMatchedJets_CSC, numMatchedJets_DT, numMatchedMuons_CSC, numMatchedMuons_DT_lowMET, numMatchedMuons_DT_highMET, numMatchedRechits_CSC, numMatchedRechits_DT;

  //csc
  int           nCscRechits;
  int           nCscRings;
  int           nDTRechits;
  int           nDtRings;


  float         matchedJetPt_CSC[N_MAX_MATCHED_JETS];
  float         matchedJetEta_CSC[N_MAX_MATCHED_JETS];
  float         matchedJetPhi_CSC[N_MAX_MATCHED_JETS];
  float         matchedJetPt_DT[N_MAX_MATCHED_JETS];
  float         matchedJetEta_DT[N_MAX_MATCHED_JETS];
  float         matchedJetPhi_DT[N_MAX_MATCHED_JETS];
  float         matchedMuonPt_CSC[N_MAX_MATCHED_MUONS];
  float         matchedMuonEta_CSC[N_MAX_MATCHED_MUONS];
  float         matchedMuonPhi_CSC[N_MAX_MATCHED_MUONS];
  float         matchedMuonPt_DT_lowMET[N_MAX_MATCHED_MUONS];
  float         matchedMuonEta_DT_lowMET[N_MAX_MATCHED_MUONS];
  float         matchedMuonPhi_DT_lowMET[N_MAX_MATCHED_MUONS];
  float         matchedMuonPt_DT_highMET[N_MAX_MATCHED_MUONS];
  float         matchedMuonEta_DT_highMET[N_MAX_MATCHED_MUONS];
  float         matchedMuonPhi_DT_highMET[N_MAX_MATCHED_MUONS];

  int           nDtRechitClusters;
  bool           dtRechitClusterOverlap[N_MAX_CSC];
  int           dtRechitClusterNSegStation1[N_MAX_CSC];
  int           dtRechitClusterNSegStation2[N_MAX_CSC];
  int           dtRechitClusterNSegStation3[N_MAX_CSC];
  int           dtRechitClusterNSegStation4[N_MAX_CSC];
  int           dtRechitClusterNOppositeSegStation1[N_MAX_CSC];
  int           dtRechitClusterNOppositeSegStation2[N_MAX_CSC];
  int           dtRechitClusterNOppositeSegStation3[N_MAX_CSC];
  int           dtRechitClusterNOppositeSegStation4[N_MAX_CSC];
  int           dtRechitClusterNHitStation1[N_MAX_CSC];
  int           dtRechitClusterNHitStation2[N_MAX_CSC];
  int           dtRechitClusterNHitStation3[N_MAX_CSC];
  int           dtRechitClusterNHitStation4[N_MAX_CSC];


  int           dtRechitCluster_match_MB1hits_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_cosmics_plus[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_cosmics_minus[N_MAX_CSC];
  int           dtRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_MB1Seg_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_RPCBx_dPhi0p5[N_MAX_CSC];
  int         dtRechitCluster_match_RPChits_dPhi0p5[N_MAX_CSC];

  int           dtRechitCluster_match_RB1_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_RB1_dPhi0p5[N_MAX_CSC];



  bool          dtRechitCluster_match_gLLP[N_MAX_CSC];
  int           dtRechitCluster_match_gLLP_index[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_minDeltaR[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_eta[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_phi[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_decay_r[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_decay_z[N_MAX_CSC];
  bool         dtRechitCluster_match_gLLP_csc[N_MAX_CSC];
  bool         dtRechitCluster_match_gLLP_dt[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_e[N_MAX_CSC];

  float         dtRechitClusterX[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterY[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterZ[N_MAX_CSC];   //[nCsc]

  int         dtRechitClusterWheel[N_MAX_CSC];

  float         dtRechitClusterEta[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterPhi[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterSize[N_MAX_CSC];
  int           dtRechitClusterNoiseHit[N_MAX_CSC];
  int           dtRechitClusterNoiseHitStation1[N_MAX_CSC];
  int           dtRechitClusterNoiseHitStation2[N_MAX_CSC];
  int           dtRechitClusterNoiseHitStation3[N_MAX_CSC];
  int           dtRechitClusterNoiseHitStation4[N_MAX_CSC];


  float         dtRechitClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterMaxStation[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterNStation10[N_MAX_CSC];
  float          dtRechitClusterAvgStation10[N_MAX_CSC];
  float         dtRechitClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterNChamber[N_MAX_CSC];

  float         dtRechitClusterJetVetoPt[N_MAX_CSC];
  float         dtRechitClusterJetVetoE[N_MAX_CSC];
  bool         dtRechitClusterJetVetoLooseId[N_MAX_CSC];
  bool         dtRechitClusterJetVetoTightId[N_MAX_CSC];


  float         dtRechitClusterMuonVetoPt[N_MAX_CSC];
  float         dtRechitClusterMuonVetoE[N_MAX_CSC];

  bool          dtRechitClusterMuonVetoTightId[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoGlobal[N_MAX_CSC];

  float         dtRechitClusterMet_dPhi[N_MAX_CSC];

  bool          dtRechitCluster_matchToProbeMuon[N_MAX_CSC];
  
  int         dtRechitClusternXY[N_MAX_CSC];
  int         dtRechitClusternZ[N_MAX_CSC];
  float         dtRechitClusterXSpread[N_MAX_CSC];
  float         dtRechitClusterYSpread[N_MAX_CSC];
  float         dtRechitClusterZSpread[N_MAX_CSC];
  float         dtRechitClusterXYSpread[N_MAX_CSC];
  float         dtRechitClusterRSpread[N_MAX_CSC];
  float         dtRechitClusterEtaPhiSpread[N_MAX_CSC];
  float         dtRechitClusterEtaSpread[N_MAX_CSC];
  float         dtRechitClusterPhiSpread[N_MAX_CSC];
  float         dtRechitClusterDeltaRSpread[N_MAX_CSC];
  float         dtRechitClusterMajorAxis[N_MAX_CSC];
  float         dtRechitClusterMinorAxis[N_MAX_CSC];
  float         dtRechitClusterSkewX[N_MAX_CSC];
  float         dtRechitClusterSkewY[N_MAX_CSC];
  float         dtRechitClusterSkewZ[N_MAX_CSC];
  float         dtRechitClusterKurtX[N_MAX_CSC];
  float         dtRechitClusterKurtY[N_MAX_CSC];
  float         dtRechitClusterKurtZ[N_MAX_CSC];


  bool           dtRechitCluster_matchToMuon1[N_MAX_CSC];
  bool           dtRechitCluster_matchToMuon2[N_MAX_CSC];
  bool          dtRechitCluster_matchToNotProbeMuon[N_MAX_CSC];
  bool          dtRechitCluster_matchToProbeAndJet[N_MAX_CSC];
  bool          dtRechitCluster_matchToHighPtJet[N_MAX_CSC];
  bool          dtRechitCluster_matchToLowPtJet[N_MAX_CSC];
  bool          dtRechitCluster_notMatched[N_MAX_CSC];
  bool          dtRechitCluster_PassTimeVeto[N_MAX_CSC];
  bool          dtRechitCluster_passForwardVeto[N_MAX_CSC];
  bool          dtRechitCluster_HLTCscCluster_Loose_Decision[N_MAX_CSC];
  bool          dtRechitCluster_HLTCscCluster_Medium_Decision[N_MAX_CSC];
  bool          dtRechitCluster_HLTCscCluster_Tight_Decision[N_MAX_CSC];
  int           dtRechitCluster_match_gParticle_id[N_MAX_CSC];
  int           dtRechitCluster_match_gParticle_mother_id[N_MAX_CSC];
  float         dtRechitCluster_match_gParticle_cluster_deltaR[N_MAX_CSC];
  float         dtRechitCluster_match_gParticle_pt[N_MAX_CSC];
  
  
  
  int           nCscRechitClusters;


  bool          cscRechitCluster_match_gLLP[N_MAX_CSC];
  int           cscRechitCluster_match_gLLP_index[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_minDeltaR[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_eta[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_phi[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_decay_r[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_decay_z[N_MAX_CSC];
  bool         cscRechitCluster_match_gLLP_csc[N_MAX_CSC];
  bool         cscRechitCluster_match_gLLP_dt[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_e[N_MAX_CSC];



  float         cscRechitClusterX[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterY[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterZ[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterTimeWeighted[N_MAX_CSC];
  float         cscRechitClusterTimeSpreadWeightedAll[N_MAX_CSC];
  float         cscRechitClusterTime[N_MAX_CSC];
  float         cscRechitClusterTimeSpread[N_MAX_CSC];
  float         cscRechitClusterEta[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterPhi[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterSize[N_MAX_CSC];
  int           cscRechitClusternXY[N_MAX_CSC];
  int           cscRechitClusternZ[N_MAX_CSC];
  float         cscRechitClusterXSpread[N_MAX_CSC];
  float         cscRechitClusterYSpread[N_MAX_CSC];
  float         cscRechitClusterZSpread[N_MAX_CSC];
  float         cscRechitClusterXYSpread[N_MAX_CSC];
  float         cscRechitClusterRSpread[N_MAX_CSC];
  float         cscRechitClusterEtaPhiSpread[N_MAX_CSC];
  float         cscRechitClusterEtaSpread[N_MAX_CSC];
  float         cscRechitClusterPhiSpread[N_MAX_CSC];
  float         cscRechitClusterDeltaRSpread[N_MAX_CSC];
  float         cscRechitClusterMajorAxis[N_MAX_CSC];
  float         cscRechitClusterMinorAxis[N_MAX_CSC];
  float         cscRechitClusterSkewX[N_MAX_CSC];
  float         cscRechitClusterSkewY[N_MAX_CSC];
  float         cscRechitClusterSkewZ[N_MAX_CSC];
  float         cscRechitClusterKurtX[N_MAX_CSC];
  float         cscRechitClusterKurtY[N_MAX_CSC];
  float         cscRechitClusterKurtZ[N_MAX_CSC];


  float         cscRechitClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterMaxStation[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterNStation10[N_MAX_CSC];
  float          cscRechitClusterAvgStation10[N_MAX_CSC];
  float         cscRechitClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterNChamber[N_MAX_CSC];

  float         cscRechitClusterJetVetoPt[N_MAX_CSC];
  bool         cscRechitClusterJetVetoLooseId[N_MAX_CSC];
  bool         cscRechitClusterJetVetoTightId[N_MAX_CSC];
  float         cscRechitClusterJetVetoE[N_MAX_CSC];

  float         cscRechitClusterMuonVetoPt[N_MAX_CSC];
  float         cscRechitClusterMuonVetoE[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoGlobal[N_MAX_CSC];

  int           cscRechitCluster_match_dtSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_RE12_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_RB1_0p4[N_MAX_CSC];


  int           cscRechitClusterNRechitChamberPlus11[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus12[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus13[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus21[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus22[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus31[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus32[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus41[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus42[N_MAX_CSC];

  int           cscRechitClusterNRechitChamberMinus11[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus12[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus13[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus21[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus22[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus31[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus32[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus41[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus42[N_MAX_CSC];
  float         cscRechitClusterMet_dPhi[N_MAX_CSC];
  bool          cscRechitCluster_matchToProbeMuon[N_MAX_CSC];
  bool           cscRechitCluster_matchToMuon1[N_MAX_CSC];
  bool          cscRechitCluster_matchToMuon2[N_MAX_CSC];
  bool          cscRechitCluster_matchToNotProbeMuon[N_MAX_CSC];
  bool          cscRechitCluster_matchToProbeAndJet[N_MAX_CSC];
  bool          cscRechitCluster_matchToHighPtJet[N_MAX_CSC];
  bool          cscRechitCluster_matchToLowPtJet[N_MAX_CSC];
  bool          cscRechitCluster_notMatched[N_MAX_CSC];
  bool          cscRechitCluster_PassTimeVeto[N_MAX_CSC];
  bool          cscRechitCluster_passME1112Veto[N_MAX_CSC];
  bool          cscRechitCluster_HLTCscCluster_Loose_Decision[N_MAX_CSC];
  bool          cscRechitCluster_HLTCscCluster_Medium_Decision[N_MAX_CSC];
  bool          cscRechitCluster_HLTCscCluster_Tight_Decision[N_MAX_CSC];
  int           cscRechitCluster_match_gParticle_id[N_MAX_CSC];
  int           cscRechitCluster_match_gParticle_mother_id[N_MAX_CSC];
  float         cscRechitCluster_match_gParticle_cluster_deltaR[N_MAX_CSC];
  float         cscRechitCluster_match_gParticle_pt[N_MAX_CSC];


  //gLLP
  int nGLLP;
  float gLLP_eta[N_MAX_LLP];
  float gLLP_phi[N_MAX_LLP];
  float gLLP_csc[N_MAX_LLP];
  float gLLP_dt[N_MAX_LLP];
  float gLLP_beta[N_MAX_LLP];
  float gLLP_e[N_MAX_LLP];
  float gLLP_pt[N_MAX_LLP];
  float gLLP_ctau[N_MAX_LLP];
  float gLLP_decay_vertex_r[N_MAX_LLP];
  float gLLP_decay_vertex_x[N_MAX_LLP];
  float gLLP_decay_vertex_y[N_MAX_LLP];
  float gLLP_decay_vertex_z[N_MAX_LLP];

  float gHiggsPt;
  float gHiggsEta;
  float gHiggsPhi;
  float gHiggsE;

  //leptons

  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int  lepPdgId[N_MAX_LEPTONS];
  float lepDZ[N_MAX_LEPTONS];

  bool lepTightId[N_MAX_LEPTONS];
  bool lepPassLooseIso[N_MAX_LEPTONS];
  bool lepPassTightIso[N_MAX_LEPTONS];
  bool lepPassVTightIso[N_MAX_LEPTONS];
  bool lepPassVVTightIso[N_MAX_LEPTONS];

  //lep variables (added from https://github.com/cms-lpc-llp/llp_analyzer/blob/master/include/LiteTreeMuonSystem.h)
  float lepTriggerSF[N_MAX_LEPTONS];
  float lepTightIdSF[N_MAX_LEPTONS];
  float lepLooseIdSF[N_MAX_LEPTONS];
  float lepTightIsoSF[N_MAX_LEPTONS];
  float lepLooseIsoSF[N_MAX_LEPTONS];
  float lepTriggerMCEfficiency[N_MAX_LEPTONS];
  float lepTightIdMCEfficiency[N_MAX_LEPTONS];
  float lepLooseIdMCEfficiency[N_MAX_LEPTONS];
  float lepTightIsoMCEfficiency[N_MAX_LEPTONS];
  float lepLooseIsoMCEfficiency[N_MAX_LEPTONS];
  float lepEff[N_MAX_LEPTONS];
  float lepSF[N_MAX_LEPTONS];
  bool lepTag[N_MAX_LEPTONS];

  bool lepLoosePassId[N_MAX_LEPTONS];
  bool lepMediumPassId[N_MAX_LEPTONS];
  bool lepTightPassId[N_MAX_LEPTONS];
  bool lepPassVetoId[N_MAX_LEPTONS];
  bool lepFromZ[N_MAX_LEPTONS];
  bool lepPassId[N_MAX_LEPTONS];
  float lepOverallSF;
  bool  lepPassProbe[N_MAX_LEPTONS];

  //Z-candidate (added from https://github.com/cms-lpc-llp/llp_analyzer/blob/master/include/LiteTreeMuonSystem.h)
  float MT;
  float ZMass1;
  float ZMass;
  float ZPt;
  float ZEta;
  float ZPhi;
  int ZleptonIndex1;
  int ZleptonIndex2;
  
  //jets
  int nJets;
  float jetE[N_MAX_JETS];
  float jetPt[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];


  bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision;
  bool HLT_CscCluster_Loose;
  bool HLT_CscCluster_Medium;
  bool HLT_CscCluster_Tight;
  bool HLT_IsoMu20;
  bool HLT_IsoMu24;
  bool L1_SingleMuShower_Nominal;
  bool L1_SingleMuShower_Tight;




  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
