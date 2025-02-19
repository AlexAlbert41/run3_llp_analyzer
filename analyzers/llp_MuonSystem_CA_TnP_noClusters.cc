//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/Selector.hh"
#include "llp_MuonSystem_CA_TnP_noClusters.h"
#include "RazorHelper.h"
#include "RazorAnalyzer_trigEff.h"
#include "TreeMuonSystemCombination_TnP.h"
#include "TreeMuonSystem_Skim_Merge_TnP.h"

#include "CACluster.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <iostream>
#include <random>
#include <bitset>
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

//using namespace fastjet;
using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;

struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  // bool passLooseId;
  // bool passMediumId;
  bool passId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;
  bool passLooseId;
  bool isGlobal;
};


struct jets
{
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  int ecalNRechits;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;
  float highPtJet; //added by Alex
};

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;


void llp_MuonSystem_CA_TnP_noClusters::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option
  //
  //
  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);
  //
  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;

  cout<<"Analysis Tag: "<<analysisTag<<endl;
  bool signalScan = int(options/10) == 1;
  int option = options%10;
  // if (options % 1){
  //   option = 1; // used when running condor
  // }
  // else{
  //   option = 0;// used when running locally
  // }

  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }
  if( signalScan )
  {
    std::cout << "[INFO]: running with Signal scan" << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running without Signal scan " << option << std::endl;
  }



  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "Razor2016_80X";

  }


  const int zh_lepton0_cut = 15;
  const int zh_lepton1_cut = 15;

  const int wh_Muon_pt_cut = 25;
  const int wh_elePt_cut = 35;



  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");


  TreeMuonSystem_Skim_Merge_TnP *MuonSystem = new TreeMuonSystem_Skim_Merge_TnP();
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*>Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  map<pair<int,int>, TH1F*> accep2D;
  map<pair<int,int>, TH1F*> accep_met2D;
  map<pair<int,int>, TH1F*> Total2D;



  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);

  TH1F *probe_cluster_minDeltaR = new TH1F("probe_cluster_minDeltaR", "probe_cluster_minDeltaR", 100, 0, 5);
  TH1F *probe_cluster_minDeltaR_DT = new TH1F("probe_cluster_minDeltaR_DT", "probe_cluster_minDeltaR_DT", 100, 0, 5);

  TH1F *ZMass_Hist = new TH1F("Z Mass", "Z Mass", 100, 50, 120);

  //JetDefinition jet_def( antikt_algorithm, .4 );
  //fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);

  //vector<fastjet::PseudoJet> input_particles;


  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1




  //--------------------------------
  //Initialize helper
  //--------------------------------
  cout << "Pre Razor Helper" << endl;
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData);
  cout << "Post Razor Helper" << endl;


  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  //float maxEvents = 100000;
  
  //variables for checking acceptance on cut-by-cut basis
  float past_HLT_IsoMu24 = 0;
  float greater_than_two_muons = 0;
  float greater_than_one_tag = 0;
  float two_satsifying_probe = 0;
  float leading_muon_pt = 0;
  float subleading_muon_pt = 0;
  float Zmumu_events_passed = 0;
  float total_events_passed = 0;
  float events_timed_noforward = 0;
  float events_timed_matched = 0;
  float events_noforward_matched = 0;
  float events_with_in_time_cluster = 0;
  float events_with_no_forward_hits = 0;
  float events_with_no_forward_hits_DT = 0;
  float events_with_cluster_matched_deltaR = 0;
  float events_with_cluster_matched_deltaR_DT = 0;
  

  //event-by-event variables for probe muon criteria
  float pass_pt = 0;
  float pass_eta = 0;
  float pass_looseId = 0;
  float pass_looseIso = 0;
  float pass_pt_iso = 0;
  float pass_id_iso = 0;
  float pass_id_pt = 0;

  float pass_pt_second = 0;
  float pass_eta_second = 0;
  float pass_looseId_second = 0;
  float pass_looseIso_second = 0;
  float num_two_probe = 0;

int events_with_dt=0;
  float total_muons = 0;
  float total_second_muons = 0;

  float events_two_tag = 0;

  bool checkProbeMuons = false; //this only counts the number of probe muons and skips rest of code!

  float maxEvents = fChain->GetEntries();
  //float maxEvents = 20000;
  cout<<maxEvents<<endl;
  for (Long64_t jentry=0; jentry<maxEvents; jentry++) {
    //begin event
    if(jentry % 1000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    // if (jentry<4000)continue;
    // cout<<jentry<<endl;
    //cout<<"About to load tree"<<endl;
    Long64_t ientry = LoadTree(jentry);
    //cout<<"Loaded tree"<<endl;
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //cout<<"Got entry"<<endl;
    //cout<<"jentry: "<<jentry<<endl;
    if (!HLT_IsoMu24) continue;
    past_HLT_IsoMu24++;
    //fill normalization histogram
    MuonSystem->InitVariables();
    //cout<<"Inited variables"<<endl;
    // std::cout << "deb1 " << jentry << std::endl;

    bool pass_pt_event = false;
    bool pass_eta_event = false;
    bool pass_looseId_event = false;
    bool pass_looseIso_event = false;

    /* Commented out by Alex for now
    if (!isData && signalScan)
    {
      string mh_substring = lheComments->substr(lheComments->find("MH-")+3);
      int mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-")+3);
      int mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-")+6);
      int ctau = stoi(ctau_substring.substr(0,ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

      pair<int,int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0){ //create file and tree
        //format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] =  MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1,0.5,1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1,0.5,1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1,0.5,1.5);



        cout << "Created new output file " << thisFileName << endl;
      }
          //Fill NEvents hist
      NEvents2D[signalPair]->Fill(1.0, genWeight);


    }
    */
    //event info
    //MuonSystem->weight=weight;
    // commented out for MC Simulation Study
    /*
    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
    }
    */
    
    MuonSystem->weight = genWeight;
    NEvents->Fill(1, genWeight);
    //MuonSystem->runNum = run;
    //MuonSystem->lumiSec = luminosityBlock;
    //MuonSystem->evtNum = event;

    //MuonSystem->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
    //MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;

  MuonSystem->runNum = run;
    MuonSystem->lumiSec = luminosityBlock;
    MuonSystem->evtNum = event;

    //check event flags 
    //old flags commented out
    /*
    MuonSystem->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
    MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
    MuonSystem->Flag2_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
    MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
    MuonSystem->Flag_BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
    MuonSystem->Flag_hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;
    MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;
    MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;

    MuonSystem->Flag_all = (Flag_HBHENoiseFilter && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter);
    */
   MuonSystem->Flag_goodVertices = Flag_goodVertices;
      MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
      MuonSystem->Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
      MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
      MuonSystem->Flag_BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
      MuonSystem->Flag_hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;
      MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;
      MuonSystem->Flag_all = Flag_eeBadScFilter && Flag_hfNoisyHitsFilter && Flag_BadPFMuonDzFilter && Flag_BadPFMuonFilter && Flag_EcalDeadCellTriggerPrimitiveFilter
                              && Flag_globalSuperTightHalo2016Filter && Flag_goodVertices;
      if (analysisTag == "Summer24") MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;

      // Flag_ecalBadCalibFilter for nanoAOD: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#ECal_BadCalibration_Filter_Flag
      if (analysisTag == "Summer24") MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
      else{
        MuonSystem->Flag_ecalBadCalibFilter = true;
        if (isData && run >= 362433 && run<=367144)
        {
          if (PuppiMET_pt > 100)
          {
            for(int i = 0; i < nJet; i++)
            {
              if (Jet_pt[i]<50) continue;
              if (!(Jet_eta[i] <= -0.1 && Jet_eta[i]>=-0.5 && Jet_phi[i] <-1.8 && Jet_phi[i]> -2.1)) continue;
              if (!(Jet_neEmEF[i] >0.9 || Jet_chEmEF[i]>0.9)) continue;
              if (deltaPhi(PuppiMET_phi, Jet_phi[i])<2.9) continue;
              Flag_ecalBadCalibFilter = false;
            }
          }
        }
      }
     // jet veto map, following selections here: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
      MuonSystem->jetVeto = true;
      
      for(int i = 0; i < nJet; i++)
      {
        if (Jet_pt[i] <= 15) continue;
        if (Jet_neEmEF[i] + Jet_chEmEF[i] >= 0.9) continue;
        std::bitset<sizeof(int) * 8> jetID(Jet_jetId[i]);
        //cout<<"Filter bit in base 10: "<<binaryNumber<<endl;
        if (!jetID.test(1)) continue;
        //if (!jetPassIDTight[i]) continue;
        //remove overlaps
        bool overlap = false;
        for(int j = 0; j < nMuon; j++)
        {
          if (!Muon_isPFcand[j])continue;
          if (RazorAnalyzer_trigEff::deltaR(Jet_eta[i],Jet_phi[i],Muon_eta[j], Muon_phi[j]) < 0.2) overlap = true;
        }
        if (overlap) continue;
        //cout<<"about to call jet veto map"<<endl;
        helper->getJetVetoMap(0,1);
        if (helper->getJetVetoMap(Jet_eta[i],Jet_phi[i])>0.0) MuonSystem->jetVeto = false;
        if (analysisTag == "Summer24" && helper->getJetVetoFpixMap(Jet_eta[i],Jet_phi[i])>0.0) MuonSystem->jetVeto = false;
      } 
      
    if (!isData){
      MuonSystem->npu = Pileup_nTrueInt;
      MuonSystem->pileupWeight = helper->getPileupWeight(Pileup_nTrueInt);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(Pileup_nTrueInt) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(Pileup_nTrueInt) / MuonSystem->pileupWeight;
    }
    
    
    /* commented out, no gLLPs in nTuples
    if (!isData)
    {
       for(int i = 0; i < 2;i++)
       {
	       MuonSystem->gLLP_eta[MuonSystem->nGLLP] = gLLP_eta[i];
         MuonSystem->gLLP_phi[MuonSystem->nGLLP] = gLLP_phi[i];
         MuonSystem->gLLP_e[MuonSystem->nGLLP] = gLLP_e[i];
         MuonSystem->gLLP_pt[MuonSystem->nGLLP] = gLLP_pt[i];

         MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
         MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = gLLP_decay_vertex_x[i];
         MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = gLLP_decay_vertex_y[i];
         MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = gLLP_decay_vertex_z[i];
         float beta = gLLP_beta[i];
         float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
         float gamma = 1.0/sqrt(1-beta*beta);
         MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex/(beta * gamma);
         MuonSystem->gLLP_beta[MuonSystem->nGLLP] = gLLP_beta[i];


           if (abs(MuonSystem->gLLP_eta[i]) < 2.4
             && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>400
             && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
           if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
             && MuonSystem->gLLP_decay_vertex_r[i] < 800
              && MuonSystem->gLLP_decay_vertex_r[i] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;

            MuonSystem->nGLLP++;
       }

    
    }//end of isData
    */
      //get NPU
      //MuonSystem->npv = nPV;
      //MuonSystem->rho = fixedGridRhoFastjetAll;
      MuonSystem->met = MET_pt;
      MuonSystem->metPhi = MET_phi;
      MuonSystem->puppiMet = PuppiMET_pt;
      MuonSystem->puppiMetPhi = PuppiMET_phi;
      /*
      if(signalScan && !isData)Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      if(signalScan && !isData)
      {
        accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      }
      else if (!isData)
      {
        accep->Fill(1.0, genWeight*MuonSystem->pileupWeight);

      }
      */

      //Triggers
      //for(int i = 0; i < NTriggersMAX; i++){
        //MuonSystem->HLTDecision[i] = HLTDecision[i];
        //MuonSystem->HLT_CscCluster_Loose = HLT_CscCluster_Loose;
        //MuonSystem->HLT_CscCluster_Medium = HLT_CscCluster_Medium;
        //MuonSystem->HLT_CscCluster_Tight = HLT_CscCluster_Tight;
        //MuonSystem->HLT_IsoMu20 = HLT_IsoMu20;
        //MuonSystem->HLT_IsoMu24 = HLT_IsoMu24;
        //MuonSystem->L1_SingleMuShower_Nominal = L1_SingleMuShower_Nominal;
        //MuonSystem->L1_SingleMuShower_Tight = L1_SingleMuShower_Tight;
      //}

      //*************************************************************************
      //Start Object Selection
      //*************************************************************************

      std::vector<leptons> Leptons;
      //-------------------------------
      //Muons
      //-------------------------------
      
      //cout<<"nMuon "<<nMuon<<endl;
      if (nMuon>=2){
        greater_than_two_muons++;
      }
      for( int i = 0; i < nMuon; i++ )  
      {
        total_muons++;
        bool pass_pt_muon = false;
        bool pass_eta_muon = false;
        bool pass_looseId_muon = false;
        bool pass_looseIso_muon = false;

        bool on_second_muon = pass_looseId_event && pass_pt_event && pass_eta_event && pass_looseIso_event;
        if (on_second_muon){
          total_second_muons++;
        }
        //if(!Muon_looseId[i]) continue;
        if(Muon_looseId[i] && !on_second_muon) {
          pass_looseId_event = true;
          pass_looseId_muon = true;
          pass_looseId++;
        }
        if(Muon_looseId[i] && on_second_muon) {
          pass_looseId_event = true;
          pass_looseId_muon = true;
          pass_looseId++;
          pass_looseId_second++;
        }
        //cout<<"Found not loose muon"<<endl;
        //cout<<"Muon_pt: "<<Muon_pt[i]<<endl;
        //if(Muon_pt[i] < 20.0) continue; //do we want this cut, only hard muons?
        if(Muon_pt[i] > 20.0 && !on_second_muon) {
          pass_pt_event = true;
          pass_pt_muon = true;
          pass_pt++;
        }

        if(Muon_pt[i]>20.0 && on_second_muon) {
          pass_pt_event = true;
          pass_pt_muon = true;
          pass_pt++;
          pass_pt_second++;
        }
        //cout<<"Found hard, not loose muon"<<endl;
        //cout<<"Muon_eta"<<Muon_eta[i]<<endl;
        //if(fabs(Muon_eta[i]) > 2.4) continue;
        if(fabs(Muon_eta[i]) < 2.4 && !on_second_muon){
          pass_eta_event = true;
          pass_eta_muon = true;
          pass_eta++;
        }
        if(fabs(Muon_eta[i])<2.4 && on_second_muon) {
          pass_eta_event = true;
          pass_eta_muon = true;
          pass_eta++;
          pass_eta_second++;
        }
        //remove overlaps
        
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzer_trigEff::deltaR(Muon_eta[i],Muon_phi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;


        leptons tmpMuon;
        tmpMuon.lepton.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i], Muon_phi[i], MU_MASS);
        tmpMuon.pdgId = 13 * -1 * Muon_charge[i];
        tmpMuon.dZ = Muon_dz[i];
        tmpMuon.passId = Muon_tightId[i];
        //ask about iso variables
        //float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / Muon_pt[i];
        float muonIso = Muon_pfRelIso04_all[i]; //added by alex, address set to Muon_pfRelIso04_all branch orignally fron nanoAODs
        tmpMuon.passLooseIso = muonIso<0.25;
        //tmpMuon.passLooseIso = muonIso<0.40; //loosen iso for probe muon
        tmpMuon.passTightIso = muonIso<0.15;
        tmpMuon.passVTightIso = muonIso<0.10;
        tmpMuon.passVVTightIso = muonIso<0.05;

        tmpMuon.passVetoId = false;
        //if (!tmpMuon.passLooseIso) continue;
        if(tmpMuon.passLooseIso && !on_second_muon) {
          pass_looseIso_event = true;
          pass_looseIso_muon = true;
          pass_looseIso++;
        }
        if(tmpMuon.passLooseIso && on_second_muon) {
          pass_looseIso_event = true;
          pass_looseIso_muon = true;
          pass_looseIso++;
          pass_looseIso_second++;
        }

        if(Muon_looseId[i] && Muon_pt[i] > 20.0){
          pass_id_pt++;
        }

        if(Muon_looseId[i] && tmpMuon.passLooseIso){
          pass_id_iso++;
        }

        if(Muon_pt[i] && tmpMuon.passLooseIso){
          pass_pt_iso++;
        }


        if (!(pass_looseId_muon && pass_eta_muon && pass_pt_muon && pass_looseIso_muon)) continue;
        Leptons.push_back(tmpMuon);
      }
      if (Leptons.size() != 2) continue;
      MuonSystem->numProbeMuons = Leptons.size();
      
      if (checkProbeMuons){
         MuonSystem->tree_->Fill();
        continue;
      }
      //-------------------------------
      //Electrons
      //-------------------------------
      
      /*COMMENTED OUT BECAUSE MERGED NTUPLES DO NOT CONTAIN ELECTRONS
      for( int i = 0; i < nElectrons; i++ )
      {
        if (!ele_passCutBasedIDVeto[i]) continue;
        if(elePt[i] < 35) continue;
        if(fabs(eleEta[i]) > 2.5) continue;

        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzer_trigEff::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;
        leptons tmpElectron;
        tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
        tmpElectron.pdgId = 11 * -1 * eleCharge[i];
        tmpElectron.dZ = ele_dZ[i];
        tmpElectron.passId = ele_passCutBasedIDTight[i];
        Leptons.push_back(tmpElectron);
      }
      */
      
      sort(Leptons.begin(), Leptons.end(), my_largest_pt);
      //cout<< "length of leptons: " << Leptons.size() << endl;
      if (Leptons.size() == 2) num_two_probe++;
      for ( auto &tmp : Leptons ) //note that these lepton variables will only include muon info
      {
        MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
        MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
        MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
        MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
        //MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
        //MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
        //MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.passId;
        //MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
        //MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
        //MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
        //MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
        //MuonSystem->nLeptons++;
      }
     //cout<<"nLeptons: "<<MuonSystem->nLeptons<<endl; 
    // the following code , until jets, is for reconsturction of the Z-peak. It is taken from https://github.com/cms-lpc-llp/llp_analyzer/blob/master/analyzers/llp_MuonSystem_TnP.cc#L489-L526
    double ZMass = -999;
    double ZPt = -999;
    double tmpDistToZPole = 9999;
    pair<uint,uint> ZCandidateLeptonIndex;
    bool foundZ = false;
    TLorentzVector ZCandidate;
    double leadingLepPt = 0.0;
    int tagCount=0;
    for( uint i = 0; i < Leptons.size(); i++ )
    {
      for( uint j = i+1; j < Leptons.size(); j++ )
      {
        if (!( Leptons[i].pdgId == -1*Leptons[j].pdgId )) continue;// same flavor opposite charge
        double tmpMass = (Leptons[i].lepton+Leptons[j].lepton).M();
        //select the pair closest to Z pole mass
        if ( fabs( tmpMass - Z_MASS) < tmpDistToZPole)
        {
          tmpDistToZPole = tmpMass;
          if (Leptons[i].pdgId > 0)
          {
            ZCandidateLeptonIndex = pair<int,int>(i,j);
          }
          else
          {
            ZCandidateLeptonIndex = pair<int,int>(j,i);
          }
          ZMass = tmpMass;
          ZPt = (Leptons[i].lepton+Leptons[j].lepton).Pt();
          ZCandidate = Leptons[i].lepton+Leptons[j].lepton;
          leadingLepPt = max(Leptons[i].lepton.Pt(),Leptons[j].lepton.Pt());
          foundZ = true;
        }
      }
    }
    if (ZMass<50||ZMass>120) continue;
    // if (foundZ  && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
    if (foundZ  && Leptons.size() == 2 )
    {
      MuonSystem->ZMass = ZMass;
      //MuonSystem->ZPt   = ZPt;
      //MuonSystem->ZEta  = ZCandidate.Eta();
      //MuonSystem->ZPhi  = ZCandidate.Phi();
      //MuonSystem->ZleptonIndex1 = ZCandidateLeptonIndex.first;
      //MuonSystem->ZleptonIndex2 = ZCandidateLeptonIndex.second;
      //MuonSystem->category = 2;
    } // endif foundZ
    /*else{
      for ( unsigned int i = Leptons.size(); i>0; --i )
      {
        int index = i-1;
        if (abs(Leptons[index].pdgId) == 13 && Leptons[index].lepton.Pt() < wh_Muon_pt_cut)
        {
           Leptons.erase(Leptons.begin() + index);
        }
        else if (abs(Leptons[index].pdgId) == 11 && Leptons[index].lepton.Pt() < wh_elePt_cut)
        {
          Leptons.erase(Leptons.begin() + index);
        }
      }
      if (Leptons.size() == 1) MuonSystem->category = 1;
      else MuonSystem->category = 0;
    }*/
    bool tag = false;
    int numTag=0;
    std::vector<float> deltaRs;
    std::vector<float> trigObjIndex;
    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      
      if (!isData)
      {
        /* Commenting out all razor helper (corrections specific to Run 2?)
        MuonSystem->lepTriggerSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
        MuonSystem->lepTightIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
        MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
        MuonSystem->lepTightIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
        MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
        MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
        MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
        MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
        MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
        MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
        */
      }


      
      if (MuonSystem->lepPassTightIso[MuonSystem->nLeptons] && MuonSystem->lepPassId[MuonSystem->nLeptons] && MuonSystem->lepPt[MuonSystem->nLeptons]>26)
      //if (MuonSystem->lepPassTightIso[MuonSystem->nLeptons] && MuonSystem->lepPassId[MuonSystem->nLeptons])
      {
        //comment out next three lines when matching to trigger object
        //MuonSystem->lepTag[MuonSystem->nLeptons] = true;
        //tag = true;
        //tagCount++;
        //now, attempt to match trigger object to tag muon. If not matched in deltaR, we set tag to false
        //cout<<"Found tag muon"<<endl;
        
        bool matchedTriggerObj = false;
        
        for (int trigObjNum=0; trigObjNum < nTrigObj; trigObjNum++)
        {

          //cout<<"In trigger object loop"<<endl;
          if (abs(TrigObj_id[trigObjNum]) != 13) continue;
          //out<<"Filter bit: "<<TrigObj_filterBits[trigObjNum]<<endl;
          //now convert filterBits into base 10 and see whether the bit for 2^3 is set
          std::bitset<sizeof(int) * 8> binaryNumber(TrigObj_filterBits[trigObjNum]);
          //cout<<"Filter bit in base 10: "<<binaryNumber<<endl;
          if (!(binaryNumber.test(3) && binaryNumber.test(0))) continue;
          //cout<<"Trigger muon with filter bit found"<<endl;
          if (deltaR(tmp.lepton.Eta(), tmp.lepton.Phi(), TrigObj_eta[trigObjNum], TrigObj_phi[trigObjNum]) < 0.1)
          {
            matchedTriggerObj = true;
            deltaRs.push_back(deltaR(tmp.lepton.Eta(), tmp.lepton.Phi(), TrigObj_eta[trigObjNum], TrigObj_phi[trigObjNum]));
            trigObjIndex.push_back(trigObjNum);
            //cout<<"Matched trigger object to tag muon"<<endl;
            break;
          }
        }
      
    
        if (matchedTriggerObj){
          numTag++;
          MuonSystem->lepTag[MuonSystem->nLeptons] = true;
          tag = true;
          tagCount++;
        }
        else{
          MuonSystem->lepTag[MuonSystem->nLeptons] = false;
          /*commenting out scale factor stuff, will compute in python
          if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
        }
        */
        }
      }
        MuonSystem->nLeptons++;
    }
        /*
        if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepTightIdSF[MuonSystem->nLeptons] * MuonSystem->lepTightIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
        
      }
      else
      {
        MuonSystem->lepTag[MuonSystem->nLeptons] = false;
        if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
      }
      */


      /*if(!isData) //variables names modified to match new makeClass, not sure about mother one, before was gParticleMotherId
      {
        for (int i=0; i < nGenPart; i++)
        { float tmpDR = deltaR(GenPart_eta[i],GenPart_phi[i],MuonSystem->lepEta[MuonSystem->nLeptons],MuonSystem->lepPhi[MuonSystem->nLeptons]);
          if ((abs(GenPart_pdgId[i]) == 13) && abs(GenPart_genPartIdxMother[i]) == 23 && tmpDR<0.4) MuonSystem->lepFromZ[MuonSystem->nLeptons] = true;
        }
      }
    
      MuonSystem->nLeptons++;
    }*/
    if (numTag==2){
      /*
	    cout<<"numTag: "<<numTag<<endl;
      for (int i=0; i<deltaRs.size(); i++){
        cout<<"deltaR: "<<deltaRs[i]<<endl;
      }
      for (int i=0; i<trigObjIndex.size(); i++){
        cout<<"trigObjIndex: "<<trigObjIndex[i]<<endl;
      }
      */
      events_two_tag++;
    }

    //compute acceptances
    //if(MuonSystem->nLeptons==2)two_satsifying_probe++;
    //if (abs(MuonSystem->lepPdgId[0])!=13)continue; not relevant, all muons
    //if (abs(MuonSystem->ZMass)<50)continue; removed to see whole Z-peak
    //if (abs(MuonSystem->lepPt[0])>35)leading_muon_pt++;
    //if (abs(MuonSystem->lepPt[1])>20)subleading_muon_pt++;
    // if(MuonSystem->ZMass<120)continue; remove to see whole Z-peak
    //if (tagCount>0)greater_than_one_tag++;


    if(MuonSystem->nLeptons!=2)continue;
    //if(MuonSystem->category!=2)continue;
    if (abs(MuonSystem->lepPdgId[0])!=13)continue;
    //if (abs(MuonSystem->ZMass)<50)continue; removed to see whole Z-peak
    //if (abs(MuonSystem->lepPt[0])<35)continue;
    if (abs(MuonSystem->lepPt[1])<20)continue;
    // if(MuonSystem->ZMass<120)continue; remove to see whole Z-peak
    if (tagCount==0) continue;
    //cout<<"at lep tag"<<endl;
    //if(MuonSystem->lepTag[0] == MuonSystem->lepTag[1]) continue;
    //cout<<"made it past lepTag cuts" << endl;
      //require one tag one probe
    
    /*
    if(!isData)
    {
      MuonSystem->lepOverallSF = 1.0 - (1.0 - MuonSystem->lepSF[0] * MuonSystem->lepEff[0]) * (1 - MuonSystem->lepSF[1] * MuonSystem->lepEff[1]);
      MuonSystem->lepOverallSF = MuonSystem->lepOverallSF / (1.0 - (1.0 - MuonSystem->lepEff[0]) * (1 - MuonSystem->lepEff[1]));
    }
    else{
      MuonSystem->lepOverallSF = 1.0;
    }
    */
    //cout<<"made it past muon cuts!" << endl;
  MuonSystem->numTag = tagCount;
  TLorentzVector met;
  met.SetPtEtaPhiE(MET_pt,0,MET_phi,MET_pt);
  if ( Leptons.size() > 0 )
  {
    TLorentzVector visible = Leptons[0].lepton;
    //MuonSystem->MT = GetMT(visible,met);
  }
   Zmumu_events_passed++;   
    //cout<<"Zmumu events passed: "<<Zmumu_events_passed<<endl;
    for (int i=0; i <tagCount; i++){
      if (isData){ZMass_Hist->Fill(MuonSystem->ZMass);}
      else{ZMass_Hist->Fill(MuonSystem->ZMass,genWeight*MuonSystem->pileupWeight);}

    }
      MuonSystem->tree_->Fill();

      }
      if(!isData && signalScan)
      {
        for(auto &filePtr : Files2D)
         {
           cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
           filePtr.second->cd();
           Trees2D[filePtr.first]->Write();
           NEvents2D[filePtr.first]->Write("NEvents");
           Total2D[filePtr.first]->Write("Total");
           accep2D[filePtr.first]->Write("acceptance");
           accep_met2D[filePtr.first]->Write("acceptance_met");
           filePtr.second->Close();

         }
      }
      else if (!isData)
      {
         cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
         cout << "Writing output trees..." << endl;
         outFile->cd();
         MuonSystem->tree_->Write();
         NEvents->Write();
         accep->Write("acceptance");
         accep_met->Write("acceptance_met");
         outFile->Close();
      }


      else
      {
        cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
        cout << "Writing output trees..." << endl;
        outFile->cd();
        MuonSystem->tree_->Write();
        Nmet200->Write();
        NmetFilter->Write();
        Nlep0->Write();
        Njet1->Write();
        NcosmicVeto->Write();
        NEvents->Write();
        probe_cluster_minDeltaR->Write();
        // outFile->Write();
        outFile->Close();
      }
      /*
      cout<<"Events Two Tagged: "<<events_two_tag<<endl;
      cout<<"Total Events Passed: "<<total_events_passed<<endl;
      //cout<<"Acceptances: "<<endl;
      cout<<"Cutflow (only two-muon events): "<<endl;
      cout<<"Pass HLT_IsoMu24 filter: "<<past_HLT_IsoMu24/maxEvents*100<<"%"<<endl;
      cout<<"At least two muons pre-filtering: "<<greater_than_two_muons/past_HLT_IsoMu24*100<<"%"<<endl;
      //cout<<"Probe Muon Cut Efficiencies (at least one muon passing criteria): "<<endl;
      cout<<"Probe Muon Selection Cuts (out of total number of muons): "<<endl;
      //cout<<"Pt > 20: "<<pass_pt/total_muons*100<<"%"<<endl;
      cout<<"Eta < 2.4: "<<pass_eta/total_muons*100<<"%"<<endl;
      cout<<"Pass Loose Id: "<<pass_looseId/total_muons*100<<"%"<<endl;
      cout<<"Pass Loose Iso: "<<pass_looseIso/total_muons*100<<"%"<<endl;
      cout<<"Pass Loose Id and Iso: "<<pass_id_iso/total_muons*100<<"%"<<endl;
      cout<<"Pass Loose Id and pT: "<<pass_id_pt/total_muons*100<<"%"<<endl;
      cout<<"Pass pT and Iso: "<<pass_pt_iso/total_muons*100<<"%"<<endl;
      cout<<"Probe Muon Selection Cuts (provided at least 1 muon already passed): "<<endl;
      cout<<"Pt > 20: "<<pass_pt_second/total_second_muons*100<<"%"<<endl;
      cout<<"Eta < 2.4: "<<pass_eta_second/total_second_muons*100<<"%"<<endl;
      cout<<"Pass Loose Id: "<<pass_looseId_second/total_second_muons*100<<"%"<<endl;
      cout<<"Pass Loose Iso: "<<pass_looseIso_second/total_second_muons*100<<"%"<<endl;
      cout<<"Exactly Two Muons Passing Probe: "<<num_two_probe/past_HLT_IsoMu24*100<<"%"<<endl;
      cout<<"Two Muons Passing Probe with Opposite Charges: "<<two_satsifying_probe/past_HLT_IsoMu24*100<<"%"<<endl;
      cout<<"At Least One Tag: "<<greater_than_one_tag/past_HLT_IsoMu24*100<<"%"<<endl;
      cout<<"Leading Muon pT>35: "<<leading_muon_pt/past_HLT_IsoMu24*100<<"%"<<endl;
      cout<<"Subleading Muon pT>20: "<<subleading_muon_pt/past_HLT_IsoMu24*100<<"%"<<endl;
      cout<<"Cumulative Efficiency - finding at least 1 TnP Pair: "<<Zmumu_events_passed/past_HLT_IsoMu24*100<<"%"<<endl;
      cout<<"Efficiency for Matching to Cluster (no ME11/12, in time, deltaR) for events with tag-probe pair"<<endl;
      cout<<"At Least 1 Cluster with no hits in ME11/12: "<<float(events_with_no_forward_hits)/Zmumu_events_passed*100<<"%"<<endl;
      cout<<"At Least 1 In-Time Cluster: "<<float(events_with_in_time_cluster)/Zmumu_events_passed*100<<"%"<<endl;
      cout<<"At Least 1 Cluster With DeltaR(Cluster, Probe)<0.4: "<<float(events_with_cluster_matched_deltaR)/Zmumu_events_passed*100<<"%"<<endl;
      cout<<"At Least 1 In-Time Cluster with no hits in ME11/12: "<<float(events_timed_noforward)/Zmumu_events_passed*100<<"%"<<endl;
      cout<<"At Least 1 In-Time Cluster With DeltaR(Cluster, Probe)<0.4: "<<float(events_timed_matched)/Zmumu_events_passed*100<<"%"<<endl;
      cout<<"At Least 1 Cluster With no hits in ME11/12 and DeltaR(Cluster, Probe)<0.4: "<<float(events_noforward_matched)/Zmumu_events_passed*100<<"%"<<endl;
      cout<<"Cumulative Matching Rate : "<<total_events_passed/Zmumu_events_passed*100<<"%"<<endl;
      cout<<"Events with dt : "<<events_with_dt<<endl;
      cout<<"Total Efficiency: "<<total_events_passed/maxEvents*100<<"%"<<endl;
      */

  
}

