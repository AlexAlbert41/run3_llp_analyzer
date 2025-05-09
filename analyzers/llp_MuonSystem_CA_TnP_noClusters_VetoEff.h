#ifndef DEF_llp_MuonSystem_CA_TnP_noClusters_VetoEff
#define DEF_llp_MuonSystem_CA_TnP_noClusters_VetoEff

#include "RazorAnalyzer_trigEff.h"

class llp_MuonSystem_CA_TnP_noClusters_VetoEff: public RazorAnalyzer_trigEff {
    public: 
        llp_MuonSystem_CA_TnP_noClusters_VetoEff(TTree *tree=0): RazorAnalyzer_trigEff(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
