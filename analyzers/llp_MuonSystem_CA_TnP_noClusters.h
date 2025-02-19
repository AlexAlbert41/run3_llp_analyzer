#ifndef DEF_llp_MuonSystem_CA_TnP_noClusters
#define DEF_llp_MuonSystem_CA_TnP_noClusters

#include "RazorAnalyzer_trigEff.h"

class llp_MuonSystem_CA_TnP_noClusters: public RazorAnalyzer_trigEff {
    public: 
        llp_MuonSystem_CA_TnP_noClusters(TTree *tree=0): RazorAnalyzer_trigEff(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
