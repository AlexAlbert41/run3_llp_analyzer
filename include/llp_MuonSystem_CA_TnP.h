#ifndef DEF_llp_MuonSystem_CA_TnP
#define DEF_llp_MuonSystem_CA_TnP

#include "RazorAnalyzer_trigEff.h"

class llp_MuonSystem_CA_TnP: public RazorAnalyzer_trigEff {
    public: 
        llp_MuonSystem_CA_TnP(TTree *tree=0): RazorAnalyzer_trigEff(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
