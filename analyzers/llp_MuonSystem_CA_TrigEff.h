#ifndef DEF_llp_MuonSystem_CA_TrigEff
#define DEF_llp_MuonSystem_CA_TrigEff

#include "RazorAnalyzer_trigEff.h"

class llp_MuonSystem_CA_TrigEff: public RazorAnalyzer_trigEff {
    public: 
        llp_MuonSystem_CA_TrigEff(TTree *tree=0): RazorAnalyzer_trigEff(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
