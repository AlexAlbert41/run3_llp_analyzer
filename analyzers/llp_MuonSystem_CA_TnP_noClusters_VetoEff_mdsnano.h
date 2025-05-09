#ifndef DEF_llp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano
#define DEF_llp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano

#include "RazorAnalyzerMerged.h"

class llp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano: public RazorAnalyzerMerged {
    public: 
        llp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano(TTree *tree=0): RazorAnalyzerMerged(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
