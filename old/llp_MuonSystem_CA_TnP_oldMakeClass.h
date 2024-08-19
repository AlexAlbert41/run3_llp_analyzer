#ifndef DEF_llp_MuonSystem_CA_TnP_oldMakeClass
#define DEF_llp_MuonSystem_CA_TnP_oldMakeClass

#include "RazorAnalyzer.h"

class llp_MuonSystem_CA_TnP_oldMakeClass: public RazorAnalyzer {
    public: 
        llp_MuonSystem_CA_TnP_oldMakeClass(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
