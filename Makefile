include Makefile.inc

DIRS = python
SRCDIR = src
INCLUDEDIR = include
ANADIR = analyzers
BINDIR = bin
INCLUDELIST= SimpleTable.h Linkdef.h

ANALYZERS = $(wildcard $(ANADIR)/*.cc)
ANALYZERSH = $(ANALYZERS:cc=h)
ANALYZERSOBJ = $(ANALYZERS:cc=o)
RUNNERS = $(addprefix $(BINDIR)/Run,$(notdir $(basename $(ANALYZERS))))
RUNNERSCC = $(addsuffix .cc,$(addprefix $(ANADIR)/,$(notdir $(RUNNERS))))
UTILS =$(SRCDIR)/RazorHelper.cc $(SRCDIR)/DBSCAN.cc  $(SRCDIR)/CACluster.cc ${SRCDIR}/TreeMuonSystemCombination.cc ${SRCDIR}/TreeMuonSystemCombination_TnP.cc ${SRCDIR}/TreeMuonSystem_Skim_Merge_TnP.cc ${SRCDIR}/TreeMuonSystem_VetoEff.cc
UTILSOBJ = $(UTILS:cc=o)
EXECUTABLES = NormalizeNtuple SkimNtuple $(RUNNERS)
#EXECUTABLES = $(RUNNERS)
HELPERSCRIPT = python/MakeAnalyzerCode.py
HELPERSCRIPT_TnP = python/MakeAnalyzerCode_trigEff.py
HELPERSCRIPT_trigEff = python/MakeAnalyzerCode_trigEff.py
HELPERSCRIPT_merged = python/MakeAnalyzerCodeMerged.py


.PHONY: clean all lxplus copy_runners copy_runners_TnP copy_runners_TrigEff copy_runners_TnP_noClusters copy_runners_TnP_mdsnano

all: copy_runners copy_runners_TnP copy_runners_TnP_noClusters copy_runners_TnP_noClusters_VetoEff copy_runners_TnP_noClusters_VetoEff_mdsnano copy_runners_TnP_mdsnano copy_runners_TrigEff copy_runners_TrigEff_mdsnano $(EXECUTABLES)

lxplus: all

print_vars:
	echo $(ANALYZERS)

clean:
	@-rm $(EXECUTABLES)
	@rm -f $(SRCDIR)/*.o $(ANADIR)/*.o

#copy_runners:
#		@for d in $(subst Run,,$(notdir $(basename $(RUNNERSCC)))); do ( if [ ! -f "src/Run"$$d".cc" ]; then echo $$d" file does not exists, copying"; $(HELPERSCRIPT) $$d; fi ) ; done


copy_runners:
		@for d in $(subst Run,,$(notdir $(basename $(RUNNERSCC)))); do \
			if [ $$d != "llp_MuonSystem_CA_TnP" ] && [ $$d != "llp_MuonSystem_CA_TrigEff" ] && [ $$d != "llp_MuonSystem_CA_TrigEff_mdsnano" ]&& [ $$d != "llp_MuonSystem_CA_TnP_noClusters" ] && [ $$d != "llp_MuonSystem_CA_TnP_mdsnano" ] && [ $$d != "llp_MuonSystem_CA_TnP_noClusters_VetoEff" ] && [ $$d != "llp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano" ]; then \
				if [ ! -f "src/Run"$$d".cc" ]; then \
					echo $$d "file does not exist, copying"; \
					echo "Running python/MakeAnalyzerCode.py $$d"; \
					$(HELPERSCRIPT) $$d; \
				fi; \
			fi; \
		done

copy_runners_TnP:		
		d="llp_MuonSystem_CA_TnP"; \
		if [ ! -f "src/Run"$$d".cc" ]; then \
		    echo $$d" file does not exists, copying"; \
			echo "Running python/MakeAnalyzerCode_trigEff.py $$d"; \
			$(HELPERSCRIPT_TnP) $$d; \
		fi

copy_runners_TnP_noClusters:		
		d="llp_MuonSystem_CA_TnP_noClusters"; \
		if [ ! -f "src/Run"$$d".cc" ]; then \
		    echo $$d" file does not exists, copying"; \
			echo "Running python/MakeAnalyzerCode_trigEff.py $$d"; \
			$(HELPERSCRIPT_TnP) $$d; \
		fi

copy_runners_TnP_noClusters_VetoEff:		
		d="llp_MuonSystem_CA_TnP_noClusters_VetoEff"; \
		if [ ! -f "src/Run"$$d".cc" ]; then \
		    echo $$d" file does not exists, copying"; \
			echo "Running python/MakeAnalyzerCode_trigEff.py $$d"; \
			$(HELPERSCRIPT_TnP) $$d; \
		fi

copy_runners_TnP_noClusters_VetoEff_mdsnano:		
		d="llp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano"; \
		if [ ! -f "src/Run"$$d".cc" ]; then \
		    echo $$d" file does not exists, copying"; \
			echo "Running python/MakeAnalyzerCode_merged.py $$d"; \
			$(HELPERSCRIPT_merged) $$d; \
		fi

copy_runners_TnP_mdsnano:		
		d="llp_MuonSystem_CA_TnP_mdsnano"; \
		if [ ! -f "src/Run"$$d".cc" ]; then \
		    echo $$d" file does not exists, copying"; \
			echo "Running python/MakeAnalyzerCode_merged.py $$d"; \
			$(HELPERSCRIPT_merged) $$d; \
		fi

copy_runners_TrigEff:		
		d="llp_MuonSystem_CA_TrigEff"; \
		if [ ! -f "src/Run"$$d".cc" ]; then \
		    echo $$d" file does not exists, copying"; \
				echo "Running python/MakeAnalyzerCode_trigEff.py $$d"; \
		    $(HELPERSCRIPT_trigEff) $$d; \
		fi

copy_runners_TrigEff_mdsnano:		
		d="llp_MuonSystem_CA_TrigEff_mdsnano"; \
		if [ ! -f "src/Run"$$d".cc" ]; then \
		    echo $$d" file does not exists, copying"; \
				echo "Running python/MakeAnalyzerCode_merged.py $$d"; \
		    $(HELPERSCRIPT_merged) $$d; \
		fi

$(INCLUDEDIR)/rootdict.cc:
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) -I$(INCLUDEDIR) $(INCLUDELIST)

$(SRCDIR)/SimpleTable.o: $(SRCDIR)/SimpleTable.cc
	$(CXX) -c $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/llp_event.o: $(SRCDIR)/llp_event.C $(INCLUDEDIR)/llp_event.h
	$(CXX) $(SRCDIR)/llp_event.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/llp_event.o $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

#added by me
$(SRCDIR)/llp_event_TnP.o: $(SRCDIR)/llp_event_TnP.C $(INCLUDEDIR)/llp_event_TnP.h
	$(CXX) $(SRCDIR)/llp_event_TnP.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

#added by me
$(SRCDIR)/RazorAnalyzer_TnP.o: $(SRCDIR)/llp_event_TnP.o $(SRCDIR)/RazorAnalyzer_TnP.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer_TnP.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

#added by me
$(SRCDIR)/llp_event_trigEff.o: $(SRCDIR)/llp_event_trigEff.C $(INCLUDEDIR)/llp_event_trigEff.h
	$(CXX) $(SRCDIR)/llp_event_trigEff.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/merged_event.o: $(SRCDIR)/merged_event.C $(INCLUDEDIR)/merged_event.h
	$(CXX) $(SRCDIR)/merged_event.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

#added by me
$(SRCDIR)/RazorAnalyzer_trigEff.o: $(SRCDIR)/llp_event_trigEff.o  $(SRCDIR)/RazorAnalyzer_trigEff.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer_trigEff.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorAnalyzerMerged.o: $(SRCDIR)/merged_event.o  $(SRCDIR)/RazorAnalyzerMerged.cc
	$(CXX) $(SRCDIR)/RazorAnalyzerMerged.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/llp_event_skim_merge.o: $(SRCDIR)/llp_event_skim_merge.C $(INCLUDEDIR)/llp_event_skim_merge.h
	$(CXX) $(SRCDIR)/llp_event_skim_merge.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorAnalyzer_Skim_Merge.o: $(SRCDIR)/llp_event_skim_merge.o $(SRCDIR)/RazorAnalyzer_Skim_Merge.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer_Skim_Merge.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(UTILSOBJ): %.o: %.cc
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(ANALYZERSOBJ): $(ANADIR)/%.o: $(ANADIR)/%.cc $(ANADIR)/%.h
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(BINDIR)/Runllp_MuonSystem: $(SRCDIR)/llp_event.o $(SRCDIR)/RazorAnalyzer.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem.o $(SRCDIR)/Runllp_MuonSystem.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(BINDIR)/Runllp_MuonSystem_CA: $(SRCDIR)/llp_event.o $(SRCDIR)/RazorAnalyzer.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA.o $(SRCDIR)/Runllp_MuonSystem_CA.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(BINDIR)/Runllp_MuonSystem_CA_TnP: $(SRCDIR)/llp_event_trigEff.o $(SRCDIR)/merged_event.o $(SRCDIR)/RazorAnalyzer_trigEff.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_TnP.o $(SRCDIR)/Runllp_MuonSystem_CA_TnP.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(BINDIR)/Runllp_MuonSystem_CA_TnP_noClusters: $(SRCDIR)/llp_event_trigEff.o $(SRCDIR)/merged_event.o $(SRCDIR)/RazorAnalyzer_trigEff.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_TnP_noClusters.o $(SRCDIR)/Runllp_MuonSystem_CA_TnP_noClusters.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(BINDIR)/Runllp_MuonSystem_CA_TnP_noClusters_VetoEff: $(SRCDIR)/llp_event_trigEff.o $(SRCDIR)/merged_event.o $(SRCDIR)/RazorAnalyzer_trigEff.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_TnP_noClusters_VetoEff.o $(SRCDIR)/Runllp_MuonSystem_CA_TnP_noClusters_VetoEff.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(BINDIR)/Runllp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano: $(SRCDIR)/merged_event.o $(SRCDIR)/RazorAnalyzerMerged.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano.o $(SRCDIR)/Runllp_MuonSystem_CA_TnP_noClusters_VetoEff_mdsnano.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
		
$(BINDIR)/Runllp_MuonSystem_CA_TnP_mdsnano: $(SRCDIR)/merged_event.o $(SRCDIR)/RazorAnalyzerMerged.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_TnP_mdsnano.o $(SRCDIR)/Runllp_MuonSystem_CA_TnP_mdsnano.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(BINDIR)/Runllp_MuonSystem_CA_TrigEff: $(SRCDIR)/llp_event_trigEff.o $(SRCDIR)/merged_event.o $(SRCDIR)/RazorAnalyzer_trigEff.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_TrigEff.o $(SRCDIR)/Runllp_MuonSystem_CA_TrigEff.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(BINDIR)/Runllp_MuonSystem_CA_TrigEff_mdsnano: $(SRCDIR)/merged_event.o $(SRCDIR)/RazorAnalyzerMerged.o $(UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_TrigEff_mdsnano.o $(SRCDIR)/Runllp_MuonSystem_CA_TrigEff_mdsnano.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

NormalizeNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/NormalizeNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

SkimNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/SkimNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
