# Instructions for Making HMT Trigger Efficiency Plots

Follow the Setup instructions in the main README on the cmslpc. The OS doesn't matter, since the code will be run on an EL7 image instead.

---

## Environment Setup

On a separate terminal, re-log in to LPC, but forward port 8888 to the remote machine:


`ssh -L localhost:8888:localhost:8888 <username>@cmslpc-el8.fnal.gov`

Set up an EL7 image:
`source /cvmfs/cms.cern.ch/cmsset_default.sh`
`cmssw-el7 --bind /uscms_data/d1/<username> # will bind to your home directory`

Within the image:
`source /cvmfs/cms.cern.ch/cmsset_default.sh`
### Source an LHC computing grid environment with which the notebooks can run
`source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh` <br>
### Launch a Jupyter notebook
`jupyter notebook --no-browser --port=8888 --ip 127.0.0.1`

The printed URLs can be copied into your browser, and notebooks can be opened.
Plotting Scripts
I've written code that can produce two types of plots that show the efficiency for a given file.
1. Coffea Hist Plots <br>
The first are plots that use the now-legacy coffea.hist package. This hist style was used for our DP note: DP2024_099.pdf. Since the code takes a very long time to run on large datasets, I think the best strategy is to store the histograms via pickle, and then make the final plots using a separate script.
The notebooks are in the `Plotting_Scripts directory`: <br>
To process the data / make the pickle: `L1_TrigEff_082726_Cleaned.ipynb` <br>
To make the coffea hists from the pickle, you can use: `Plot_Effs_From_Pkl.ipynb` <br>
2. ROOT Histograms <br>
The other set of plots are ROOT histograms that we use directly for our analysis. These are probably more important. For this, you can run the notebook: <br>
`Make_Efficiency_ROOT_Files.ipynb` <br>
Output Files <br>
The output files I got for 2024 are in the directory: <br>
`2024_trigger_efficiencies_latertest` <br>
