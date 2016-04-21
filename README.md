mc
==============

Code for private production of FullSim MC points.

#### Generate gridpack
In lxplus, check out [genproductions](https://github.com/cms-sw/genproductions) 

    git clone https://github.com/cms-sw/genproductions
    cd genproductions/bin/MadGraph5_aMCatNLO/cards/production/13TeV/
    mkdir TChiHH_mChi-400
    
Place the `cards/TChiHH_mChi-400` folder into `genproductions/bin/MadGraph5_aMCatNLO/cards/production/13TeV/`, and send a job to the lxplus batch as instructed in the [aMCatNLO twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuickGuideMadGraph5aMCatNLO#Quick_tutorial_on_how_to_produce):

    kinit -Af manuelf@CERN.CH
    cd ~/work/code/genproductions/bin/MadGraph5_aMCatNLO
    ./submit_gridpack_generation.sh 12000 12000 1nw TChiHH_mChi-400 cards/production/13TeV/TChiHH_mChi-400/ 8nh

This generates a gridback in the file `TChiHH_mChi-400_tarball.tar.xz`.


#### Generate MINIAOD from the gridpack
Run the following command

    ./run_TChiHH_mChi-400.sh <nevents> <random_seed> <gridpack_folder>
    
This will produce `TChiHH_mChi-400_RS-1_madgraphMLM-pythia8_RunIISpring15MiniAODv2-FastAsympt25ns_74X_MINIAODSIM.root`
