#!/bin/bash

NEVENTS=$1
RANDOM_SEED=$2
OUTDIR=$3

GRIDPACK_PATH=$PWD
# GRIDPACK_PATH=/home/users/ana/gridpacks
PROCESS="SMS-TChiHH_4b_mChi-400"
MODEL="# model SMS-TChiHH_4b_400_1"
TAG=""
CHI_DECAY=\
"DECAY   1000023     1.00000000E+00   # Chi_20 decays
     1.00000000E+00    2     1000022         25
DECAY   1000025     1.00000000E+00   # Chi_30 decays
     1.00000000E+00    2     1000022         25"
QCUT=84

source /code/osgcode/cmssoft/cmsset_default.sh
# on lxplus only
# export X509_USER_PROXY=$HOME/private/personal/voms_proxy.cert

export SCRAM_ARCH=slc6_amd64_gcc481
if [ -r CMSSW_7_1_20_patch3/src ] ; then 
 echo release CMSSW_7_1_20_patch3 already exists
else
scram p CMSSW CMSSW_7_1_20_patch3
fi
cd CMSSW_7_1_20_patch3/src
eval `scram runtime -sh`

scram b

tar -xaf ${GRIDPACK_PATH}/${PROCESS}_tarball.tar.xz

echo "Running event generation..."
./runcmsgrid.sh $NEVENTS $RANDOM_SEED $(getconf _NPROCESSORS_ONLN)
cd ../../

awk '/<MGGenerationInfo/{n++}{print >"test"n".lhe"}' CMSSW_7_1_20_patch3/src/cmsgrid_final.lhe
rm -rf CMSSW_7_1_20_patch3/src/*
rm test.lhe

echo \
"<LesHouchesEvents version=\"3.0\">
<header>
<MGVersion>
# MG/ME version    : 2.2.2
</MGVersion>
<MG5ProcCard>
import model mssm
define p u u~ d d~ s s~ c c~ b b~ g
define j p
define q p
generate p p > go go  @1
add process p p >  go go j  @2
add process p p > go go j j @3
# Specify process(es) to run
#slepton-slepton
<MGProcCard>
# Begin PROCESS # This is TAG. Do not modify this line
pp>T1tttt    @0       # First Process
# End PROCESS  # This is TAG. Do not modify this line
</MGProcCard>
</MG5ProcCard>
<slha>
#
BLOCK MASS  # Mass Spectrum
# PDG code           mass       particle
   1000001     1.00000000E+05   # ~d_L
   2000001     1.00000000E+05   # ~d_R
   1000002     1.00000000E+05   # ~u_L
   2000002     1.00000000E+05   # ~u_R
   1000003     1.00000000E+05   # ~s_L
   2000003     1.00000000E+05   # ~s_R
   1000004     1.00000000E+05   # ~c_L
   2000004     1.00000000E+05   # ~c_R
   1000005     1.00000000E+05   # ~b_1
   2000005     1.00000000E+05   # ~b_2
   1000006     1.00000000E+05   # ~t_1
   2000006     1.00000000E+05   # ~t_2
   1000011     1.00000000E+04   # ~e_L
   2000011     1.00000000E+04   # ~e_R
   1000012     1.00000000E+04   # ~nu_eL
   1000013     1.00000000E+04   # ~mu_L
   2000013     1.00000000E+04   # ~mu_R
   1000014     1.00000000E+04   # ~nu_muL
   1000015     1.00000000E+04   # ~tau_1
   2000015     1.00000000E+04   # ~tau_2
   1000016     1.00000000E+04   # ~nu_tauL
   1000021     1.00000000E+04   # ~g
   1000022     1.               # ~chi_10
   1000023     400.             # ~chi_20
   1000025     400.             # ~chi_30
   1000035     1.00000000E+04   # ~chi_40
   1000024     1.00000000E+04   # ~chi_1+
   1000037     1.00000000E+04   # ~chi_2+
#
#
#
#         PDG            Width
DECAY   1000006     0.00000000E+00   # stop1 decays
DECAY   2000006     0.00000000E+00   # stop2 decays
DECAY   1000005     0.00000000E+00   # sbottom1 decays
DECAY   2000005     0.00000000E+00   # sbottom2 decays
#
#         PDG            Width
DECAY   1000011     0.00000000E+00   # selectron_L decays
DECAY   2000011     0.00000000E+00   # selectron_R decays
DECAY   1000013     0.00000000E+00   # smuon_L decays
DECAY   2000013     0.00000000E+00   # smuon_R decays
DECAY   1000015     0.00000000E+00   # stau_1 decays
DECAY   2000015     0.00000000E+00   # stau_2 decays
#
#         PDG            Width
DECAY   1000012     0.00000000E+00   # snu_elL decays
DECAY   1000014     0.00000000E+00   # snu_muL decays
DECAY   1000016     0.00000000E+00   # snu_tauL decays
$CHI_DECAY
DECAY   1000021     0.00000000E+00   # gluino decays
DECAY   1000022     0.00000000E+00   # neutralino1 decays
DECAY   1000024     0.00000000E+00   # chargino1 decays

</slha>

<MGRunCard>
   1   = ickkw   ! turning matching on/off for multi-jet sample
  30   = xqcut   ! minimum kt jet measure between partons
</MGRunCard>
" >> undecayed.lhe

cat test1.lhe >> undecayed.lhe

sed -i "s|<scales pt_clust_1|$MODEL\n<scales pt_clust_1|g" undecayed.lhe

echo "Running pLHE step..."
cmsDriver.py \
    step1 --filein file:undecayed.lhe \
    --fileout file:${PROCESS}_plhe.root \
    --mc --eventcontent LHE --datatier LHE \
    --conditions MCRUN2_71_V1::All \
    --step NONE --python_filename plhe_driver.py \
    --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1 || exit $? ; 
cmsRun -e -j plhe_driver_rt.xml plhe_driver.py || exit $? ; 



####################################   GEN-SIM  ####################################
#========== https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_test/SUS-RunIIWinter15GS-00160

echo "============  Making GENSIM: showering..."

cd CMSSW_7_1_20_patch3/src

mkdir -p Configuration/GenProduction/python

echo \
"import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

generator = cms.EDFilter(\"Pythia8HadronizerFilter\",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring(
            'JetMatching:setMad = off',
            'JetMatching:scheme = 1',
            'JetMatching:merge = on',
            'JetMatching:jetAlgorithm = 2',
            'JetMatching:etaJetMax = 5.',
            'JetMatching:coneRadius = 1.',
            'JetMatching:slowJetPower = 1',
            'JetMatching:qCut = $QCUT.', #this is the actual merging scale                        
            'JetMatching:nQmatch = 5', #4 corresponds to 4-flavour scheme (no matching of b-quarks), 5 for 5-flavour scheme
            'JetMatching:nJetMax = 2', #number of partons in born matrix element for highest multiplicity         
            'JetMatching:doShowerKt = off', #off for MLM matching, turn on for shower-kT matching  
            '25:m0 = 125.0',      # Set the mass of the Higgs to 125 GeV
            '25:onMode = off',    # Turn off all higgs decays 
            '25:onIfAny = 5 -5',  # Turn on decays that include a b or a bbar
        ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'processParameters',
                                    )
    )
)" > Configuration/GenProduction/python/genfragment.py

scram b
cd ../../
cmsDriver.py Configuration/GenProduction/python/genfragment.py \
        --filein file:${PROCESS}_plhe.root \
	--fileout file:${PROCESS}_RS-${RANDOM_SEED}_GEN.root \
	--mc --eventcontent RAWSIM \
        --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring \
        --customise_command "process.source.setRunNumber = cms.untracked.uint32($RANDOM_SEED)" \
	--datatier GEN-SIM \
	--conditions MCRUN2_71_V1::All --beamspot NominalCollision2015 --step GEN --magField 38T_PostLS1 \
	--python_filename ${PROCESS}_gen_cfg.py \
	--no_exec -n -1 || exit $? ; 

cmsRun -e -j SUS-RunIIWinter15GS-00160_rt.xml ${PROCESS}_gen_cfg.py || exit $? ; 

#===================================================================================================

## Copy output to hadoop
lcg-cp -b -D srmv2 --vo cms -t 2400 --verbose file:${PROCESS}_RS-${RANDOM_SEED}_GEN.root srm://bsrm-3.t2.ucsd.edu:8443/srm/v2/server?SFN=${OUTDIR}/${PROCESS}_RS-${RANDOM_SEED}_GEN.root

rm -rf CMSSW_7_1_20_patch3
rm *root *lhe *xml *py

echo "Bye."
