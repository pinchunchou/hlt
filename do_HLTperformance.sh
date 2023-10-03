#!/bin/sh
g++ HLTperformance_HI2023.cpp -o Execute_HLTperformance `root-config --cflags --libs` -lASImage

#./Execute_HLTperformance --trigType 0 > logs/PbPb0.log 2> logs/PbPb0.err < /dev/null &
#./Execute_HLTperformance --trigType 1 > logs/PbPb1.log 2> logs/PbPb1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 > logs/PbPb2.log 2> logs/PbPb2.err < /dev/null &
#./Execute_HLTperformance --trigType 3 > logs/PbPb3.log 2> logs/PbPb3.err < /dev/null &
#./Execute_HLTperformance --trigType 4 > logs/PbPb4.log 2> logs/PbPb4.err < /dev/null &
#./Execute_HLTperformance --trigType 5 > logs/PbPb5.log 2> logs/PbPb5.err < /dev/null &
#./Execute_HLTperformance --trigType 6 > logs/PbPb6.log 2> logs/PbPb6.err < /dev/null &
#./Execute_HLTperformance --trigType 7 > logs/PbPb7.log 2> logs/PbPb7.err < /dev/null &
#./Execute_HLTperformance --trigType 8 > logs/PbPb8.log 2> logs/PbPb8.err < /dev/null &

#./Execute_HLTperformance --trigType 0 --isPbPb false > logs/pp0.log 2> logs/pp0.err < /dev/null &
#./Execute_HLTperformance --trigType 1 --isPbPb false > logs/pp1.log 2> logs/pp1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 --isPbPb false > logs/pp2.log 2> logs/pp2.err < /dev/null &
#./Execute_HLTperformance --trigType 3 --isPbPb false > logs/pp3.log 2> logs/pp3.err < /dev/null &
#./Execute_HLTperformance --trigType 4 --isPbPb false > logs/pp4.log 2> logs/pp4.err < /dev/null &
#./Execute_HLTperformance --trigType 5 --isPbPb false > logs/pp5.log 2> logs/pp5.err < /dev/null &
#./Execute_HLTperformance --trigType 6 --isPbPb false > logs/pp6.log 2> logs/pp6.err < /dev/null &
#./Execute_HLTperformance --trigType 7 --isPbPb false > logs/pp7.log 2> logs/pp7.err < /dev/null &

#./Execute_HLTperformance --trigType 0 --noL1 true > logs/PbPb0_noL1.log 2> logs/PbPb0_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 1 --noL1 true > logs/PbPb1_noL1.log 2> logs/PbPb1_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 --noL1 true > logs/PbPb2_noL1.log 2> logs/PbPb2_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 3 --noL1 true > logs/PbPb3_noL1.log 2> logs/PbPb3_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 8 --noL1 true > logs/PbPb8_noL1.log 2> logs/PbPb8_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 0 --isPbPb false --noL1 true > logs/pp0_noL1.log 2> logs/pp0_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 1 --isPbPb false --noL1 true > logs/pp1_noL1.log 2> logs/pp1_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 --isPbPb false --noL1 true > logs/pp2_noL1.log 2> logs/pp2_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 3 --isPbPb false --noL1 true > logs/pp3_noL1.log 2> logs/pp3_noL1.err < /dev/null &

#dataDir=/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/PbPb2023_run374289_HIPhysicsRawPrime0_withDFinder_2023-09-26/0000/
#dataDir=/eos/cms/store/group/phys_heavyions/jmijusko/run3RapidValidation/PbPb2023_run374322_PhysicsHIPhysicsRawPrime0_withDFinder_2023-09-28/0000/
#dataDir=/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/PbPb2023_run374345_HIPhysicsRawPrime0_triggerObjects_2023-09-29/CRAB_UserFiles/crab_PbPb2023_run374345_HIPhysicsRawPrime0_triggerObjects_2023-09-29/230930_011704/0000/
#dataDir="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime[3,4,5,6,8,9]/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime[3,4,5,6,8,9]_374354_Dpt2trk1/*/0000/"
#dataDir="/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/HIPhysicsRawPrime0_HIRun2023A-PromptReco-v1_run374322_2023-09-30/0000/"

dataDir=\"/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime9/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime9_374354_Dpt2trk1/231001_011507/0000/*\",
dataDir+=\"/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime8/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime8_374354_Dpt2trk1/231001_011439/0000/*\",
dataDir+=\"/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime6/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime6_374354_Dpt2trk1/231001_011111/0000/*\",
dataDir+=\"/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime5/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime5_374354_Dpt2trk1/231001_010540/0000/*\",
dataDir+=\"/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime4/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime4_374354_Dpt2trk1/231001_010449/0000/*\",
dataDir+=\"/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime3/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime3_374354_Dpt2trk1/231001_010436/0000/*\",
dataDir+=\"/eos/cms/store/group/phys_heavyions/mstojano/run3RapidValidation/PbPb2023_run374345_HIPhysicsRawPrime1_2023-09-29/HIPhysicsRawPrime1/crab_PbPb2023_run374345_HIPhysicsRawPrime1_2023-09-29/230930_151533/0000/*\"
#dataDir+=\"/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/HIPhysicsRawPrime0_HIRun2023A-PromptReco-v1_run374322_2023-09-30/0000/*\"


runtext="Run 374322,374345,374354 (5.36 TeV)"
#folder="figs/20230928/Run374289_HIPhysicsRawPrime0"
folder="figs/20231003/Run374322_345_354_PromptReco"
logpath="logs/Run2023_374322_345_354_pmt"
suffixText="_374322_345_354_prompt"
drMax=0.5

./Execute_HLTperformance --trigType 0 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 4,4,5 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_0.log 2>${logpath}_0.err < /dev/null &
./Execute_HLTperformance --trigType 1 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1700,10,1 --LPSvec 3,3,3 --L1ID 2,2,2 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_1.log 2>${logpath}_1.err < /dev/null &
./Execute_HLTperformance --trigType 2 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 4,4,5 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_2.log 2>${logpath}_2.err < /dev/null &
./Execute_HLTperformance --trigType 3 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 650,5,1   --LPSvec 3,3,3 --L1ID 2,2,2 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_3.log 2>${logpath}_3.err < /dev/null &
./Execute_HLTperformance --trigType 4 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 20,7,1    --LPSvec 3,3,1 --L1ID 2,2,3 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_4.log 2>${logpath}_4.err < /dev/null &
./Execute_HLTperformance --trigType 5 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 3,4,4 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_5.log 2>${logpath}_5.err < /dev/null &
./Execute_HLTperformance --trigType 6 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1,1,1     --LPSvec 3,3,3 --L1ID 6,6,6 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_6.log 2>${logpath}_6.err < /dev/null &
./Execute_HLTperformance --trigType 7 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1,1,1     --LPSvec 3,3,3 --L1ID 6,6,6 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_7.log 2>${logpath}_7.err < /dev/null &
./Execute_HLTperformance --trigType 8 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "1" --PSvec 1,1,1     --LPSvec 3,3,3 --L1ID 6,6,6 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_8.log 2>${logpath}_8.err < /dev/null &