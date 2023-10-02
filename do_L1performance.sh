g++ L1performance_HI2023.cpp -o Execute_L1performance `root-config --cflags --libs` -lASImage

#dataDir=/eos/cms/store/group/phys_heavyions/jmijusko/run3RapidValidation/PbPb2023_run374322_PhysicsHIPhysicsRawPrime0_withDFinder_2023-09-28/0000/
#dataDir=/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/PbPb2023_run374354_HIPhysicsRawPrime0_triggerObjects_2023-09-29/CRAB_UserFiles/crab_PbPb2023_run374354_HIPhysicsRawPrime0_triggerObjects_2023-09-29/230930_011735/0000/
dataDir="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime[3,4,5,6,8,9]/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime[3,4,5,6,8,9]_374354_Dpt2trk1/*/0000/"

runtext="Run 374354 (5.36 TeV)"
folder="figs/20231001/Run374354_HIPhysicsRawPrime345689"
suffixText="_374354_prompt"
drMax=0.5
logpath="logs/Run2023_374354_L1pmt"

#./Execute_L1performance --trigType 0 --isL1Obj true --nocut false --isMC true --folder $folder > logs/test_MC_0.log 2> logs/test_MC_0.err < /dev/null &

./Execute_L1performance --trigType 0  --isL1Obj true --nocut false --isMC false --trigsuf "_BptxAND" --PSvec 1,1,1 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_0.log  2>${logpath}_0.err  < /dev/null &
./Execute_L1performance --trigType 1  --isL1Obj true --nocut false --isMC false --trigsuf "_BptxAND" --PSvec 3 	   --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_1.log  2>${logpath}_1.err  < /dev/null &
./Execute_L1performance --trigType 2  --isL1Obj true --nocut false --isMC false --trigsuf "_BptxAND" --PSvec 3 	   --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_2.log  2>${logpath}_2.err  < /dev/null &
./Execute_L1performance --trigType -1 --isL1Obj true --nocut false --isMC false --trigsuf "" 		 --PSvec 631   --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_m1.log 2>${logpath}_m1.err < /dev/null &
