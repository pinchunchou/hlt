#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/CorrelationTuple/EventMatcher.h"
#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/TreeHeaders/ggHiNtuplizerTree.h"
#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/TreeHeaders/l1ObjectTree.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"

#include <TFile.h>
#include <TTree.h>
#include <TPad.h>
#include <TROOT.h>
#include <cmath>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TChain.h>
#include <TLine.h>
#include <TColor.h>
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TSystem.h>

#include <iostream>
#include <iomanip>

#include <glob.h>

#include "/afs/cern.ch/user/p/pchou/PhysicsHIZHadron2022/CommonCode/include/CommandLine.h"

using namespace std;

std::vector<int> GetPrimaryColors();
std::vector<std::string> my_glob(const char *pattern);

//void L1performance_HI2023(int trigType)
int main(int argc, char *argv[])
{

  CommandLine CL(argc, argv);

  int trigType    = CL.GetInt("trigType", 0);
  bool isMC       = CL.GetBool("isMC", true);
  bool isPbPb     = CL.GetBool("isPbPb", true);
  bool noL1       = CL.GetBool("noL1", false);

  bool isL1Obj    = CL.GetBool("isL1Obj", false);

  bool nocut      = CL.GetBool("nocut", false);
  
  string folder   = CL.Get("folder", "figs/20230929/");
  string dataDir  = CL.Get("dataDir", "/eos/cms/store/group/phys_heavyions/wangj/RECO2023/miniaod_HIExpressRawPrime_374322/");
  string filename = CL.Get("filename", "L1Eff_");
  string runtext  = CL.Get("runtext", "Run #### data");
  string suffix   = CL.Get("suffix", "");
  string trigsuf  = CL.Get("trigsuf", "");
  double dRmax    = CL.GetDouble("dRmax", 0.5);
  double GendRmax = CL.GetDouble("GendRmax", 0.1);

  vector<int> PSvec = CL.GetIntVector("PSvec",vector<int>{1,1,1,1});

  vector<string> triggers;      
  vector<string> trigger_names; 

  if(trigType==-1){
    triggers = CL.GetStringVector("triggers",vector<string>{"L1_SingleEG7_BptxAND","L1_SingleEG15_BptxAND","L1_SingleEG21_BptxAND","L1_SingleEG30_BptxAND"});
    trigger_names = CL.GetStringVector("trigger_names",vector<string>{"L1 Single EG7 BptxAND","L1 SingleEG 15 BptxAND","L1 SingleEG 21 BptxAND","L1 SingleEG 30 BptxAND"});
  }

  double endcap_legendpp       = CL.GetDouble("endcap_legendpp",0);
  double doublee_legshift      = CL.GetDouble("doublee_legshift",0);

  bool isEB           = CL.GetBool("isEB", false);
  bool isEC           = CL.GetBool("isEC", false);
  bool isEle          = CL.GetBool("isEle", false);
  bool isDoubleEle    = CL.GetBool("isDoubleEle", false);
  bool isMass50       = CL.GetBool("isMass50", false);
  bool isDoublePhoton = CL.GetBool("isDoublePhoton", false);

  int pt_min    = CL.GetInt("pt_min", 5);
  int pt_max    = CL.GetInt("pt_max", 50);
  int pt_bin    = CL.GetInt("pt_bin", 25);
  int ax_min    = CL.GetInt("ax_min", 5);
  int ax_max    = CL.GetInt("ax_max", 55);

  if(isPbPb&&isMC)
    lumi_sqrtS = "PbPb, #sqrt{s_{NN}} = 5.36 TeV, Run 3 Simulation";
  else if(isPbPb==false&&isMC)
    lumi_sqrtS = "ppref, #sqrt{s_{NN}} = 5.02 TeV, Run 3 Simulation";
  else if(isPbPb&&isMC==false)
    lumi_sqrtS = "PbPb data, "  + runtext;
  else if(isPbPb==false&&isMC==false)
    lumi_sqrtS = "ppref data, " + runtext;
  else
    lumi_sqrtS = runtext;


  if(isMC)
    filename += "MC_";
  else
    filename += "data_";

  if(isPbPb)
    filename += "PbPb_";
  else
    filename += "ppref_";

  setTDRStyle();
  gStyle->SetLegendBorderSize(0);

  vector<int> Colors = GetPrimaryColors();
  int marks[10] = {20,21,33,22,23,24,25,27,26,32};

  writeExtraText = true;       // if extra text
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=11;

  string pp_HI;
  if(isPbPb)
    pp_HI = "HI";
  else
    pp_HI = "PPRef";


 if(trigType==-1){
    triggers = CL.GetStringVector("triggers",vector<string>{});
    trigger_names = CL.GetStringVector("trigger_names",vector<string>{});
  }

  if(trigType==0){
    //triggers.push_back("L1_SingleEG7" + trigsuf);
    triggers.push_back("L1_SingleEG15" + trigsuf);
    triggers.push_back("L1_SingleEG21" + trigsuf); 
    triggers.push_back("L1_SingleEG30" + trigsuf); 
    //trigger_names.push_back("L1 SingleEG 7");
    trigger_names.push_back("L1 SingleEG 15");
    trigger_names.push_back("L1 SingleEG 21");
    trigger_names.push_back("L1 SingleEG 30");
    filename+="SingleEG";
  }else if(trigType==1){
    triggers.push_back("L1_DoubleEG5" + trigsuf);
    trigger_names.push_back("L1 DoubleEG 5");
    filename+="DoubleEG";
    isDoublePhoton=true;
  }else if(trigType==2){
    triggers.push_back("L1_SingleEG7" + trigsuf);
    trigger_names.push_back("L1 SingleEG 7");
    filename+="SingleEG7";
  }else{
    triggers.push_back("L1_ZeroBias");
    trigger_names.push_back("L1 ZeroBias");
    filename+="ZeroBias";
  }

  if(noL1&&isEle==false&&isMC)
    filename+="_noL1";

  if(isEB)
    filename+="_barrel";

  if(isEC)
    filename+="_endcap";

  filename+=suffix;

  string frtfile;

  if(isEle==false&&isMC&&isPbPb)
    frtfile = "~/eos/HIRECO_MINIAOD_CMSSW1321/PbPb_MC_QCDPhoton/Pythia8_QCDPhoton15_HydjetEmbedded_TuneCP5_13_0_12/PbPb_MC_QCDPhoton_13_0_12_mAOD_new/230830_145820/0000/HiForestMiniAOD_*.root";
  else if(isEle&&isMC&&isPbPb)
    frtfile = "~/eos/HIRECO_MINIAOD_CMSSW1321/PbPb_MC_Zee/Pythia8_Ze10e10_Embedded_TuneCP5_13_0_12/PbPb_MC_Zee_13_0_12_mAOD_new/230830_145820/0000/HiForestMiniAOD_*.root";
  else if(isEle==false&&isMC&&isPbPb==false)
    frtfile = "~/eos/HIRECO_MINIAOD_CMSSW1321/ppref_MC_QCDPhoton/QCDPhoton_TuneCP5_5p36TeV_ppref-pythia8_final_200M/QCDPhoton_TuneCP5_5p36TeV_ppref-pythia8_new_mAOD/230909_162228/0000/*.root";
  else if(isEle&&isMC&&isPbPb==false)
    frtfile = "~/eos/HIRECO_MINIAOD_CMSSW1321/ppref_MC_ZeePU/Zee_TuneCP5_5p36TeV_ppref-pythia8_final_2M/Zee_TuneCP5_5p36TeV_ppref-pythia8_new_mAOD/230821_123934/0000/*.root";
  else
    frtfile = dataDir + "/*.root";

  string hltfile;

  if(noL1==false&&isEle==false&&isMC&&isPbPb)
    hltfile = "~/eos_base/HLT_DIGI_CMSSW1321/PbPb_MC_QCDPhoton/Macro/CRAB_UserFiles/PbPb_MC_QCDPhoton_1320V33_Macro/230912_133559/0000/*.root";
  else if(noL1&&isEle==false&&isMC&&isPbPb)
    hltfile = "~/eos_base/HLT_DIGI_CMSSW1321/PbPb_MC_QCDPhoton/Macro/CRAB_UserFiles/PbPb_MC_QCDPhoton_1320V33_Macro_noL1/230912_134947/0000/*.root";
  else if(isEle&&isMC&&isPbPb)
    hltfile = "~/eos_base/HLT_DIGI_CMSSW1321/PbPb_MC_Zee/Macro/CRAB_UserFiles/PbPb_MC_Zee_1320V33_Macro/230912_134817/0000/*.root";
  
  else if(noL1==false&&isEle==false&&isMC&&isPbPb==false)
    hltfile = "~/eos_base/HLT_DIGI_CMSSW1321/ppref_MC_QCDPhoton_1320V32/Macro/CRAB_UserFiles/ppref_MC_QCDPhoton_1320V32_tag132X2023_Macro/230911_121651/0000/*.root";
  else if(noL1&&isEle==false&&isMC&&isPbPb==false)
    hltfile = "~/eos/HLT_DIGI_CMSSW1321/ppref_MC_QCDPhoton_1320V32/Macro/CRAB_UserFiles/ppref_MC_QCDPhoton_1320V32_tag132X2023_noL1_Macro/230911_161027/0000/*.root";
  else if(isEle&&isMC&&isPbPb==false)
    hltfile = "~/eos_base/HLT_DIGI_CMSSW1321/ppref_MC_Zee_1320V32/Macro/CRAB_UserFiles/ppref_MC_Zee_1320V32_tag132X2023_Macro/230911_112331/0000/*.root";  
  else
    hltfile = dataDir + "/*.root";

  std::cout << "Reading files..." << std::endl;

  TChain *HltTree = new TChain("hltanalysis/HltTree");
  for (const auto &filename : my_glob(hltfile.c_str()))
    HltTree->Add(filename.c_str());

  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");

  for (const auto &filename : my_glob(frtfile.c_str())){
    EventTree->Add(filename.c_str());
    HiTree->Add(filename.c_str());
  }

  TChain *L1ObjTree = new TChain("l1object/L1UpgradeFlatTree");
  if(isL1Obj){
    for (const auto &filename : my_glob(hltfile.c_str()))
      L1ObjTree->Add(filename.c_str());
    //L1ObjTree->Add(hltfile.c_str());
    L1ObjTree->SetBranchStatus("*",0);
    L1ObjTree->SetBranchStatus("egEta", 1);
    L1ObjTree->SetBranchStatus("egPhi", 1);
    L1ObjTree->SetBranchStatus("egEt", 1);
  }
  //vector <TChain *> Trees;
  vector <TH1F *> hnums;
  vector <TH1F *> hdenoms;

  //TH1F* hdenom   = new TH1F("hdenom","",pt_bin,pt_min,pt_max);

  HltTree->SetBranchStatus("*",0);
  HltTree->SetBranchStatus("Event", 1);
  HltTree->SetBranchStatus("LumiBlock", 1);
  HltTree->SetBranchStatus("Run", 1);

  //Bool_t HLT_Trig[100];
  Int_t HLT_Trig[100], PSs[100];

  for(long unsigned int iT = 0; iT < triggers.size(); iT++){
    hdenoms.push_back(new TH1F(("hdenum"+ to_string(iT)).c_str(),"", pt_bin, pt_min, pt_max));
    hnums.push_back(new TH1F(("hnum"+ to_string(iT)).c_str(),"", pt_bin, pt_min, pt_max));

    std::cout<<"triggers["<<iT<<"] = "<<triggers[iT]<<std::endl;
    HltTree->SetBranchStatus(triggers[iT].c_str(), 1);
    HltTree->SetBranchAddress(triggers[iT].c_str(), &HLT_Trig[iT]);
    HltTree->SetBranchStatus((triggers[iT]+"_Prescl").c_str(), 1);
    HltTree->SetBranchAddress((triggers[iT]+"_Prescl").c_str(), &PSs[iT]);
  }

  HiTree->SetBranchStatus("*",0);
  HiTree->SetBranchStatus("hiBin",1);
  HiTree->SetBranchStatus("hiHF",1);
  if(isMC)
    HiTree->SetBranchStatus("weight",1);

  ULong64_t       hlt_event;
  Int_t           hlt_lumi, hlt_run, hiBin;
  Float_t         hiHF, pthat_weight;

  HiTree->SetBranchAddress("hiBin", &hiBin);
  HiTree->SetBranchAddress("hiHF", &hiHF);
  if(isMC)
    HiTree->SetBranchAddress("weight", &pthat_weight);

  HltTree->SetBranchAddress("Event", &hlt_event);
  HltTree->SetBranchAddress("LumiBlock", &hlt_lumi);
  HltTree->SetBranchAddress("Run", &hlt_run);

  EventTree->SetBranchStatus("*",0);
  EventTree->SetBranchStatus("run",1);
  EventTree->SetBranchStatus("event",1);
  EventTree->SetBranchStatus("lumis",1);
  
  if(isEle==false){
    EventTree->SetBranchStatus("nPho",1);
    if(isMC)
      EventTree->SetBranchStatus("pho_genMatchedIndex",1); 
  
    EventTree->SetBranchStatus("phoEt",1);
    EventTree->SetBranchStatus("phoEta",1);
    EventTree->SetBranchStatus("phoHoverE",1);
    EventTree->SetBranchStatus("phoPhi",1);
    EventTree->SetBranchStatus("phoSCEta",1);
    EventTree->SetBranchStatus("phoSCPhi",1);
    EventTree->SetBranchStatus("phoSigmaIEtaIEta_2012",1);
  
    EventTree->SetBranchStatus("pfpIso3subUE",1);
    EventTree->SetBranchStatus("pfcIso3subUE",1);
    EventTree->SetBranchStatus("pfnIso3subUE",1);
  
    if(isPbPb){
      EventTree->SetBranchStatus("pho_swissCrx",1);
      EventTree->SetBranchStatus("pho_seedTime",1);
      EventTree->SetBranchStatus("pho_ecalClusterIsoR3",1);
      EventTree->SetBranchStatus("pho_hcalRechitIsoR3",1);
      EventTree->SetBranchStatus("pho_trackIsoR3PtCut20",1);
    }

  }else{
    EventTree->SetBranchStatus("nEle",1);

    EventTree->SetBranchStatus("elePt",1);
    EventTree->SetBranchStatus("eleEta",1);
    EventTree->SetBranchStatus("eleHoverE",1);
    EventTree->SetBranchStatus("elePhi",1);
  
    EventTree->SetBranchStatus("eleSigmaIEtaIEta_2012",1);
    EventTree->SetBranchStatus("eledEtaSeedAtVtx",1);
    EventTree->SetBranchStatus("eledPhiAtVtx",1);
    EventTree->SetBranchStatus("eleHoverEBc",1);
    EventTree->SetBranchStatus("eleEoverPInv",1);
    EventTree->SetBranchStatus("eleMissHits",1);
    EventTree->SetBranchStatus("eleIP3D",1);
  
    EventTree->SetBranchStatus("eleSCEta",1);
    EventTree->SetBranchStatus("eleSCRawEn",1);
  }
   
  if(isMC){ 
    EventTree->SetBranchStatus("mcPID",1);
    EventTree->SetBranchStatus("mcCalIsoDR04",1);
    EventTree->SetBranchStatus("mcIndex",1);
    EventTree->SetBranchStatus("mcEta",1);
    EventTree->SetBranchStatus("mcPhi",1);
  }


  EventMatcher* emTrig = new EventMatcher();
  EventMatcher* em = new EventMatcher();

  if(isMC){
    Long64_t duplicateEntriesHlt = 0;
    Long64_t entriesAnalyzedHlt = 0;
  
    Long64_t entriesHlt = HltTree->GetEntries();
    std::cout << "entries in HLT = " << entriesHlt << std::endl;
    std::cout << "Loop over the file with HLT emulation output..." << std::endl;
  
    for (Long64_t j_entry = 0; j_entry < entriesHlt; ++j_entry){
      if (j_entry % 10000 == 0)  {
        std::cout << "current entry = " <<j_entry<< " out of " <<entriesHlt<< " : " <<std::setprecision(2)<<(double)j_entry/entriesHlt*100<< " %" << std::endl;
      }
      HltTree->GetEntry(j_entry);
  
      bool eventAdded = emTrig->addEvent(hlt_run, hlt_lumi, hlt_event, j_entry);
      if(!eventAdded) // this event is duplicate, skip this one.
      {
        duplicateEntriesHlt++;
        continue;
      }
      entriesAnalyzedHlt++;
    }
  
    std::cout << "###" << std::endl;
    std::cout << "Loop HLT ENDED" << std::endl;
    std::cout << "entries HLT          = " << entriesHlt << std::endl;
    std::cout << "duplicateEntries HLT = " << duplicateEntriesHlt << std::endl;
    std::cout << "entriesAnalyzed HLT  = " << entriesAnalyzedHlt << std::endl;
    std::cout << "###" << std::endl;
  }

  std::cout << "Loop over the file with offline reco output..." << std::endl;

  ggHiNtuplizer ggHi;
  ggHi.setupTreeForReading(EventTree);

  Long64_t entriesTmp = EventTree->GetEntries();
  std::cout << "entries in File = " << entriesTmp << std::endl;

  Long64_t entryTrig = 0;
  Int_t entriesNotFoundinTrigger=0;
  Int_t duplicateEntries=0, entriesAnalyzed=0;

  Int_t counts[100]={0};
  Int_t count0 = 0;

  for (Long64_t j_entry = 0; j_entry < entriesTmp; ++j_entry){
    if (j_entry % 10000 == 0)  {
        std::cout << "current entry = " <<j_entry<< " out of " <<entriesTmp<< " : " <<std::setprecision(2)<<(double)j_entry/entriesTmp*100<< " %" << std::endl;
    }

    EventTree->GetEntry(j_entry);
    HiTree->GetEntry(j_entry);


    bool eventAdded = em->addEvent(ggHi.run, ggHi.lumis, ggHi.event, j_entry);
    if(!eventAdded) // this event is duplicate, skip this one.
    {
        duplicateEntries++;
        continue;
    }

    if(isMC){
      entryTrig = emTrig->getEntry(ggHi.run, ggHi.lumis, ggHi.event);

      if (entryTrig < 0) {
        entriesNotFoundinTrigger++;
        continue;
      }

      emTrig->removeEvent(ggHi.run, ggHi.lumis, ggHi.event);
    }else{
      entryTrig = j_entry;
    }

    entriesAnalyzed++;

    //cout<<"entryTrig = "<<entryTrig<<endl;
    HltTree->GetEntry(entryTrig);

    Float_t weighting;

    if(isPbPb&&isMC) 
      weighting = pthat_weight*findNcoll(hiBin);
    else if(isPbPb==false&&isMC) 
      weighting = pthat_weight;
    else
      weighting = 1;
    

    float sumIso, sumIso2, maxPt=0, max2Pt=0;
    int maxPt_i = -1, max2Pt_i = -1;

    //cout<<"Loop into the photons in the event."<<endl;
    //float dRmax = 0.5;
    //float GendRmax = 0.1;
    //float dRmax = 100000000;

    if(isMass50)
      dRmax = 100000000;

    float dr00;

    std::vector<int> ele_genMatchedIndex;
    if(isEle&&isMC){
      for(long unsigned int iele=0;iele<size(*ggHi.elePt);iele++){
        int matchedIndex = -1;
        float minDR = GendRmax;
        for (unsigned igen = 0; igen < size(*ggHi.mcEt); ++igen) {
          if ( abs((*ggHi.mcPID)[igen]) != 11)
            continue;
          float deta = TMath::Abs((*ggHi.eleEta)[iele] - (*ggHi.mcEta)[igen]);
          float dphi = TMath::Abs((*ggHi.elePhi)[iele] - (*ggHi.mcPhi)[igen]);
          if (dphi > TMath::Pi()) 
            dphi = TMath::Abs(dphi - 2*TMath::Pi());
          if (TMath::Sqrt(dphi*dphi + deta*deta) < minDR){ 
            minDR = TMath::Sqrt(dphi*dphi + deta*deta);
            matchedIndex = igen;
          }
        }
        ele_genMatchedIndex.push_back(matchedIndex);
      }
    }

    Int_t num_of_PE;
    if(isEle)
      num_of_PE = ggHi.nEle;
    else
      num_of_PE = ggHi.nPho;

    //std::cout<<"num_of_PE = "<<num_of_PE<<std::endl;
    //std::cout<<"size(*ggHi.phoEt) = "<<size(*ggHi.phoEt)<<std::endl;
    for(int i=0; i < num_of_PE; ++i) {

      if(isMC){
        int64_t genID;
        if(isEle){
          genID = ele_genMatchedIndex[i];
          if (genID == -1)  continue; 
          auto pid = (*ggHi.mcPID)[genID];
          auto mpid = (*ggHi.mcMomPID)[genID];
          if (std::abs(pid) != 11 || std::abs(mpid) != 23) continue;
        }else{
          genID = (*ggHi.pho_genMatchedIndex)[i];
          if (genID == -1) { continue; }
          auto pid = (*ggHi.mcPID)[genID];
          auto mpid = (*ggHi.mcMomPID)[genID];
          if (pid != 22 || (std::abs(mpid) > 22 && mpid != -999)) continue;
          dr00 = sqrt(((*ggHi.mcEta)[genID]-(*ggHi.phoEta)[i])*((*ggHi.mcEta)[genID]-(*ggHi.phoEta)[i])+((*ggHi.mcPhi)[genID]-(*ggHi.phoPhi)[i])*((*ggHi.mcPhi)[genID]-(*ggHi.phoPhi)[i]));
          if(dr00 > GendRmax) continue;
        }      
      }

      float Et;
      if(isEle)
        Et = (*ggHi.elePt)[i];
      else
        Et = (*ggHi.phoEt)[i];

      //std::cout<<"Et = "<<Et<<", maxPt = "<<maxPt<<std::endl;
      if(Et>maxPt){
        maxPt = Et; 
        maxPt_i = i;
      }else if(Et>max2Pt&&Et<maxPt){
        max2Pt = Et; 
        max2Pt_i = i;
      }

    }

    if(maxPt_i==-1) continue;
    if(max2Pt_i==-1&&(isDoubleEle||isDoublePhoton))  continue;
    sumIso2 = -10000;

    if(isEle==false&&isPbPb&&isDoublePhoton==false)
      sumIso = (*ggHi.pho_ecalClusterIsoR3)[maxPt_i]+(*ggHi.pho_hcalRechitIsoR3)[maxPt_i]+(*ggHi.pho_trackIsoR3PtCut20)[maxPt_i];
    else if(isEle==false&&isPbPb==false&&isDoublePhoton==false)
      sumIso = (*ggHi.pfpIso3subUE)[maxPt_i]+(*ggHi.pfcIso3subUE)[maxPt_i]+(*ggHi.pfnIso3subUE)[maxPt_i];
    else if(isPbPb&&isDoublePhoton){
      sumIso = (*ggHi.pho_ecalClusterIsoR3)[maxPt_i]+(*ggHi.pho_hcalRechitIsoR3)[maxPt_i]+(*ggHi.pho_trackIsoR3PtCut20)[maxPt_i];
      sumIso2 = (*ggHi.pho_ecalClusterIsoR3)[max2Pt_i]+(*ggHi.pho_hcalRechitIsoR3)[max2Pt_i]+(*ggHi.pho_trackIsoR3PtCut20)[max2Pt_i];
    }else if(isPbPb==false&&isDoublePhoton){
      sumIso = (*ggHi.pfpIso3subUE)[maxPt_i]+(*ggHi.pfcIso3subUE)[maxPt_i]+(*ggHi.pfnIso3subUE)[maxPt_i];
      sumIso2 = (*ggHi.pfpIso3subUE)[max2Pt_i]+(*ggHi.pfcIso3subUE)[max2Pt_i]+(*ggHi.pfnIso3subUE)[max2Pt_i];
    }

    //std::cout<<"Photon Selection"<<std::endl;

    bool sel_cut = false, sel_cut2 = false;

    if(isPbPb){
      if(isEle==false&&isDoublePhoton==false){
        if(isEC==false && abs((*ggHi.phoSCEta)[maxPt_i]) < 1.44 && (*ggHi.phoHoverE)[maxPt_i] < 0.247665 && (*ggHi.phoSigmaIEtaIEta_2012)[maxPt_i] < 0.012186 && (*ggHi.pho_swissCrx)[maxPt_i] < 0.9 && abs((*ggHi.pho_seedTime)[maxPt_i]) < 3 && sumIso <= 11.697505) sel_cut=true;
        else if(isEB==false && abs((*ggHi.phoSCEta)[maxPt_i]) > 1.57 && abs((*ggHi.phoSCEta)[maxPt_i]) < 2.1 && (*ggHi.phoHoverE)[maxPt_i] < 0.398866 && (*ggHi.phoSigmaIEtaIEta_2012)[maxPt_i] < 0.044998 && (*ggHi.pho_swissCrx)[maxPt_i] < 0.9 && abs((*ggHi.pho_seedTime)[maxPt_i]) < 3 && sumIso < 20.911811) sel_cut=true;
        else sel_cut = false;
        sel_cut2 = true;
      }else if(isDoublePhoton){
        if(isEC==false && abs((*ggHi.phoSCEta)[maxPt_i]) < 1.44 && (*ggHi.phoHoverE)[maxPt_i] < 0.247665 && (*ggHi.phoSigmaIEtaIEta_2012)[maxPt_i] < 0.012186 && (*ggHi.pho_swissCrx)[maxPt_i] < 0.9 && abs((*ggHi.pho_seedTime)[maxPt_i]) < 3 && sumIso <= 11.697505) sel_cut=true;
        else if(isEB==false && abs((*ggHi.phoSCEta)[maxPt_i]) > 1.57 && abs((*ggHi.phoSCEta)[maxPt_i]) < 2.1 && (*ggHi.phoHoverE)[maxPt_i] < 0.398866 && (*ggHi.phoSigmaIEtaIEta_2012)[maxPt_i] < 0.044998 && (*ggHi.pho_swissCrx)[maxPt_i] < 0.9 && abs((*ggHi.pho_seedTime)[maxPt_i]) < 3 && sumIso < 20.911811) sel_cut=true;
        else sel_cut = false;
      
        if(isEC==false && abs((*ggHi.phoSCEta)[max2Pt_i]) < 1.44 && (*ggHi.phoHoverE)[max2Pt_i] < 0.247665 && (*ggHi.phoSigmaIEtaIEta_2012)[max2Pt_i] < 0.012186 && (*ggHi.pho_swissCrx)[max2Pt_i] < 0.9 && abs((*ggHi.pho_seedTime)[max2Pt_i]) < 3 && sumIso2 <= 11.697505) sel_cut2=true;
        else if(isEB==false && abs((*ggHi.phoSCEta)[max2Pt_i]) > 1.57 && abs((*ggHi.phoSCEta)[max2Pt_i]) < 2.1 && (*ggHi.phoHoverE)[max2Pt_i] < 0.398866 && (*ggHi.phoSigmaIEtaIEta_2012)[max2Pt_i] < 0.044998 && (*ggHi.pho_swissCrx)[max2Pt_i] < 0.9 && abs((*ggHi.pho_seedTime)[max2Pt_i]) < 3 && sumIso2 < 20.911811) sel_cut2=true;
        else sel_cut2 = false;
    
      }else if(isEle&&isDoubleEle==false){
        if(hiBin>=0 && hiBin < 60){
          if(isEC==false && abs((*ggHi.eleEta)[maxPt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0135 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0038 && abs((*ggHi.eledPhiAtVtx)[maxPt_i] )< 0.0376 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1616 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0177  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = true;
          else if(isEB==false && abs((*ggHi.eleEta)[maxPt_i]) > 1.57 && abs((*ggHi.eleEta)[maxPt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0466 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0063 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.1186 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1317 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0201  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = true;
          else sel_cut = false;
        }else if(hiBin>=60 && hiBin < 200){
          if(isEC==false && abs((*ggHi.eleEta)[maxPt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0107 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0035 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.0327 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1268 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0774  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = true;
          else if(isEB==false && abs((*ggHi.eleEta)[maxPt_i]) > 1.57 && abs((*ggHi.eleEta)[maxPt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0339 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0067 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.0838 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.0977 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0193  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = true;
          else sel_cut = false;
        }else{
          sel_cut = false;
        }
        sel_cut2 = true;
      }else if(isDoubleEle){
        if(hiBin>=0 && hiBin < 60){
          if(isEC==false && abs((*ggHi.eleEta)[maxPt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0135 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0038 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.0376 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1616 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0177  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = true;
          else if(isEB==false && abs((*ggHi.eleEta)[maxPt_i]) > 1.57 && abs((*ggHi.eleEta)[maxPt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0466 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0063 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.1186 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1317 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0201  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = true;
          else sel_cut = false;
        }else if(hiBin>=60 && hiBin < 200){
          if(isEC==false && abs((*ggHi.eleEta)[maxPt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0107 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0035 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.0327 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1268 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0774  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = true;
          else if(isEB==false && abs((*ggHi.eleEta)[maxPt_i]) > 1.57 && abs((*ggHi.eleEta)[maxPt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0339 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0067 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.0838 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.0977 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0193  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = true;
          else sel_cut = false;
        }else{
          sel_cut = false;
        }
      
        if(hiBin>=0 && hiBin < 60){
          if(isEC==false && abs((*ggHi.eleEta)[max2Pt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.0135 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.0038 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.0376 && (*ggHi.eleHoverEBc)[max2Pt_i] < 0.1616 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.0177  && (*ggHi.eleMissHits)[max2Pt_i] <= 1 && (*ggHi.eleIP3D)[max2Pt_i] < 0.03) sel_cut2 = true;
          else if(isEB==false && abs((*ggHi.eleEta)[max2Pt_i]) > 1.57 && abs((*ggHi.eleEta)[max2Pt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.0466 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.0063 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.1186 && (*ggHi.eleHoverEBc)[max2Pt_i] < 0.1317 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.0201  && (*ggHi.eleMissHits)[max2Pt_i] <= 1 && (*ggHi.eleIP3D)[max2Pt_i] < 0.03) sel_cut2 = true;
          else sel_cut2 = false;
        }else if(hiBin>=60 && hiBin < 200){
          if(isEC==false && abs((*ggHi.eleEta)[max2Pt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.0107 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.0035 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.0327 && (*ggHi.eleHoverEBc)[max2Pt_i] < 0.1268 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.0774  && (*ggHi.eleMissHits)[max2Pt_i] <= 1 && (*ggHi.eleIP3D)[max2Pt_i] < 0.03) sel_cut2 = true;
          else if(isEB==false && abs((*ggHi.eleEta)[max2Pt_i]) > 1.57 && abs((*ggHi.eleEta)[max2Pt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.0339 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.0067 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.0838 && (*ggHi.eleHoverEBc)[max2Pt_i] < 0.0977 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.0193  && (*ggHi.eleMissHits)[max2Pt_i] <= 1 && (*ggHi.eleIP3D)[max2Pt_i] < 0.03) sel_cut2 = true;
          else sel_cut2 = false;
        }else{
          sel_cut2 = false;
        }      
      }else{
        sel_cut = true;
        sel_cut2 = true;
      }
    }else{
      if(isEle==false){
        if(isEC==false && abs((*ggHi.phoSCEta)[maxPt_i]) < 1.44 && (*ggHi.phoHoverE)[maxPt_i] < 0.072266 && (*ggHi.phoSigmaIEtaIEta_2012)[maxPt_i] < 0.010806 && sumIso < 0.416894) sel_cut=true;
        else if(isEB==false && abs((*ggHi.phoSCEta)[maxPt_i]) > 1.57 && abs((*ggHi.phoSCEta)[maxPt_i]) < 2.1 && (*ggHi.phoHoverE)[maxPt_i] < 0.032548 && (*ggHi.phoSigmaIEtaIEta_2012)[maxPt_i] < 0.027323 && sumIso < 0.970591) sel_cut=true;
        else sel_cut = false;
        sel_cut2 = true;
      }else if(isEle&&isDoubleEle==false){
        if(isEC==false && abs((*ggHi.eleEta)[maxPt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.01020 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.00327 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.06055 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.52688 ) sel_cut = true;
        else if(isEB==false && abs((*ggHi.eleEta)[maxPt_i]) > 1.57 && abs((*ggHi.eleEta)[maxPt_i]) < 2.1 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.02957 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.00490 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.09622 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.14600 ) sel_cut = true;
        else sel_cut = false;
        sel_cut2 = true;
      }else if(isDoubleEle){
        if(isEC==false && abs((*ggHi.eleEta)[maxPt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.01020 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.00327 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.06055 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.52688 ) sel_cut = true;
        else if(isEB==false && abs((*ggHi.eleEta)[maxPt_i]) > 1.57 && abs((*ggHi.eleEta)[maxPt_i]) < 2.1 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.02957 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.00490 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.09622 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.14600 ) sel_cut = true;
        else sel_cut = false;

        if(isEC==false && abs((*ggHi.eleEta)[max2Pt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.01020 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.00327 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.06055 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.52688 ) sel_cut2 = true;
        else if(isEB==false && abs((*ggHi.eleEta)[max2Pt_i]) > 1.57 && abs((*ggHi.eleEta)[max2Pt_i]) < 2.1 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.02957 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.00490 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.09622 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.14600 ) sel_cut2 = true;
        else sel_cut2 = false;
      }else{
        sel_cut = true;
        sel_cut2 = true;
      }
    }

    if(nocut){
      sel_cut=true;
      sel_cut2=true;
    }

    double emass = 0.0005111;

    if(isDoubleEle){
      TLorentzVector ev1, ev2;  
      ev1.SetPtEtaPhiM((*ggHi.elePt)[maxPt_i] ,(*ggHi.eleEta)[maxPt_i] ,(*ggHi.elePhi)[maxPt_i] ,emass);
      ev2.SetPtEtaPhiM((*ggHi.elePt)[max2Pt_i],(*ggHi.eleEta)[max2Pt_i],(*ggHi.elePhi)[max2Pt_i],emass);
    
      TLorentzVector Zv12 = ev1+ev2;
      double Zmass = Zv12.M();
    
      if(Zmass<60 || Zmass>120) sel_cut = false;
    }

    //std::cout<<"a"<<std::endl;
/*
    if(sel_cut&&sel_cut2){
      //std::cout<<"(*ggHi.phoEt)[maxPt_i] = "<<(*ggHi.phoEt)[maxPt_i]<<std::endl;
      if(isEle==false&&isDoublePhoton==false)
        hdenom->Fill((*ggHi.phoEt)[maxPt_i],weighting); 
      else if(isEle&&isDoubleEle==false)
        hdenom->Fill((*ggHi.elePt)[maxPt_i],weighting); 
      else if(isDoubleEle)
        hdenom->Fill((*ggHi.elePt)[max2Pt_i],weighting); 
      else if(isDoublePhoton)
        hdenom->Fill((*ggHi.phoEt)[max2Pt_i],weighting);  
      count0++;
    }
*/
    //std::cout<<"b"<<std::endl;

    l1Object hObj;
    if(isL1Obj){ 
      hObj.setupTreeForReading(L1ObjTree);
      L1ObjTree->GetEntry(entryTrig);
    }

    for(long unsigned int iT = 0; iT < triggers.size(); iT++){
      
      //std::cout<<"c"<<std::endl;

      float drTrig;

      if(isMC==false&&PSs[iT]!=PSvec[iT])
        continue;

      if(sel_cut&&sel_cut2){
        float n_pho_ho = 0;
        float drtemp;
        drTrig=10000;

        if(isEle==false&&isDoublePhoton==false)
          hdenoms[iT]->Fill((*ggHi.phoEt)[maxPt_i],weighting); 
        else if(isEle&&isDoubleEle==false)
          hdenoms[iT]->Fill((*ggHi.elePt)[maxPt_i],weighting); 
        else if(isDoubleEle)
          hdenoms[iT]->Fill((*ggHi.elePt)[max2Pt_i],weighting); 
        else if(isDoublePhoton)
          hdenoms[iT]->Fill((*ggHi.phoEt)[max2Pt_i],weighting);  
        
        count0++;

        //std::cout<<"triggers[iT] = "<<triggers[iT]<<std::endl;
        //std::cout<<"HLT_Trig[iT] = "<<HLT_Trig[iT]<<std::endl;

        if((*ggHi.phoEt)[maxPt_i]>64&&HLT_Trig[iT]==0){
          cout<<"Run:"<<hlt_run<<", LumiBlock: "<<hlt_lumi<<", Event: "<<hlt_event<<endl;
          std::cout<<triggers[iT]<<" = "<<HLT_Trig[iT]<<std::endl;
          std::cout<<"phoEt = "<<(*ggHi.phoEt)[maxPt_i]<<", phoEta = "<<(*ggHi.phoEta)[maxPt_i]<<", phoPhi = "<<(*ggHi.phoPhi)[maxPt_i]<<std::endl;
          cout<<"phoHoverE = "<<(*ggHi.phoHoverE)[maxPt_i]<<", phoSigmaIEtaIEta_2012 = "<<(*ggHi.phoSigmaIEtaIEta_2012)[maxPt_i]<<endl;
          cout<<"Other photon/electrons:"<<endl;
          for(int i=0; i < num_of_PE; ++i) {
            std::cout<<"phoEt = "<<(*ggHi.phoEt)[i]<<", phoEta = "<<(*ggHi.phoEta)[i]<<", phoPhi = "<<(*ggHi.phoPhi)[i]<<std::endl;
            cout<<"phoHoverE = "<<(*ggHi.phoHoverE)[i]<<", phoSigmaIEtaIEta_2012 = "<<(*ggHi.phoSigmaIEtaIEta_2012)[i]<<endl;
          }

          cout<<"------------"<<endl;
        }

        if(HLT_Trig[iT]>0){
          if(isL1Obj){
          n_pho_ho = size(*hObj.egEta);
            for(int i=0;i<n_pho_ho;i++){
              if(isEle==false&&isDoublePhoton==false)
                drtemp = sqrt(((*hObj.egEta)[i]-(*ggHi.phoEta)[maxPt_i])*((*hObj.egEta)[i]-(*ggHi.phoEta)[maxPt_i])+((*hObj.egPhi)[i]-(*ggHi.phoPhi)[maxPt_i])*((*hObj.egPhi)[i]-(*ggHi.phoPhi)[maxPt_i]));
              else if(isEle&&isDoubleEle==false)
                drtemp = sqrt(((*hObj.egEta)[i]-(*ggHi.eleEta)[maxPt_i])*((*hObj.egEta)[i]-(*ggHi.eleEta)[maxPt_i])+((*hObj.egPhi)[i]-(*ggHi.elePhi)[maxPt_i])*((*hObj.egPhi)[i]-(*ggHi.elePhi)[maxPt_i]));
              else if(isDoubleEle)
                drtemp = sqrt(((*hObj.egEta)[i]-(*ggHi.eleEta)[max2Pt_i])*((*hObj.egEta)[i]-(*ggHi.eleEta)[max2Pt_i])+((*hObj.egPhi)[i]-(*ggHi.elePhi)[max2Pt_i])*((*hObj.egPhi)[i]-(*ggHi.elePhi)[max2Pt_i]));
              else if(isDoublePhoton)
                drtemp = sqrt(((*hObj.egEta)[i]-(*ggHi.phoEta)[max2Pt_i])*((*hObj.egEta)[i]-(*ggHi.phoEta)[max2Pt_i])+((*hObj.egPhi)[i]-(*ggHi.phoPhi)[max2Pt_i])*((*hObj.egPhi)[i]-(*ggHi.phoPhi)[max2Pt_i]));
          
              if (drtemp < drTrig) drTrig=drtemp;
            }
          }

          

          if(drTrig < dRmax){
            if(isEle==false&&isDoublePhoton==false)
              hnums[iT]->Fill((*ggHi.phoEt)[maxPt_i],weighting); 
            else if(isEle&&isDoubleEle==false)
              hnums[iT]->Fill((*ggHi.elePt)[maxPt_i],weighting); 
            else if(isDoubleEle)
              hnums[iT]->Fill((*ggHi.elePt)[max2Pt_i],weighting); 
            else if(isDoublePhoton)
              hnums[iT]->Fill((*ggHi.phoEt)[max2Pt_i],weighting); 
            
            counts[iT]++;
          }
        }

      }

    }
  }

  std::cout << "count0 = "<<count0;
  for(long unsigned int iT = 0; iT < triggers.size(); iT++){
    std::cout << ", counts["<<iT<<"] = "<<counts[iT];
  }
  std::cout << std::endl;

  std::cout << "###" << std::endl;
  std::cout << "Loop Forest ENDED"  << std::endl;
  std::cout << "entries Forest          = " << entriesTmp << std::endl;
  std::cout << "duplicateEntries Forest = " << duplicateEntries << std::endl;
  std::cout << "entriesAnalyzed Forest  = " << entriesAnalyzed << std::endl;
  std::cout << "entriesNotFoundinTrigger  = " << entriesNotFoundinTrigger << std::endl;
  std::cout << "###" << std::endl;

  TCanvas *canvas1 = new TCanvas("canvas1","",1000,800);
  canvas1->SetMargin(0.15, 0.05, 0.15, 0.08);
  canvas1->GetFrame()->SetBorderSize(12);

  vector <TGraphAsymmErrors *> hratios;

  TLegend leg(0.58+doublee_legshift,0.18+endcap_legendpp,0.98,0.4+endcap_legendpp);

  for(long unsigned int iT = 0; iT < triggers.size(); iT++){
    hratios.push_back(new TGraphAsymmErrors());
    hratios[iT]->SetName(("hratio" + to_string(iT)).c_str());
    hratios[iT]->BayesDivide(hnums[iT], hdenoms[iT]);

    if(isMC==false){
      std::cout<<"Trigger: "<<triggers[iT]<<std::endl;
      std::cout<<"PS:"<<PSvec[iT]<<std::endl;
      hratios[iT]->Scale(PSvec[iT]);
    }

    if(isEle)
      hratios[iT]->GetXaxis()->SetTitle("p^{e}_{T} [GeV/c]");
    else
      hratios[iT]->GetXaxis()->SetTitle("p^{#gamma}_{T} [GeV/c]");
    
    hratios[iT]->GetYaxis()->SetTitle("L1 Efficiency");
    hratios[iT]->SetMarkerStyle(marks[iT]);
    hratios[iT]->SetMarkerColor(Colors[iT]);
    hratios[iT]->SetMarkerSize(2);
    hratios[iT]->SetLineColor(Colors[iT]);
    hratios[iT]->SetMinimum(0);
    hratios[iT]->SetMaximum(1.2);
    hratios[iT]->GetXaxis()->SetLimits(ax_min,ax_max);
    
    if(iT==0)
      hratios[0]->Draw("AP");
    else
      hratios[iT]->Draw("P");

    leg.AddEntry(hratios[iT] ,trigger_names[iT].c_str(),"ep");
  }

  leg.SetFillColorAlpha(kWhite,0);
  leg.Draw();

  TLine *l1 = new TLine(ax_min,1,ax_max,1);
  l1->SetLineWidth(2);
  l1->SetLineStyle(2);
  l1->Draw();

  if(isEB){
    TLatex *pt = new TLatex(0.78,0.42+endcap_legendpp,"|#eta^{#gamma}| < 1.44");
    pt->SetTextFont(42);
    pt->SetNDC(true);
    pt->Draw();
  }

  if(isEC){
    TLatex *pt = new TLatex(0.78,0.42+endcap_legendpp,"|#eta^{#gamma}| > 1.44");
    pt->SetTextFont(42);
    pt->SetNDC(true);
    pt->Draw();
  }

  CMS_lumi( canvas1, iPeriod, iPos );

  gSystem->Exec(("mkdir -p " + folder).c_str());
  canvas1->SaveAs( (folder + "/" + filename + ".png").c_str()); //_cent10 _EB _cent30_80

  return 0; 
}

std::vector<int> GetPrimaryColors()
{
   static std::vector<int> Colors;
   if(Colors.size() > 0)
      return Colors;
   //Colors.push_back(TColor::GetColor("#E74C3C")); // Alizarin
   //Colors.push_back(TColor::GetColor("#3498DB")); // Peter River
   //Colors.push_back(TColor::GetColor("#F1C40F")); // Sum Flower
   //Colors.push_back(TColor::GetColor("#2ECC71")); // Emerald
   //Colors.push_back(TColor::GetColor("#7F8C8D")); // Asbestos
   //Colors.push_back(TColor::GetColor("#8E44AD")); // Wisteria
   //Colors.push_back(TColor::GetColor("#2C3E50")); // Green Sea (dark)
   //Colors.push_back(TColor::GetColor("#16A085")); // Green Sea
   //Colors.push_back(TColor::GetColor("#E67E22")); // Carrot

   Colors.push_back(TColor::GetColor("#377eb8")); // (Blue)
   Colors.push_back(TColor::GetColor("#ff7f00")); // (Orange)
   Colors.push_back(TColor::GetColor("#4daf4a")); // (Green)
   Colors.push_back(TColor::GetColor("#a65628")); // (Brown)
   Colors.push_back(TColor::GetColor("#984ea3")); // (Purple)
   Colors.push_back(TColor::GetColor("#e41a1c")); // (Red)
   Colors.push_back(TColor::GetColor("#f781bf")); // (Pink)

   return Colors;
}

std::vector<std::string> my_glob(const char *pattern) {
    glob_t g;
    glob(pattern, GLOB_TILDE, nullptr, &g); // one should ensure glob returns 0!
    std::vector<std::string> filelist;
    filelist.reserve(g.gl_pathc);
    for (size_t i = 0; i < g.gl_pathc; ++i) {
        filelist.emplace_back(g.gl_pathv[i]);
    }
    globfree(&g);
    return filelist;
}