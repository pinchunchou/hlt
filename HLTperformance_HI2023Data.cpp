#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/CorrelationTuple/EventMatcher.h"
#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/TreeHeaders/ggHiNtuplizerTree.h"
#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/TreeHeaders/hltObjectTree.h"
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
#include "HI_EGMCut2023.h"

#include "/afs/cern.ch/user/p/pchou/PhysicsHIZHadron2022/CommonCode/include/CommandLine.h"

using namespace std;

std::vector<int> GetPrimaryColors();
std::vector<std::string> my_glob(const char *pattern);

//void HLTperformance_HI2023Data(int trigType)
int main(int argc, char *argv[])
{

  CommandLine CL(argc, argv);

  int trigType    = CL.GetInt("trigType", 0);
  int idCut       = CL.GetInt("idCut", 0);
  bool isMC       = CL.GetBool("isMC", true);
  bool isPbPb     = CL.GetBool("isPbPb", true);
  bool noL1       = CL.GetBool("noL1", false);

  bool isHLTObj   = CL.GetBool("isHLTObj", false);

  bool nocut      = CL.GetBool("nocut", false);

  bool isL1denom  = CL.GetBool("isL1denom", false);
  
  string folder   = CL.Get("folder", "figs/20230926/");
  //string dataDir  = CL.Get("dataDir", "/eos/cms/store/group/phys_heavyions/wangj/RECO2023/miniaod_HIExpressRawPrime_374322/");
  string filename = CL.Get("filename", "HLTEff_");
  string runtext  = CL.Get("runtext", "Run #### data");
  string suffix   = CL.Get("suffix", "");
  string trigsuf  = CL.Get("trigsuf", "");
  double dRmax    = CL.GetDouble("dRmax", 0.5);
  double GendRmax = CL.GetDouble("GendRmax", 0.1);
  double Fraction = CL.GetDouble("Fraction", 1);

  vector<string> dataDir = CL.GetStringVector("dataDir",vector<string>{});

  vector<double> PSvec = CL.GetDoubleVector("PSvec",vector<double>{1,1,1});
  vector<double> LPSvec = CL.GetDoubleVector("LPSvec",vector<double>{1,1,1});
  vector<int> L1ID = CL.GetIntVector("L1ID",vector<int>{4,4,5});

  vector<string> triggers;      
  vector<string> triggerobjs;
  vector<string> trigger_names; 

  string l1list[7] = {
    "L1_ZeroBias", //0
    "L1_MinimumBiasHF1_AND_BptxAND", //1
    "L1_SingleEG7_BptxAND", //2
    "L1_SingleEG15_BptxAND", //3
    "L1_SingleEG21_BptxAND", //4
    "L1_SingleEG30_BptxAND", //5
    "L1_DoubleEG5_BptxAND" //6
  };

  if(trigType==-1){
    triggers = CL.GetStringVector("triggers",vector<string>{"HLT_HIGEDPhoton40_v8","HLT_HIGEDPhoton50_v8","HLT_HIGEDPhoton60_v8"});
    trigger_names = CL.GetStringVector("trigger_names",vector<string>{"HI GED Photon 40","HI GED Photon 50","HI GED Photon 60"});
  }

  double endcap_legendpp       = CL.GetDouble("endcap_legendpp",0);
  double doublee_legshift      = CL.GetDouble("doublee_legshift",0);

  bool isEB           = CL.GetBool("isEB", false);
  bool isEC           = CL.GetBool("isEC", false);
  bool isEle          = CL.GetBool("isEle", false);
  bool isDoubleEle    = CL.GetBool("isDoubleEle", false);
  bool isMass50       = CL.GetBool("isMass50", false);
  bool isDoublePhoton = CL.GetBool("isDoublePhoton", false);

  int pt_min    = CL.GetInt("pt_min", 10);
  int pt_max    = CL.GetInt("pt_max", 90);
  int pt_bin    = CL.GetInt("pt_bin", 28);
  int ax_min    = CL.GetInt("ax_min", 10);
  int ax_max    = CL.GetInt("ax_max",105);

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
  int marks[10] = {20,21,33,24,25,27,22,23,26,32};

  writeExtraText = true;       // if extra text
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=11;

  string pp_HI;
  if(isPbPb)
    pp_HI = "HI";
  else
    pp_HI = "PPRef";

  if(trigType==0){
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton40_v");
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton50_v");
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton60_v"); 
    trigger_names.push_back(pp_HI +" GED Photon 40");
    trigger_names.push_back(pp_HI +" GED Photon 50");
    trigger_names.push_back(pp_HI +" GED Photon 60");
    filename+="QCDPhoton";
  }else if(trigType==1){
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton10_v");
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton20_v");
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton30_v"); 
    trigger_names.push_back(pp_HI +" GED Photon 10");
    trigger_names.push_back(pp_HI +" GED Photon 20");
    trigger_names.push_back(pp_HI +" GED Photon 30");
    pt_min = 5; pt_max=50; pt_bin=25; ax_min = 5; ax_max=55;
    filename+="QCDPhoton123";
  }else if(trigType==2){
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton40_EB_v");
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton50_EB_v");
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton60_EB_v"); 
    trigger_names.push_back(pp_HI +" GED Photon 40 EB");
    trigger_names.push_back(pp_HI +" GED Photon 50 EB");
    trigger_names.push_back(pp_HI +" GED Photon 60 EB");
    isEB = true;
    filename+="QCDPhoton_EB";
  }else if(trigType==3){
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton10_EB_v");
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton20_EB_v");
    triggers.push_back("HLT_" + pp_HI + "GEDPhoton30_EB_v"); 
    trigger_names.push_back(pp_HI +" GED Photon 10 EB");
    trigger_names.push_back(pp_HI +" GED Photon 20 EB");
    trigger_names.push_back(pp_HI +" GED Photon 30 EB");
    isEB = true;
    pt_min = 5; pt_max=50; pt_bin=25; ax_min = 5; ax_max=55;
    filename+="QCDPhoton123_EB";
  }else if(trigType==4){
    triggers.push_back("HLT_" + pp_HI + "Ele10Gsf_v");
    triggers.push_back("HLT_" + pp_HI + "Ele15Gsf_v");
    triggers.push_back("HLT_" + pp_HI + "Ele20Gsf_v"); 
    trigger_names.push_back(pp_HI +" Ele 10 Gsf");
    trigger_names.push_back(pp_HI +" Ele 15 Gsf");
    trigger_names.push_back(pp_HI +" Ele 20 Gsf");
    isEle=true;
    pt_min = 5; pt_max=50; pt_bin=25; ax_min = 5; ax_max=55;
    filename+="Zee";
  }else if(trigType==5){
    triggers.push_back("HLT_" + pp_HI + "Ele30Gsf_v");
    triggers.push_back("HLT_" + pp_HI + "Ele40Gsf_v");
    triggers.push_back("HLT_" + pp_HI + "Ele50Gsf_v"); 
    trigger_names.push_back(pp_HI +" Ele 30 Gsf");
    trigger_names.push_back(pp_HI +" Ele 40 Gsf");
    trigger_names.push_back(pp_HI +" Ele 50 Gsf");
    isEle=true;
    pt_min = 5; pt_max=80; pt_bin=40; ax_min = 5; ax_max=85;
    filename+="Zee345";
  }else if(trigType==6){
    triggers.push_back("HLT_" + pp_HI + "DoubleEle10Gsf_v" );
    triggers.push_back("HLT_" + pp_HI + "Ele15Ele10Gsf_v");
    triggers.push_back("HLT_" + pp_HI + "DoubleEle15Gsf_v"); 
    trigger_names.push_back(pp_HI +" Double Ele 10 Gsf");
    trigger_names.push_back(pp_HI +" Ele 15 Ele 10 Gsf");
    trigger_names.push_back(pp_HI +" Double Ele 15 Gsf");
    isEle=true;
    isDoubleEle=true;
    doublee_legshift=-0.1;
    pt_min = 5; pt_max=50; pt_bin=25; ax_min = 5; ax_max=55;
    filename+="doublee";
  }else if(trigType==7){
    triggers.push_back("HLT_" + pp_HI + "DoubleEle10GsfMass50_v");
    triggers.push_back("HLT_" + pp_HI + "Ele15Ele10GsfMass50_v");
    triggers.push_back("HLT_" + pp_HI + "DoubleEle15GsfMass50_v"); 
    trigger_names.push_back(pp_HI +" Double Ele 10 Gsf Mass50");
    trigger_names.push_back(pp_HI +" Ele 15 Ele 10 Gsf Mass50");
    trigger_names.push_back(pp_HI +" Double Ele 15 Gsf Mass50");
    isEle=true;
    isDoubleEle=true;
    isMass50=true;
    doublee_legshift=-0.2;
    pt_min = 5; pt_max=50; pt_bin=25; ax_min = 5; ax_max=55;
    filename+="doublee_mass50";
  }else if(trigType==8){
    triggers.push_back("HLT_" + pp_HI + "DoubleGEDPhoton20_v"); 
    trigger_names.push_back(pp_HI +" Double GED Photon 20");
    isDoublePhoton = true;
    doublee_legshift=-0.1;
    pt_min = 5; pt_max=90; pt_bin=17; ax_min = 5; ax_max=100;
    filename+="diphoton";
  }

  if(noL1&&isEle==false&&isMC)
    filename+="_noL1";

  if(isEB)
    filename+="_barrel";

  if(isEC)
    filename+="_endcap";

  filename+=suffix;

  std::cout << "Reading files..." << std::endl;

  // For loop over file list
  //TChain *HltTree = new TChain("hltanalysis/HltTree");  
  //TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  //TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");

  TTree* HltTree;
  TTree* EventTree;
  TTree* HiTree;

  
  vector <TH1F *> hnums;
  vector <TH1F *> hdenoms;
  vector <TTree *> Trees;

  Int_t counts[100]={0};
  Int_t count0 = 0;
  double PSs[100]={0};
  Int_t N_F[4]={0}, N_M[4]={0};

  Int_t entriesTot = 0;

  for(long unsigned int iT = 0; iT < triggers.size(); iT++){
    hdenoms.push_back(new TH1F(("hdenum"+ to_string(iT)).c_str(),"", pt_bin, pt_min, pt_max));
    hnums.push_back(new TH1F(("hnum"+ to_string(iT)).c_str(),"", pt_bin, pt_min, pt_max));
  }

  cout<<"dataDir.size() = "<<dataDir.size()<<endl;
  for(int iF = 0; iF < dataDir.size(); iF++){
    //cout<<"iF = "<<iF<<endl;
    int file_count = 0;
    cout<<"dataDir["<<iF<<"] = "<<dataDir[iF]<<endl;
    for (const auto &filename : my_glob(dataDir[iF].c_str())){
      Trees.clear();
      //cout<<"file_count = "<<file_count<<endl;
      //cout<<"filename = "<<filename<<endl;
      file_count++;

      TFile file_data(filename.c_str(),"READ");
      TTree* HltTree = (TTree*) file_data.Get("hltanalysis/HltTree");
      TTree* EventTree = (TTree*) file_data.Get("ggHiNtuplizer/EventTree");
      TTree* HiTree = (TTree*) file_data.Get("hiEvtAnalyzer/HiTree");

      //HltTree->Add(filename.c_str());
      //EventTree->Add(filename.c_str());
      //HiTree->Add(filename.c_str());

      //cout<<"A0"<<endl;

      TObjArray *toa = HltTree->GetListOfBranches();

      HltTree->SetBranchStatus("*",0);
      HltTree->SetBranchStatus("Event", 1);
      HltTree->SetBranchStatus("LumiBlock", 1);
      HltTree->SetBranchStatus("Run", 1);

      //cout<<"A1"<<endl;

      //Bool_t HLT_Trig[100];
      Int_t HLT_Trig[100], HLT_PSNum[100], HLT_PSDnm[100], L1_Trig[100], L1_PS[100];

      for(int iP = 0;iP<size(l1list); iP++){
        HltTree->SetBranchStatus((l1list[iP]).c_str(), 1);
        HltTree->SetBranchAddress((l1list[iP]).c_str(), &L1_Trig[iP]);
        HltTree->SetBranchStatus((l1list[iP]+"_Prescl").c_str(), 1);
        HltTree->SetBranchAddress((l1list[iP]+"_Prescl").c_str(), &L1_PS[iP]);
      }

      //cout<<"A2"<<endl;

      for(long unsigned int iT = 0; iT < triggers.size(); iT++){
        //hdenoms.push_back(new TH1F(("hdenum"+ to_string(iT)).c_str(),"", pt_bin, pt_min, pt_max));
        //hnums.push_back(new TH1F(("hnum"+ to_string(iT)).c_str(),"", pt_bin, pt_min, pt_max));

        //cout<<"AA1"<<endl;
        if(isHLTObj){
          Trees.push_back((TTree*) file_data.Get(("hltobject/" + triggers[iT]).c_str()));
          //cout<<"AA2"<<endl;
          //cout<<"triggers["<<iT<<"] = "<<triggers[iT]<<endl;
          //cout<<"hltobject entries = "<<Trees[iT]->GetEntries()<<endl;
          Trees[iT]->SetBranchStatus("*",0);
          Trees[iT]->SetBranchStatus("eta", 1);
          Trees[iT]->SetBranchStatus("phi", 1);
          Trees[iT]->SetBranchStatus("pt", 1);
          //cout<<"AA3"<<endl;
        }

        //cout<<"A3"<<endl;

        for(int trig_ver = 0; trig_ver<100; trig_ver++){
          if(toa->Contains((triggers[iT] + to_string(trig_ver)).c_str())){
            trigsuf = to_string(trig_ver);
            break;
          }
        }

        //std::cout<<"triggers["<<iT<<"] = "<<triggers[iT]<<trigsuf<<std::endl;
        HltTree->SetBranchStatus((triggers[iT] + trigsuf).c_str(), 1);
        HltTree->SetBranchAddress((triggers[iT] + trigsuf).c_str(), &HLT_Trig[iT]);    

        //cout<<"a"<<endl;

        if(isMC==false){
          HltTree->SetBranchStatus(((triggers[iT] + trigsuf)+"_PrescaleNumerator").c_str(), 1);
          HltTree->SetBranchAddress(((triggers[iT] + trigsuf)+"_PrescaleNumerator").c_str(), &HLT_PSNum[iT]);

          HltTree->SetBranchStatus(((triggers[iT] + trigsuf)+"_PrescaleDenominator").c_str(), 1);
          HltTree->SetBranchAddress(((triggers[iT] + trigsuf)+"_PrescaleDenominator").c_str(), &HLT_PSDnm[iT]);
        }

        //cout<<"b"<<endl;
      }

      //cout<<"a1"<<endl;

      HiTree->SetBranchStatus("*",0);
      HiTree->SetBranchStatus("hiBin",1);
      HiTree->SetBranchStatus("hiHF",1);
      HiTree->SetBranchStatus("vz",1);
      if(isMC)
        HiTree->SetBranchStatus("weight",1);

      ULong64_t       hlt_event;
      Int_t           hlt_lumi, hlt_run, hiBin;
      Float_t         hiHF, pthat_weight, vz;

      HiTree->SetBranchAddress("hiBin", &hiBin);
      HiTree->SetBranchAddress("hiHF", &hiHF);
      HiTree->SetBranchAddress("vz", &vz);

      if(isMC)
        HiTree->SetBranchAddress("weight", &pthat_weight);

      HltTree->SetBranchAddress("Event", &hlt_event);
      HltTree->SetBranchAddress("LumiBlock", &hlt_lumi);
      HltTree->SetBranchAddress("Run", &hlt_run);

      EventTree->SetBranchStatus("*",0);
      EventTree->SetBranchStatus("run",1);
      EventTree->SetBranchStatus("event",1);
      EventTree->SetBranchStatus("lumis",1);

      //cout<<"a2"<<endl;

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
          EventTree->SetBranchStatus("phoR9",1);

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

      //cout<<"a3"<<endl;

      EventMatcher* em = new EventMatcher();

      //std::cout << "Loop over the file with offline reco output..." << std::endl;

      ggHiNtuplizer ggHi;
      ggHi.setupTreeForReading(EventTree);

      Long64_t entriesTmp = EventTree->GetEntries();
      //std::cout << "entries in File = " << entriesTmp << std::endl;

      Long64_t entryTrig = 0;
      Int_t entriesNotFoundinTrigger=0;
      Int_t duplicateEntries=0, entriesAnalyzed=0;


      for (Long64_t j_entry = 0; j_entry < entriesTmp*Fraction; ++j_entry){
        if (entriesTot % 10000 == 0)  {
          cout << "current entry = " <<entriesTot<<endl;
          //std::cout << "current entry = " <<j_entry<< " out of " <<entriesTmp<< " : " <<std::setprecision(2)<<(double)j_entry/entriesTmp*100<< " %" << std::endl;
        }
        entriesTot++;

        EventTree->GetEntry(j_entry);
        HiTree->GetEntry(j_entry);
        
        bool eventAdded = em->addEvent(ggHi.run, ggHi.lumis, ggHi.event, j_entry);
        if(!eventAdded) // this event is duplicate, skip this one.
        {
          duplicateEntries++;
          continue;
        }

        entryTrig = j_entry;

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


        float maxPt=0, max2Pt=0;
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
          //cout<<"ele matchedIndex"<<endl;
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
            //cout<<"ele matchedIndex push_back"<<endl;
          }
          //cout<<"ele matchedIndex done"<<endl;
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

        //std::cout<<"Photon Selection"<<std::endl;

        bool sel_cut = false, sel_cut2 = false;

        sel_cut = PassCuts(ggHi, maxPt_i, idCut, isEC, isEB, isPbPb, isEle, hiBin);

        if(isDoublePhoton||isDoubleEle)
          sel_cut2 = PassCuts(ggHi, max2Pt_i, idCut, isEC, isEB, isPbPb, isEle, hiBin);
        else
          sel_cut2 = true;

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

        //cout<<"sel_cut = "<<sel_cut<<", sel_cut2 = "<<sel_cut2<<endl;

        float Et, EG_eta;
        if(isEle){
          Et = (*ggHi.elePt)[maxPt_i];
          EG_eta = (*ggHi.eleEta)[maxPt_i];
        }else{
          Et = (*ggHi.phoEt)[maxPt_i];
          EG_eta = (*ggHi.phoSCEta)[maxPt_i];
        }


        if(isEle==false&&sel_cut&&sel_cut2&&Et>50&&HLT_Trig[0]==0&&abs(EG_eta) < 1.442){
          cout<<"Run:"<<hlt_run<<", LumiBlock: "<<hlt_lumi<<", Event: "<<hlt_event<<endl;

          for(long unsigned int iT = 0; iT < triggers.size(); iT++)
            std::cout<<triggers[iT]<<trigsuf<<" = "<<HLT_Trig[iT]<<", ";

          cout<<endl;

          for(int iL=3;iL<6;iL++)
            cout<<l1list[iL]<<" = "<< L1_Trig[iL] <<", ";

          cout<<endl;

          for(int iC=0;iC<4;iC++)
            cout<<"Pass "<<iC<<" cut: "<<PassCuts(ggHi, maxPt_i, iC, isEC, isEB, isPbPb, isEle, hiBin)<<", ";

          cout<<endl;

          std::cout<<"phoEt = "<<(*ggHi.phoEt)[maxPt_i]<<", phoEta = "<<(*ggHi.phoEta)[maxPt_i]<<", phoPhi = "<<(*ggHi.phoPhi)[maxPt_i]<<std::endl;
          cout<<"phoHoverE = "<<(*ggHi.phoHoverE)[maxPt_i]<<", phoSigmaIEtaIEta_2012 = "<<(*ggHi.phoSigmaIEtaIEta_2012)[maxPt_i]<<endl;
          cout<<"phoR9 = "<<(*ggHi.phoR9)[maxPt_i]<<", pho_swissCrx = "<<(*ggHi.pho_swissCrx)[maxPt_i]<<endl;
          cout<<"hiHF = "<<hiHF<<", hiBin = "<<hiBin<<", vz = "<<vz<<endl;

          cout<<endl;

          cout<<"Other photon/electrons:"<<endl;
          for(int i=0; i < num_of_PE; ++i) {
            std::cout<<"phoEt = "<<(*ggHi.phoEt)[i]<<", phoEta = "<<(*ggHi.phoEta)[i]<<", phoPhi = "<<(*ggHi.phoPhi)[i]<<std::endl;
            cout<<"phoHoverE = "<<(*ggHi.phoHoverE)[i]<<", phoSigmaIEtaIEta_2012 = "<<(*ggHi.phoSigmaIEtaIEta_2012)[i]<<endl;
            cout<<"phoR9 = "<<(*ggHi.phoR9)[i]<<", pho_swissCrx = "<<(*ggHi.pho_swissCrx)[i]<<endl;
          }

          cout<<"------------"<<endl;
        }


        for(int iC=0;iC<4;iC++){
      
      //Is Barrel

          //cout<<"iC = "<<iC<<endl;

          //cout<<"EG_eta = "<<EG_eta<<endl;
          //cout<<"Et = "<<Et<<endl;
          if(abs(EG_eta) >= 1.442||Et<50)
            continue;

          //cout<<"maxPt_i = "<<maxPt_i<<", isEC = "<<isEC<<", isEB = "<<isEB<<", isPbPb = "<<isPbPb<<", isEle = "<<isEle<<", hiBin = "<<hiBin<<endl;

          bool pass_or_not = PassCuts(ggHi, maxPt_i, iC, isEC, isEB, isPbPb, isEle, hiBin);

          //cout<<"Pass "<<iC<<" cut: "<<pass_or_not<<", ";

          if(!pass_or_not)
            continue;

          if(HLT_Trig[0]>0)
            N_F[iC]++;
          else
            N_M[iC]++;
        }

        for(long unsigned int iT = 0; iT < triggers.size(); iT++){
          hltObject hObj;
          if(isHLTObj){ 
            hObj.setupTreeForReading(Trees[iT]);
            Trees[iT]->GetEntry(entryTrig);
          }
          //std::cout<<"c"<<std::endl;

          float drTrig;

          PSs[iT] = (double) HLT_PSNum[iT]/HLT_PSDnm[iT];

          if(isMC==false&&(PSs[iT]!=PSvec[iT]||L1_PS[L1ID[iT]]!=LPSvec[iT])){
            if(PSs[iT]!=0 && L1_PS[L1ID[iT]]!=0)
              cout<<"PSs["<<iT<<"] = "<<PSs[iT]<<", PSvec["<<iT<<"] = "<<PSvec[iT]<<", L1_PS["<<L1ID[iT]<<"] = "<<L1_PS[L1ID[iT]]<<", LPSvec["<<iT<<"] = "<<LPSvec[iT]<<endl;
            continue;
          }

          if(isL1denom&&L1_Trig[L1ID[iT]]==0)
            continue;

          if(sel_cut&&sel_cut2){

            //cout<<"Pass offline cut"<<endl;
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


            if(HLT_Trig[iT]>0){
              if(isHLTObj){
                n_pho_ho = size(*hObj.eta);
                for(int i=0;i<n_pho_ho;i++){
                  if(isEle==false&&isDoublePhoton==false)
                    drtemp = sqrt(((*hObj.eta)[i]-(*ggHi.phoEta)[maxPt_i])*((*hObj.eta)[i]-(*ggHi.phoEta)[maxPt_i])+((*hObj.phi)[i]-(*ggHi.phoPhi)[maxPt_i])*((*hObj.phi)[i]-(*ggHi.phoPhi)[maxPt_i]));
                  else if(isEle&&isDoubleEle==false)
                    drtemp = sqrt(((*hObj.eta)[i]-(*ggHi.eleEta)[maxPt_i])*((*hObj.eta)[i]-(*ggHi.eleEta)[maxPt_i])+((*hObj.phi)[i]-(*ggHi.elePhi)[maxPt_i])*((*hObj.phi)[i]-(*ggHi.elePhi)[maxPt_i]));
                  else if(isDoubleEle)
                    drtemp = sqrt(((*hObj.eta)[i]-(*ggHi.eleEta)[max2Pt_i])*((*hObj.eta)[i]-(*ggHi.eleEta)[max2Pt_i])+((*hObj.phi)[i]-(*ggHi.elePhi)[max2Pt_i])*((*hObj.phi)[i]-(*ggHi.elePhi)[max2Pt_i]));
                  else if(isDoublePhoton)
                    drtemp = sqrt(((*hObj.eta)[i]-(*ggHi.phoEta)[max2Pt_i])*((*hObj.eta)[i]-(*ggHi.phoEta)[max2Pt_i])+((*hObj.phi)[i]-(*ggHi.phoPhi)[max2Pt_i])*((*hObj.phi)[i]-(*ggHi.phoPhi)[max2Pt_i]));

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
              }else{
                std::cout<<"n_pho_ho = "<<n_pho_ho<<std::endl;
                std::cout<<"drTrig = "<<drTrig<<std::endl;
              }
            }
          }
        }
      }
      file_data.Close();
    }
  }


  std::cout << "count0 = "<<count0;
  for(long unsigned int iT = 0; iT < triggers.size(); iT++){
    std::cout << ", counts["<<iT<<"] = "<<counts[iT];
  }

  std::cout << std::endl;

  cout<<"N_F[0] = "<<N_F[0]<<", N_M[0] = "<<N_M[0]<<endl;
  for(int iC=1;iC<4;iC++){
    cout<<"N_F["<<iC<<"] = "<<N_F[iC]<<", N_M["<<iC<<"] = "<<N_M[iC]<<endl;
    cout<<"N_F["<<iC<<"]/N_F[0] = "<<(float) N_F[iC]/N_F[0]<<", N_M["<<iC<<"]/N_M[0] = "<<(float) N_M[iC]/N_M[0]<<endl;
  }

  std::cout << std::endl;

  std::cout << "###" << std::endl;
  std::cout << "Loop Forest ENDED"  << std::endl;
  std::cout << "entries Forest          = " << entriesTot << std::endl;
  //std::cout << "duplicateEntries Forest = " << duplicateEntries << std::endl;
  //std::cout << "entriesAnalyzed Forest  = " << entriesAnalyzed << std::endl;
  //std::cout << "entriesNotFoundinTrigger  = " << entriesNotFoundinTrigger << std::endl;
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
      std::cout<<"Trigger: "<<triggers[iT]<<trigsuf<<std::endl;
      std::cout<<"HLT PS:"<<PSvec[iT]<<std::endl;
      std::cout<<"L1 PS:"<<LPSvec[iT]<<std::endl;
      hratios[iT]->Scale(PSvec[iT]*LPSvec[iT]);
    }

    if(isEle)
      hratios[iT]->GetXaxis()->SetTitle("p^{e}_{T} [GeV/c]");
    else
      hratios[iT]->GetXaxis()->SetTitle("p^{#gamma}_{T} [GeV/c]");

    hratios[iT]->GetYaxis()->SetTitle("HLT Efficiency");
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
   Colors.push_back(TColor::GetColor("#f781bf")); // (Pink)
   Colors.push_back(TColor::GetColor("#a65628")); // (Brown)
   Colors.push_back(TColor::GetColor("#984ea3")); // (Purple)
   Colors.push_back(TColor::GetColor("#e41a1c")); // (Red)

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