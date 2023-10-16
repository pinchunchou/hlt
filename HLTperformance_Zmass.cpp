#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/CorrelationTuple/EventMatcher.h"
#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/TreeHeaders/ggHiNtuplizerTree.h"
#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/TreeHeaders/hltObjectTree.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"

#include <glob.h>
void style(){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetTitleFont(43);
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetStatFont(43);

  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetOptStat(0); /*don't show statistics box*/
  gStyle->SetOptTitle(0); /*don't show histogram titles*/
  gStyle->SetTitleSize(48, "xyz");
  gStyle->SetTitleOffset(1, "xyz");
  gStyle->SetLabelSize(36, "xyz");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  gStyle->SetLineScalePS(1.5);

  gROOT->ForceStyle();
}
std::vector<std::string> my_glob(const char *pattern);
void HLTperformance_Zmass(){

  //style();
  setTDRStyle();
  gStyle->SetLegendBorderSize(0);

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  int iPos=11;

  lumi_sqrtS = "Run 374810 (5.36 TeV)";

  std::cout << "Reading files..." << std::endl;

  //string frtfile = "/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime[3,4,5,6,8,9]/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime[3,4,5,6,8,9]_374354_Dpt2trk1/*/0000/*.root";
  string frtfile = "/eos/cms/store/group/phys_heavyions/mstojano/run3RapidValidation/PbPb2023_run374810_HIPhysicsRawPrime0_Forest_2023-10-09/CRAB_UserFiles/crab_PbPb2023_run374810_HIPhysicsRawPrime0_Forest_2023-10-09/231009_162633/000*/*.root";

  //1320V33
  string hltfile = frtfile;

  TChain *HltTree = new TChain("hltanalysis/HltTree");

  for (const auto &filename : my_glob(hltfile.c_str()))
    HltTree->Add(filename.c_str());

  TChain *EventTree = new TChain("ggHiNtuplizer/EventTree");
  TChain *HiTree = new TChain("hiEvtAnalyzer/HiTree");

  for (const auto &filename : my_glob(frtfile.c_str())){
    //cout<<"filename = "<<filename<<endl;
    EventTree->Add(filename.c_str());
    HiTree->Add(filename.c_str());
  }
  

  int pt_min = 60, pt_max=120, pt_bin=30;
  int ax_min = 60, ax_max=120;

  TH1F* h_Zmass10   = new TH1F("h_Zmass10","",pt_bin,pt_min,pt_max);
  TH1F* h_Zmass15   = new TH1F("h_Zmass15","",pt_bin,pt_min,pt_max);
  TH1F* h_Zmass20   = new TH1F("h_Zmass20","",pt_bin,pt_min,pt_max);

  HltTree->SetBranchStatus("*",0);
  HltTree->SetBranchStatus("Event", 1);
  HltTree->SetBranchStatus("LumiBlock", 1);
  HltTree->SetBranchStatus("Run", 1);

  HltTree->SetBranchStatus("HLT_HIEle10Gsf_v10", 1);
  HltTree->SetBranchStatus("HLT_HIEle15Gsf_v10", 1);
  HltTree->SetBranchStatus("HLT_HIEle20Gsf_v10", 1);

  HltTree->SetBranchStatus("HLT_HIEle30Gsf_v10", 1);
  HltTree->SetBranchStatus("HLT_HIEle40Gsf_v10", 1);
  HltTree->SetBranchStatus("HLT_HIEle50Gsf_v10", 1);


  HltTree->SetBranchStatus("HLT_HIGEDPhoton10_v10", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton20_v10", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton30_v10", 1);

  HltTree->SetBranchStatus("HLT_HIGEDPhoton40_v10", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton50_v10", 1);
  HltTree->SetBranchStatus("HLT_HIGEDPhoton60_v10", 1);

  HiTree->SetBranchStatus("*",0);
  HiTree->SetBranchStatus("hiBin",1);
  HiTree->SetBranchStatus("hiHF",1);

  ULong64_t       hlt_event;
  Int_t           hlt_lumi, hlt_run, hiBin;
  Float_t         hiHF;

  Int_t HLT_HIEle10Gsf, HLT_HIEle15Gsf, HLT_HIEle20Gsf;

  HiTree->SetBranchAddress("hiBin", &hiBin);
  HiTree->SetBranchAddress("hiHF", &hiHF);

  HltTree->SetBranchAddress("Event", &hlt_event);
  HltTree->SetBranchAddress("LumiBlock", &hlt_lumi);
  HltTree->SetBranchAddress("Run", &hlt_run);

  HltTree->SetBranchAddress("HLT_HIGEDPhoton10_v10", &HLT_HIEle10Gsf);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton20_v10", &HLT_HIEle15Gsf);
  HltTree->SetBranchAddress("HLT_HIGEDPhoton30_v10", &HLT_HIEle20Gsf);

  EventTree->SetBranchStatus("*",0);
  EventTree->SetBranchStatus("run",1);
  EventTree->SetBranchStatus("event",1);
  EventTree->SetBranchStatus("lumis",1);

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

  //EventTree->SetBranchStatus("mcEta",1);
  //EventTree->SetBranchStatus("mcPhi",1);

  EventTree->SetBranchStatus("eleSCRawEn",1);


  Long64_t entriesHlt = HltTree->GetEntries();
  std::cout << "entries in HLT = " << entriesHlt << std::endl;


  std::cout << "Loop over the file with offline reco output..." << std::endl;

  ggHiNtuplizer ggHi;
  ggHi.setupTreeForReading(EventTree);



  Long64_t entriesTmp = EventTree->GetEntries();
  std::cout << "entries in File = " << entriesTmp << std::endl;

  Long64_t entryTrig = 0;
  Int_t entriesNotFoundinTrigger=0;
  Int_t duplicateEntries=0, entriesAnalyzed=0;

  Int_t count0=0, count10=0, count15=0, count20=0;
  Int_t count0_0=0, count0_00=0;

  for (Long64_t j_entry = 0; j_entry < entriesTmp; ++j_entry){
    if (j_entry % 10000 == 0)  {
        std::cout << "current entry = " <<j_entry<< " out of " <<entriesTmp<< " : " <<std::setprecision(2)<<(double)j_entry/entriesTmp*100<< " %" << std::endl;
    }

    EventTree->GetEntry(j_entry);
    HiTree->GetEntry(j_entry);


    entriesAnalyzed++;

    //cout<<"entryTrig = "<<entryTrig<<endl;
    HltTree->GetEntry(j_entry);

    float sumIso, maxPt=0, max2Pt=0;
    int maxPt_i = -1, max2Pt_i = -1;

    float dRmax = 0.1;
    float GendRmax = 0.1;

    float dr00;
    for(int i=0; i < ggHi.nEle; ++i) {

      if((*ggHi.elePt)[i]>maxPt){
        maxPt = (*ggHi.elePt)[i]; 
        maxPt_i = i;
      }else if((*ggHi.elePt)[i]>max2Pt&&(*ggHi.elePt)[i]<maxPt){
        max2Pt = (*ggHi.elePt)[i]; 
        max2Pt_i = i;
      }
    }

    if(maxPt_i==-1) continue;
    count0_0++;
    if(max2Pt_i==-1)  continue;


    float dr10, dr15, dr20;

     bool sel_cut = kTRUE, sel_cut2 = kTRUE;
    /*
    if(hiBin>=0 && hiBin < 60){
      if(abs((*ggHi.eleEta)[maxPt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0135 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0038 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.0376 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1616 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0177  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = kTRUE;
      else if(abs((*ggHi.eleEta)[maxPt_i]) > 1.57 && abs((*ggHi.eleEta)[maxPt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0466 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0063 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.1186 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1317 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0201  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = kTRUE;
      else sel_cut = kFALSE;
    }else if(hiBin>=60 && hiBin < 200){
      if(abs((*ggHi.eleEta)[maxPt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0107 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0035 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.0327 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.1268 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0774  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = kTRUE;
      else if(abs((*ggHi.eleEta)[maxPt_i]) > 1.57 && abs((*ggHi.eleEta)[maxPt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[maxPt_i] < 0.0339 && abs((*ggHi.eledEtaSeedAtVtx)[maxPt_i]) < 0.0067 && abs((*ggHi.eledPhiAtVtx)[maxPt_i]) < 0.0838 && (*ggHi.eleHoverEBc)[maxPt_i] < 0.0977 && abs((*ggHi.eleEoverPInv)[maxPt_i]) < 0.0193  && (*ggHi.eleMissHits)[maxPt_i] <= 1 && (*ggHi.eleIP3D)[maxPt_i] < 0.03) sel_cut = kTRUE;
      else sel_cut = kFALSE;
    }else{sel_cut = kFALSE;}


    if(hiBin>=0 && hiBin < 60){
      if(abs((*ggHi.eleEta)[max2Pt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.0135 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.0038 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.0376 && (*ggHi.eleHoverEBc)[max2Pt_i] < 0.1616 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.0177  && (*ggHi.eleMissHits)[max2Pt_i] <= 1 && (*ggHi.eleIP3D)[max2Pt_i] < 0.03) sel_cut2 = kTRUE;
      else if(abs((*ggHi.eleEta)[max2Pt_i]) > 1.57 && abs((*ggHi.eleEta)[max2Pt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.0466 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.0063 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.1186 && (*ggHi.eleHoverEBc)[max2Pt_i] < 0.1317 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.0201  && (*ggHi.eleMissHits)[max2Pt_i] <= 1 && (*ggHi.eleIP3D)[max2Pt_i] < 0.03) sel_cut2 = kTRUE;
      else sel_cut2 = kFALSE;
    }else if(hiBin>=60 && hiBin < 200){
      if(abs((*ggHi.eleEta)[max2Pt_i]) < 1.442 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.0107 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.0035 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.0327 && (*ggHi.eleHoverEBc)[max2Pt_i] < 0.1268 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.0774  && (*ggHi.eleMissHits)[max2Pt_i] <= 1 && (*ggHi.eleIP3D)[max2Pt_i] < 0.03) sel_cut2 = kTRUE;
      else if(abs((*ggHi.eleEta)[max2Pt_i]) > 1.57 && abs((*ggHi.eleEta)[max2Pt_i]) < 2.4 && (*ggHi.eleSigmaIEtaIEta_2012)[max2Pt_i] < 0.0339 && abs((*ggHi.eledEtaSeedAtVtx)[max2Pt_i]) < 0.0067 && abs((*ggHi.eledPhiAtVtx)[max2Pt_i]) < 0.0838 && (*ggHi.eleHoverEBc)[max2Pt_i] < 0.0977 && abs((*ggHi.eleEoverPInv)[max2Pt_i]) < 0.0193  && (*ggHi.eleMissHits)[max2Pt_i] <= 1 && (*ggHi.eleIP3D)[max2Pt_i] < 0.03) sel_cut2 = kTRUE;
      else sel_cut2 = kFALSE;
    }else{sel_cut2 = kFALSE;}
*/
    double emass = 0.0005111;
    TLorentzVector ev1, ev2;  
    ev1.SetPtEtaPhiM((*ggHi.elePt)[maxPt_i] ,(*ggHi.eleEta)[maxPt_i] ,(*ggHi.elePhi)[maxPt_i] ,emass);
    ev2.SetPtEtaPhiM((*ggHi.elePt)[max2Pt_i],(*ggHi.eleEta)[max2Pt_i],(*ggHi.elePhi)[max2Pt_i],emass);

    TLorentzVector Zv12 = ev1+ev2;
    double Zmass = Zv12.M();

    if(Zmass<60 || Zmass>120) sel_cut = kFALSE;

    count0_00++;

    if(sel_cut&&sel_cut2){
      count0++;

      if(HLT_HIEle10Gsf>0){h_Zmass10->Fill(Zmass); count10++;}
      if(HLT_HIEle15Gsf>0){h_Zmass15->Fill(Zmass); count15++;}
      if(HLT_HIEle20Gsf>0){h_Zmass20->Fill(Zmass); count20++;}
    }

  }

  std::cout << "count0 = "<<count0<<", count10 = "<<count10<<", count15 = "<<count15<<", count20 = "<<count20 << std::endl;
  std::cout << "count0_0 = "<<count0_0<<", count0_00 = "<<count0_00<< std::endl;

  std::cout << "###" << std::endl;
  std::cout << "Loop Forest ENDED"  << std::endl;
  std::cout << "entries Forest          = " << entriesTmp << std::endl;
  std::cout << "entriesAnalyzed Forest  = " << entriesAnalyzed << std::endl;
  std::cout << "###" << std::endl;

  TCanvas *canvas1 = new TCanvas("canvas1","Comparison of two histograms",1200,800);
  canvas1->SetMargin(0.15, 0.05, 0.15, 0.08);
  canvas1->GetFrame()->SetBorderSize(12);


  h_Zmass10->SetName("hratio10");
  h_Zmass15->SetName("hratio15");
  h_Zmass20->SetName("hratio20");


  h_Zmass10->GetXaxis()->SetTitle("p^{e}_{T} [GeV/c]");
  h_Zmass10->GetYaxis()->SetTitle("HLT Efficiency");
  h_Zmass10->SetMarkerStyle(kFullCircle);
  h_Zmass10->SetMarkerColor(kBlack);
  h_Zmass10->SetMarkerSize(2);
  h_Zmass10->SetLineColor(kBlack);

  h_Zmass15->SetMarkerStyle(kFullTriangleUp);
  h_Zmass15->SetMarkerColor(kBlue);
  h_Zmass15->SetMarkerSize(2);
  h_Zmass15->SetLineColor(kBlue);

  h_Zmass20->SetMarkerStyle(kFullSquare);
  h_Zmass20->SetMarkerColor(kRed);
  h_Zmass20->SetMarkerSize(2);
  h_Zmass20->SetLineColor(kRed);

  h_Zmass10->SetMinimum(0);
  h_Zmass10->SetMaximum(1.2);

  h_Zmass10->GetXaxis()->SetLimits(ax_min,ax_max);

  h_Zmass10->Draw("P");
  h_Zmass15->Draw("P same");
  h_Zmass20->Draw("P same");

  float endcap_legendpp = 0;//0.3;

  //HLT_PPRefEle10Gsf_v

  TLegend leg(0.58,0.18+endcap_legendpp,0.98,0.4+endcap_legendpp);
  leg.AddEntry(h_Zmass10 ,"HLT_HIGEDPhoton10","ep");
  leg.AddEntry(h_Zmass15 ,"HLT_HIGEDPhoton20","ep");
  leg.AddEntry(h_Zmass20 ,"HLT_HIGEDPhoton30","ep");


  leg.SetFillColorAlpha(kWhite,0);
  leg.Draw();

  TLine *l1 = new TLine(ax_min,1,ax_max,1);
  l1->SetLineWidth(2);
  l1->SetLineStyle(2);
  l1->Draw();

  CMS_lumi( canvas1, iPeriod, iPos );

  gSystem->Exec("mkdir -p figs/20231012");

  //canvas1->SaveAs("figs/20230913/HLTEff_Zee_PbPb_1320V33.pdf"); //
  canvas1->SaveAs("figs/20231012/HLTEff_ZMass_374810_Photon123.png"); //_gt2022

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