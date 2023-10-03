#include "/eos/cms/store/group/phys_heavyions_ops/pchou/ElectroWeak-Jet-Track-Analyses/TreeHeaders/ggHiNtuplizerTree.h"

bool PassCuts(ggHiNtuplizer& ggHi, int idx, int idCut = 0, bool isEC = false, bool isEB = false, bool isPbPb = true, bool isEle = false, int hiBin=-1);

bool PassPhotonCuts_PbPb(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB);
bool PassElectronCuts_PbPb(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB, int hiBin);
bool PassPhotonCuts_pp(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB);
bool PassElectronCuts_pp(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB);


// The cut values can refer to 
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiEgamma2023

//2018 PbPb photon cut values from Kaya
// https://indico.cern.ch/event/787343/contributions/3522682/attachments/1891055/3118657/190807_photonID_and_regression_2018pbpb_2017pp_status.pdf
// See https://twiki.cern.ch/twiki/bin/view/CMS/HiEgamma2023#Information_for_the_photon_branc 
// for the meanings of each variables.


// We use "Set 1: HI iso" here.
//loose (90%), medium (80%), tight (70%), extra tight (60%)
// Barrel region
float phoHoverECut_Barrel_PbPb[4]        = {0.247995,0.238094,0.164101,0.119947};
float phoSigmaIEtaIEtaCut_Barrel_PbPb[4] = {0.012186,0.011024,0.010784,0.010392};
float sumIsoCut_Barrel_PbPb[4]           = {11.697505,6.787425,3.509457,2.099277};

float pho_swissCrxCut                    = 0.9;
float pho_seedTimeCut                    = 3;
// The above two variables are used to reject spikes in data.
// Their cut values can be found in Bharad's slide:
// https://indico.cern.ch/event/1297349/contributions/5454266/attachments/2668534/4624985/2023_06_19_Bharadwaj_EGM_2018_PbPb_ID.pdf


//Endcap region
float phoHoverECut_Endcap_PbPb[4]        = {0.398866,0.394781,0.394773,0.300082};
float phoSigmaIEtaIEtaCut_Endcap_PbPb[4] = {0.044998,0.044996,0.035869,0.031438};
float sumIsoCut_Endcap_PbPb[4]           = {20.911811,12.42409,9.631826,8.543906};

//2017 pp photon cut values from Kaya
// Link same as above
// We use "Set 2: PF iso" here, as we do not find HI iso variables in ppref samples.
//loose (90%), medium (80%), tight (70%), extra tight (60%)
// Barrel region
float phoHoverECut_Barrel_pp[4]          = {0.072266,0.016965,0.014784,0.006778};
float phoSigmaIEtaIEtaCut_Barrel_pp[4]   = {0.010806,0.010015,0.009658,0.009555};
float sumIsoCut_Barrel_pp[4]             = {0.416894,0.160915,0.073941,-0.010780};

//Endcap region
float phoHoverECut_Endcap_pp[4]          = {0.032548,0.014136,0.007197,0.006656};
float phoSigmaIEtaIEtaCut_Endcap_pp[4]   = {0.027323,0.026741,0.026241,0.026134};
float sumIsoCut_Endcap_pp[4]             = {0.970591,0.514475,0.325098,0.063265};

// 2018 PbPb electron cut values from Bharad
// https://docs.google.com/presentation/d/1W1TkKf2rzeAcOtTByQaAFK5OaWHlUVCF0v_T0aF_fQY/edit?usp=sharing

// Veto, Loose, Medium, Tight
// Barrel region
float eleSigmaIEtaIEtaCut_Barrel_PbPb_cent[4] = {0.0147,0.0135,0.0116,0.0104};
float eledEtaSeedAtVtxCut_Barrel_PbPb_cent[4] = {0.0041,0.0038,0.0037,0.0029};
float eledPhiAtVtxCut_Barrel_PbPb_cent[4]     = {0.0853,0.0376,0.0224,0.0206};
float eleHoverEBcCut_Barrel_PbPb_cent[4]      = {0.2733,0.1616,0.1589,0.1459};
float eleEoverPInvCut_Barrel_PbPb_cent[4]     = {0.0367,0.0177,0.0173,0.0105};

float eleSigmaIEtaIEtaCut_Barrel_PbPb_pher[4] = {0.0113,0.0107,0.0101,0.0099};
float eledEtaSeedAtVtxCut_Barrel_PbPb_pher[4] = {0.0037,0.0035,0.0033,0.0026};
float eledPhiAtVtxCut_Barrel_PbPb_pher[4]     = {0.1280,0.0327,0.0210,0.0170};
float eleHoverEBcCut_Barrel_PbPb_pher[4]      = {0.1814,0.1268,0.0311,0.0067};
float eleEoverPInvCut_Barrel_PbPb_pher[4]     = {0.1065,0.0774,0.0701,0.0077};

// Endcap region
float eleSigmaIEtaIEtaCut_Endcap_PbPb_cent[4] = {0.048,0.0466,0.0418,0.0358};
float eledEtaSeedAtVtxCut_Endcap_PbPb_cent[4] = {0.0097,0.0063,0.0062,0.0051};
float eledPhiAtVtxCut_Endcap_PbPb_cent[4]     = {0.2348,0.1186,0.0373,0.0266};
float eleHoverEBcCut_Endcap_PbPb_cent[4]      = {0.1898,0.1317,0.1092,0.0925};
float eleEoverPInvCut_Endcap_PbPb_cent[4]     = {0.0300,0.0201,0.0133,0.0065};

float eleSigmaIEtaIEtaCut_Endcap_PbPb_pher[4] = {0.0376,0.0339,0.0316,0.0288};
float eledEtaSeedAtVtxCut_Endcap_PbPb_pher[4] = {0.0074,0.0067,0.0051,0.0044};
float eledPhiAtVtxCut_Endcap_PbPb_pher[4]     = {0.2085,0.0838,0.0384,0.0266};
float eleHoverEBcCut_Endcap_PbPb_pher[4]      = {0.1138,0.0977,0.0810,0.0655};
float eleEoverPInvCut_Endcap_PbPb_pher[4]     = {0.0237,0.0193,0.0192,0.0123};

// Common region
float eleMissHitsCut[4]                       = {3,1,1,1};
float eleIP3DCut                              = 0.03;

// 2017 pp electron cut values from Bharad
// Link same as above

// Veto, Loose, Medium, Tight
// Barrel region
float eleSigmaIEtaIEtaCut_Barrel_pp[4] = {0.01020,0.00997,0.00948,0.00922};
float eledEtaSeedAtVtxCut_Barrel_pp[4] = {0.00327,0.00272,0.00272,0.00268};
float eledPhiAtVtxCut_Barrel_pp[4]     = {0.06055,0.03419,0.03415,0.02609};
float eleEoverPInvCut_Barrel_pp[4]     = {0.52688,0.03101,0.02762,0.02711};


// Endcap region
float eleSigmaIEtaIEtaCut_Endcap_pp[4] = {0.02957,0.02895,0.02685,0.02679};
float eledEtaSeedAtVtxCut_Endcap_pp[4] = {0.00490,0.00488,0.00379,0.00377};
float eledPhiAtVtxCut_Endcap_pp[4]     = {0.09622,0.02983,0.02812,0.02762};
float eleEoverPInvCut_Endcap_pp[4]     = {0.14600,0.01873,0.01614,0.00472};


bool PassPhotonCuts_PbPb(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB){

  float sumIso = (*ggHi.pho_ecalClusterIsoR3)[idx]+(*ggHi.pho_hcalRechitIsoR3)[idx]+(*ggHi.pho_trackIsoR3PtCut20)[idx];

  if(isEC == false && abs((*ggHi.phoSCEta)[idx]) < 1.442 ){ // Barrel region
    if ((*ggHi.phoHoverE)[idx] > phoHoverECut_Barrel_PbPb[idCut])
      return false;
    else if ((*ggHi.phoSigmaIEtaIEta_2012)[idx] > phoSigmaIEtaIEtaCut_Barrel_PbPb[idCut])
      return false;
    else if ((*ggHi.pho_swissCrx)[idx] >= pho_swissCrxCut )
      return false;
    else if (abs((*ggHi.pho_seedTime)[idx]) >= pho_seedTimeCut )
      return false;
    else if (sumIso > sumIsoCut_Barrel_PbPb[idCut])
      return false;
    else
      return true;
  }else if(isEB == false && abs((*ggHi.phoSCEta)[idx]) > 1.57 && abs((*ggHi.phoSCEta)[idx]) < 2.1){
    if ((*ggHi.phoHoverE)[idx] > phoHoverECut_Endcap_PbPb[idCut])
      return false;
    else if ((*ggHi.phoSigmaIEtaIEta_2012)[idx] > phoSigmaIEtaIEtaCut_Endcap_PbPb[idCut])
      return false;
    else if ((*ggHi.pho_swissCrx)[idx] >= pho_swissCrxCut )
      return false;
    else if (abs((*ggHi.pho_seedTime)[idx]) >= pho_seedTimeCut )
      return false;
    else if (sumIso > sumIsoCut_Endcap_PbPb[idCut])
      return false;
    else
      return true;
  }else{
    return true;
  }
}


bool PassPhotonCuts_pp(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB){

  float sumIso = (*ggHi.pfpIso3subUE)[idx]+(*ggHi.pfcIso3subUE)[idx]+(*ggHi.pfnIso3subUE)[idx];

  if(isEC == false && abs((*ggHi.phoSCEta)[idx]) < 1.44 ){ // Barrel region
    if ((*ggHi.phoHoverE)[idx] > phoHoverECut_Barrel_pp[idCut])
      return false;
    else if ((*ggHi.phoSigmaIEtaIEta_2012)[idx] > phoSigmaIEtaIEtaCut_Barrel_pp[idCut])
      return false;
    else if (sumIso > sumIsoCut_Barrel_pp[idCut])
      return false;
    else
      return true;
  }else if(isEB == false && abs((*ggHi.phoSCEta)[idx]) > 1.57 && abs((*ggHi.phoSCEta)[idx]) < 2.1){
    if ((*ggHi.phoHoverE)[idx] > phoHoverECut_Endcap_pp[idCut])
      return false;
    else if ((*ggHi.phoSigmaIEtaIEta_2012)[idx] > phoSigmaIEtaIEtaCut_Endcap_pp[idCut])
      return false;
    else if (sumIso > sumIsoCut_Endcap_pp[idCut])
      return false;
    else
      return true;
  }else{
    return true;
  }
}

bool PassElectronCuts_PbPb(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB, int hiBin){

  if(hiBin>=0 && hiBin < 60){
    if(isEC == false && abs((*ggHi.eleEta)[idx]) < 1.442 ){ // Barrel region
      if ((*ggHi.eleSigmaIEtaIEta_2012)[idx] > eleSigmaIEtaIEtaCut_Barrel_PbPb_cent[idCut])
        return false;
      else if (abs((*ggHi.eledEtaSeedAtVtx)[idx]) > eledEtaSeedAtVtxCut_Barrel_PbPb_cent[idCut])
        return false;
      else if (abs((*ggHi.eledPhiAtVtx)[idx]) > eledPhiAtVtxCut_Barrel_PbPb_cent[idCut] )
        return false;
      else if ((*ggHi.eleHoverEBc)[idx] > eleHoverEBcCut_Barrel_PbPb_cent[idCut] )
        return false;
      else if (abs((*ggHi.eleEoverPInv)[idx]) > eleEoverPInvCut_Barrel_PbPb_cent[idCut])
        return false;
      else if ((*ggHi.eleMissHits)[idx] > eleMissHitsCut[idCut])
        return false;
      else if ((*ggHi.eleIP3D)[idx] >= eleIP3DCut)
        return false;
      else
        return true;
    }else if(isEB == false && abs((*ggHi.eleEta)[idx]) > 1.57 && abs((*ggHi.eleEta)[idx]) < 2.4){
      if ((*ggHi.eleSigmaIEtaIEta_2012)[idx] > eleSigmaIEtaIEtaCut_Endcap_PbPb_cent[idCut])
        return false;
      else if (abs((*ggHi.eledEtaSeedAtVtx)[idx]) > eledEtaSeedAtVtxCut_Endcap_PbPb_cent[idCut])
        return false;
      else if (abs((*ggHi.eledPhiAtVtx)[idx]) > eledPhiAtVtxCut_Endcap_PbPb_cent[idCut] )
        return false;
      else if ((*ggHi.eleHoverEBc)[idx] > eleHoverEBcCut_Endcap_PbPb_cent[idCut] )
        return false;
      else if (abs((*ggHi.eleEoverPInv)[idx]) > eleEoverPInvCut_Endcap_PbPb_cent[idCut])
        return false;
      else if ((*ggHi.eleMissHits)[idx] > eleMissHitsCut[idCut])
        return false;
      else if ((*ggHi.eleIP3D)[idx] >= eleIP3DCut)
        return false;
      else
        return true;
    }else{
      return true;
    }
  }else if(hiBin>=60 && hiBin < 200){
    if(isEC == false && abs((*ggHi.eleEta)[idx]) < 1.442 ){ // Barrel region
      if ((*ggHi.eleSigmaIEtaIEta_2012)[idx] > eleSigmaIEtaIEtaCut_Barrel_PbPb_pher[idCut])
        return false;
      else if (abs((*ggHi.eledEtaSeedAtVtx)[idx]) > eledEtaSeedAtVtxCut_Barrel_PbPb_pher[idCut])
        return false;
      else if (abs((*ggHi.eledPhiAtVtx)[idx]) > eledPhiAtVtxCut_Barrel_PbPb_pher[idCut] )
        return false;
      else if ((*ggHi.eleHoverEBc)[idx] > eleHoverEBcCut_Barrel_PbPb_pher[idCut] )
        return false;
      else if (abs((*ggHi.eleEoverPInv)[idx]) > eleEoverPInvCut_Barrel_PbPb_pher[idCut])
        return false;
      else if ((*ggHi.eleMissHits)[idx] > eleMissHitsCut[idCut])
        return false;
      else if ((*ggHi.eleIP3D)[idx] >= eleIP3DCut)
        return false;
      else
        return true;
    }else if(isEB == false && abs((*ggHi.eleEta)[idx]) > 1.57 && abs((*ggHi.eleEta)[idx]) < 2.4){
      if ((*ggHi.eleSigmaIEtaIEta_2012)[idx] > eleSigmaIEtaIEtaCut_Endcap_PbPb_pher[idCut])
        return false;
      else if (abs((*ggHi.eledEtaSeedAtVtx)[idx]) > eledEtaSeedAtVtxCut_Endcap_PbPb_pher[idCut])
        return false;
      else if (abs((*ggHi.eledPhiAtVtx)[idx]) > eledPhiAtVtxCut_Endcap_PbPb_pher[idCut] )
        return false;
      else if ((*ggHi.eleHoverEBc)[idx] > eleHoverEBcCut_Endcap_PbPb_pher[idCut] )
        return false;
      else if (abs((*ggHi.eleEoverPInv)[idx]) > eleEoverPInvCut_Endcap_PbPb_pher[idCut])
        return false;
      else if ((*ggHi.eleMissHits)[idx] > eleMissHitsCut[idCut])
        return false;
      else if ((*ggHi.eleIP3D)[idx] >= eleIP3DCut)
        return false;
      else
        return true;
    }else{
      return true;
    }
  }
  else{
    return true;
  }
}


bool PassElectronCuts_pp(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB){

  if(isEC == false && abs((*ggHi.eleEta)[idx]) < 1.442 ){ // Barrel region
    if ((*ggHi.eleSigmaIEtaIEta_2012)[idx] > eleSigmaIEtaIEtaCut_Barrel_pp[idCut])
      return false;
    else if (abs((*ggHi.eledEtaSeedAtVtx)[idx]) > eledEtaSeedAtVtxCut_Barrel_pp[idCut])
      return false;
    else if (abs((*ggHi.eledPhiAtVtx)[idx]) > eledPhiAtVtxCut_Barrel_pp[idCut] )
      return false;
    else if (abs((*ggHi.eleEoverPInv)[idx]) > eleEoverPInvCut_Barrel_pp[idCut])
      return false;
    else
      return true;
  }else if(isEB == false && abs((*ggHi.eleEta)[idx]) > 1.57 && abs((*ggHi.eleEta)[idx]) < 2.4){
    if ((*ggHi.eleSigmaIEtaIEta_2012)[idx] > eleSigmaIEtaIEtaCut_Endcap_pp[idCut])
      return false;
    else if (abs((*ggHi.eledEtaSeedAtVtx)[idx]) > eledEtaSeedAtVtxCut_Endcap_pp[idCut])
      return false;
    else if (abs((*ggHi.eledPhiAtVtx)[idx]) > eledPhiAtVtxCut_Endcap_pp[idCut] )
      return false;
    else if (abs((*ggHi.eleEoverPInv)[idx]) > eleEoverPInvCut_Endcap_pp[idCut])
      return false;
    else
      return true;
  }else{
    return true;
  }

}

bool PassCuts(ggHiNtuplizer& ggHi, int idx, int idCut, bool isEC, bool isEB, bool isPbPb, bool isEle, int hiBin){
  if(isPbPb){
    if(isEle)
      return PassElectronCuts_PbPb(ggHi,idx,idCut,isEC,isEB,hiBin);
    else
      return PassPhotonCuts_PbPb(ggHi,idx,idCut,isEC,isEB);
  }else{
    if(isEle)
      return PassElectronCuts_pp(ggHi,idx,idCut,isEC,isEB);
    else
      return PassPhotonCuts_pp(ggHi,idx,idCut,isEC,isEB);
  }
}

