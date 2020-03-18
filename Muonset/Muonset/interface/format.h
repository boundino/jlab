// vim:set ts=4 sw=4 fdm=marker et:
#ifndef _XBFRAMEFORMAT_H_
#define _XBFRAMEFORMAT_H_

#define MAX_MUON     10000

#include <memory>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <list>
#include <iomanip>
#include <cmath>

//#include "Bfinder/Bfinder/interface/TriggerBooking.h"
#include <TLorentzVector.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TString.h>
#include <TNtuple.h>
#include <TVector3.h>

class MuonInfoBranches
{
public:
  int RunNo;
  int EvtNo;
  int LumiNo;
  int size;
  int charge[MAX_MUON];
  float pt[MAX_MUON];
  float eta[MAX_MUON];
  float phi[MAX_MUON];
  // bool isTrackerMuon[MAX_MUON];
  // bool isGlobalMuon[MAX_MUON];
  int i_nStripLayer[MAX_MUON];
  int i_nPixelLayer[MAX_MUON];
  float dzPV[MAX_MUON];
  float dxyPV[MAX_MUON];

  void regTree(TTree *root)
  {
    root->Branch("EvtInfo.RunNo", &RunNo, "EvtInfo.RunNo/I");
    root->Branch("EvtInfo.EvtNo", &EvtNo, "EvtInfo.EvtNo/I");
    root->Branch("EvtInfo.LumiNo", &LumiNo, "EvtInfo.LumiNo/I");
    root->Branch("MuonInfo.size", &size, "MuonInfo.size/I");
    root->Branch("MuonInfo.charge", charge, "MuonInfo.charge[MuonInfo.size]/I");
    root->Branch("MuonInfo.pt", pt, "MuonInfo.pt[MuonInfo.size]/F");
    root->Branch("MuonInfo.eta", eta, "MuonInfo.eta[MuonInfo.size]/F");
    root->Branch("MuonInfo.phi", phi, "MuonInfo.phi[MuonInfo.size]/F");
    // root->Branch("MuonInfo.isTrackerMuon", isTrackerMuon, "MuonInfo.isTrackerMuon[MuonInfo.size]/O");
    // root->Branch("MuonInfo.isGlobalMuon", isGlobalMuon, "MuonInfo.isGlobalMuon[MuonInfo.size]/O");
    root->Branch("MuonInfo.i_nStripLayer", i_nStripLayer, "MuonInfo.i_nStripLayer[MuonInfo.size]/I");
    root->Branch("MuonInfo.i_nPixelLayer", i_nPixelLayer, "MuonInfo.i_nPixelLayer[MuonInfo.size]/I");
    root->Branch("MuonInfo.dzPV", dzPV, "MuonInfo.dzPV[MuonInfo.size]/F");
    root->Branch("MuonInfo.dxyPV", dxyPV, "MuonInfo.dxyPV[MuonInfo.size]/F");

  }

  void setbranchadd(TTree *root)
  {
    root->SetBranchAddress("EvtInfo.RunNo", &RunNo);
    root->SetBranchAddress("EvtInfo.EvtNo", &EvtNo);
    root->SetBranchAddress("EvtInfo.LumiNo", &LumiNo);
    root->SetBranchAddress("MuonInfo.size", &size);
    root->SetBranchAddress("MuonInfo.charge", charge);
    root->SetBranchAddress("MuonInfo.pt", pt);
    root->SetBranchAddress("MuonInfo.eta", eta);
    root->SetBranchAddress("MuonInfo.phi", phi);
    // root->SetBranchAddress("MuonInfo.isTrackerMuon", isTrackerMuon);
    // root->SetBranchAddress("MuonInfo.isGlobalMuon", isGlobalMuon);
    root->SetBranchAddress("MuonInfo.i_nStripLayer", i_nStripLayer);
    root->SetBranchAddress("MuonInfo.i_nPixelLayer", i_nPixelLayer);
    root->SetBranchAddress("MuonInfo.dzPV", dzPV);
    root->SetBranchAddress("MuonInfo.dxyPV", dxyPV);
  }

};

#endif
