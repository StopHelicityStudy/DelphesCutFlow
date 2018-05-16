/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l Cutflow.C'("stop_a/stop_a/stop_a_1.root")' -b
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <vector>

//------------------------------------------------------------------------------


TLorentzVector dropJet(TLorentzVector TL_hadTop, TClonesArray *branchJet, TLorentzVector TL_lepton, bool verbose)
{
  Jet *hadSubJet;
  int nJets;
  double minMass= 10000; 
  double TopMass = 172.5;
  TLorentzVector droppedJet;
  std::cout << "drop jet \n";

  for(int i=0; i < branchJet->GetEntries(); i++){
    hadSubJet = (Jet*) branchJet->At(i);
    if (hadSubJet->PT < 25) continue;
    if (TL_lepton.DeltaR(hadSubJet->P4()) < 1.5) continue;

    double deltaR= TL_hadTop.DeltaR(hadSubJet->P4());
    if(verbose) std::cout << "\t mass of top " << (TL_hadTop-hadSubJet->P4()).M() << " " << deltaR << std::endl;

    minMass = min(minMass, (TL_hadTop-hadSubJet->P4()).M());

    if (minMass == (TL_hadTop-hadSubJet->P4()).M() ){
      droppedJet = hadSubJet->P4();
      minMass = (TL_hadTop-droppedJet).M();
      std::cout << "got here \n";
    }
    
    nJets = nJets + 1;

  }
  std::cout << nJets << " " << minMass << std::endl;

  TLorentzVector newTopCandidate = TL_hadTop - droppedJet;

  //if no jets left, return original top

  //if no jets left, return original top
  if(nJets ==0)
  {
    return TL_hadTop;
  }
  //if original top is better than new top, and we are less than 4 jets, return original top
  if (nJets < 3 && abs(newTopCandidate.M() - TopMass) > abs(TL_hadTop.M() - TopMass) )
  {
    return TL_hadTop;
  }

  //if original top is worse or we have to many jets, then run again
  return dropJet(newTopCandidate, branchJet, TL_lepton, verbose);

}

void Cutflow(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMet = treeReader->UseBranch("MissingET");
  TClonesArray *branchGen = treeReader->UseBranch("GenParticle");

  //
  Muon *muon;
  MissingET *MET;
  Jet *lepJet;
  Jet *hadSubJet;
  Electron *electron;

  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0);
  TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
  TH1 *lepJetMass = new TH1F("lepJetMass", "lepJetMass", 100, 0.0, 300);
  TH1 *lepJetPT = new TH1F("lepJetPT", "lepJetPT", 100, 0.0, 1000);
  TH1 *lepPt = new TH1F("lepPt", "lepPt", 100, 0, 1000);
  TH1 *lepEta = new TH1F("lepEta", "lepEta", 100, -4, 4);
  TH1 *lepRatio = new TH1F("lepRatio", "jet lepRatio{T}", 100, 0.0, 1.0);
  TH1 *met = new TH1F("met", "met", 100, 0.0, 1000);
  TH1 *TH_cutflow = new TH1F("TH_cutflow", "cutflow", 10, -.5, 9.5);
  TCanvas c1;

  int muonCount = 0;
  int totalCount = 0;
  bool verbose = true;
  int cutflow_count = 0;
  int nBJets_near_lep = 0;
  int nBJets_int_top = 0;

  std::vector<TLorentzVector> TL_v;

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {


    TLorentzVector TL_lepton;
    TLorentzVector TL_leptJet;
    TLorentzVector TL_hadTop;


    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(verbose) std::cout << "new event\n" << std::endl;
    if(verbose) std::cout << "\t nJets: " << branchJet->GetEntries() << " n leptons (electron or muon) " << branchMuon->GetEntries() + branchElectron->GetEntries() << std::endl;

    //number of events
    totalCount = totalCount +1;

    //get lepton
    TL_lepton.SetPtEtaPhiM(0,0,0,0);
    if(branchMuon->GetEntries() > 0)
    {
      muon = (Muon*) branchMuon->At(0);
      TL_lepton.SetPtEtaPhiM(muon->P4().Pt(),muon->P4().Eta(),muon->P4().Phi(),muon->P4().M()) ;
    } 
    if(branchElectron->GetEntries() > 0)
    {
      electron = (Electron*) branchElectron->At(0);
      TL_lepton.SetPtEtaPhiM(electron->P4().Pt(),electron->P4().Eta(),electron->P4().Phi(),electron->P4().M()) ;
    }
    if(verbose) std::cout << "\t lept pt " <<  TL_lepton.Pt() << std::endl;


    //get met
    MET = (MissingET*) branchMet->At(0);


    // select for highest bjet near muon
    if(branchJet->GetEntries() > 0 and TL_lepton.Pt() > 0)
    {
      nBJets_near_lep = 0;
      for(int i=0; i < branchJet->GetEntries(); i++){
        lepJet = (Jet*) branchJet->At(i);
        if (lepJet->PT < 25) continue;
        if (lepJet->BTag == 0) continue;
        if (TL_lepton.DeltaR(lepJet->P4()) > 1.5) continue;

        TL_leptJet = lepJet->P4();

        if(verbose) std::cout << "\t lept jet DeltaR, btag " << TL_lepton.DeltaR(TL_leptJet)  << " " << lepJet->BTag<< std::endl;
        nBJets_near_lep = nBJets_near_lep + 1;
      }
    }
    if(verbose) std::cout << "\t selected: lept jet DeltaR, btag, nBJets_near_lep " << TL_lepton.DeltaR(TL_leptJet) << " " << nBJets_near_lep<< " " <<  std::endl;

    TL_v.clear();
    // select for top subjets
    if(branchJet->GetEntries() > 0 and TL_lepton.Pt() > 0)
    {
      nBJets_int_top = 0;
      for(int i=0; i < branchJet->GetEntries(); i++){
        hadSubJet = (Jet*) branchJet->At(i);
        if (hadSubJet->PT < 25) continue;
        if (TL_lepton.DeltaR(hadSubJet->P4()) < 1.5) continue;

        TL_hadTop = TL_hadTop + hadSubJet->P4();

        TL_v.push_back(hadSubJet->P4());
        if(verbose) std::cout << "\t mass of top " << TL_hadTop.M() << std::endl;
        nBJets_int_top = nBJets_int_top + 1;
      }
      std::cout << "\n";

      //TL_hadTop = dropJet(TL_hadTop, branchJet, TL_lepton, verbose);
    }



/*

    cutflow_count = 0;
    TH_cutflow->Fill(cutflow_count);
    //TH_cutflow->GetXaxis()->SetBinLabel(cutflow_count,"all muon");
    cutflow_count = cutflow_count + 1;


    // cut on muon pt and jet pt
    if (muon->PT < 30 || lepJet->PT < 30 ) continue;


    //fill hists
    lepEta->Fill(muon->Eta);
    lepPt->Fill(muon->PT);
    lepRatio->Fill(lepJet->PT/(lepJet->PT + TL_lepton.Pt()) );
    lepJetPT->Fill(lepJet->PT);
    lepJetMass->Fill( (lepJet->P4() + TL_lepton).M() );
    met->Fill(MET->MET);

    if(verbose) std::cout << "end event\n" << std::endl;
    */
  }


  //end job
  std::cout << "muon events vs total " << muonCount << " " << totalCount << std::endl;


  TH_cutflow->Draw();
  c1.SaveAs("output/cutflow.png");

  // Show resulting histograms
  histJetPT->Draw();
  c1.SaveAs("output/test.png");

  lepPt->Draw();
  c1.SaveAs("output/lepPt.png");

  lepJetMass->Draw();
  c1.SaveAs("output/lepJetMass.png");

  lepEta->Draw();
  c1.SaveAs("output/lepEta.png");

  lepRatio->Draw();
  c1.SaveAs("output/lepRatio.png");

  lepJetPT->Draw();
  c1.SaveAs("output/lepJetPT.png");

  met->Draw();
  c1.SaveAs("output/met.png");

}

