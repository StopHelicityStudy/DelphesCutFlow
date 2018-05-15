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

//------------------------------------------------------------------------------

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

  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0);
  TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
  TH1 *lepJetMass = new TH1F("lepJetMass", "lepJetMass", 100, 0.0, 300);
  TH1 *lepJetPT = new TH1F("lepJetPT", "lepJetPT", 100, 0.0, 1000);
  TH1 *lepPt = new TH1F("lepPt", "lepPt", 100, 0, 1000);
  TH1 *lepEta = new TH1F("lepEta", "lepEta", 100, -4, 4);
  TH1 *lepRatio = new TH1F("lepRatio", "jet lepRatio{T}", 100, 0.0, 1.0);
  TH1 *met = new TH1F("met", "met", 100, 0.0, 1000);
  TCanvas c1;

  int muonCount = 0;
  int totalCount = 0;
  bool verbose = true;


  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(verbose) std::cout << "new event\n" << std::endl;
    if(verbose) std::cout << "\t nJets: " << branchJet->GetEntries() << " n leptons (electron or muon) " << branchMuon->GetEntries() + branchElectron->GetEntries() << std::endl;

    //number of events

    totalCount = totalCount +1;
    //require at least 2 jets and one muon
    if ( (branchJet->GetEntries() < 2 || branchMuon->GetEntries() < 1) ) continue;

    muonCount = muonCount +1;
    muon = (Muon*) branchMuon->At(0);


    MET = (MissingET*) branchMet->At(0);
    if(verbose) std::cout << "\t muon pt " <<  muon->P4().Pt() << std::endl;


  
    // If event contains at least 1 jet

    
    if(branchJet->GetEntries() > 0)
    {

      for(int i=0; i < branchJet->GetEntries(); i++){
        lepJet = (Jet*) branchJet->At(i);
        if(verbose) std::cout << "\t muon jet DeltaR, btag " << muon->P4().DeltaR(lepJet->P4())  << " " << lepJet->BTag<< std::endl;
        if(lepJet->PT > 30 && muon->P4().DeltaR(lepJet->P4()) < 1.4 && lepJet->BTag == 1) break;
      }
    }


    if (muon->PT < 30 || lepJet->PT < 30 ) continue;


    lepEta->Fill(muon->Eta);
    lepPt->Fill(muon->PT);
    lepRatio->Fill(lepJet->PT/(lepJet->PT + muon->P4().Pt()) );

    lepJetPT->Fill(lepJet->PT);

    lepJetMass->Fill( (lepJet->P4() + muon->P4()).M() );


    met->Fill(MET->MET);


    /*GenParticle *genPart;

    for(int i=0;i<branchGen->GetEntries();i++){
      genPart = (GenParticle*) branchGen->At(i);
      std::cout << "gen " << genPart->PT << std::endl;
    }*/

     if(verbose) std::cout << "end event\n" << std::endl;
  }


  std::cout << "muon events vs total " << muonCount << " " << totalCount << std::endl;

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

