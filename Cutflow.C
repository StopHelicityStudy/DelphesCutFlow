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
#include <math.h>

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
  Jet *hadSubJet2;
  Jet *hadSubJet3;
  Electron *electron;

  TFile f("Stop.root","recreate");


  // Book histograms

  TH1 *TH_had_top_mass = new TH1F("had_top_mass", "had_top_mass", 100, 100, 250);
  TH1 *TH_had_top_pt = new TH1F("had_top_pt", "had_top_pt", 100, 0.0, 1000);
  TH1 *TH_had_top_eta = new TH1F("had_top_eta", "had_top_eta", 100, -4, 4);

  TH1 *TH_lept_top_mass = new TH1F("lept_top_mass", "lept_top_mass", 100, 100, 250);
  TH1 *TH_lept_top_pt = new TH1F("lept_top_pt", "lept_top_pt", 100, 0.0, 1000);
  TH1 *TH_lept_top_eta = new TH1F("lept_top_eta", "lept_top_eta", 100, -4, 4);

  TH1 *TH_lep_plus_JetMass = new TH1F("lep_JetMass", "lep_JetMass", 100, 0.0, 300);
  TH1 *TH_lepJetPT = new TH1F("lepJetPT", "lepJetPT", 100, 0.0, 1000);
  TH1 *TH_lepJetEta = new TH1F("lepJetEta", "lepJetEta", 100, -4, 4);

  TH1 *TH_lepPt = new TH1F("lepPt", "lepPt", 100, 0, 1000);
  TH1 *TH_lepEta = new TH1F("lepEta", "lepEta", 100, -4, 4);

  TH1 *TH_lepRatio = new TH1F("lepRatio", "jet lepRatio", 10, 0.0, 1.0);
  TH1 *TH_recLepRatio = new TH1F("recLepRatio", "jet recLepRatio", 10, 0.0, 1.0);
  TH1 *TH_met = new TH1F("met", "met", 100, 0.0, 1000);

  TH1 *TH_cos_lep_MET = new TH1F("cos_lep_MET", "cos_lep_MET", 100, -1, 1);
  TH1 *TH_W_guess_mass = new TH1F("W_guess_mass", "W_guess_mass", 100, 50, 400);
  TH1 *TH_DeltaR_lepJet_lepton = new TH1F("DeltaR_lepJet_lepton", "DeltaR_lepJet_lepton", 100, 0.0, 1.5);



  TH1 *TH_noCuts_had_top_mass = new TH1F("noCuts_had_top_mass", "had_top_mass", 100, 100, 250);
  TH1 *TH_noCuts_had_top_pt = new TH1F("noCuts_had_top_pt", "had_top_pt", 100, 0.0, 1000);
  TH1 *TH_noCuts_had_top_eta = new TH1F("noCuts_had_top_eta", "had_top_eta", 100, -4, 4);

  TH1 *TH_noCuts_lept_top_mass = new TH1F("noCuts_lept_top_mass", "lept_top_mass", 100, 100, 250);
  TH1 *TH_noCuts_lept_top_pt = new TH1F("noCuts_lept_top_pt", "lept_top_pt", 100, 0.0, 1000);
  TH1 *TH_noCuts_lept_top_eta = new TH1F("noCuts_lept_top_eta", "lept_top_eta", 100, -4, 4);

  TH1 *TH_noCuts_lep_plus_JetMass = new TH1F("noCuts_lep_JetMass", "lep_JetMass", 100, 0.0, 300);
  TH1 *TH_noCuts_lepJetPT = new TH1F("noCuts_lepJetPT", "lepJetPT", 100, 0.0, 1000);
  TH1 *TH_noCuts_lepJetEta = new TH1F("noCuts_lepJetEta", "lepJetEta", 100, -4, 4);

  TH1 *TH_noCuts_lepPt = new TH1F("noCuts_lepPt", "lepPt", 100, 0, 1000);
  TH1 *TH_noCuts_lepEta = new TH1F("noCuts_lepEta", "lepEta", 100, -4, 4);

  TH1 *TH_noCuts_lepRatio = new TH1F("noCuts_lepRatio", "jet lepRatio", 10, 0.0, 1.0);
  TH1 *TH_noCuts_recLepRatio = new TH1F("noCuts_recLepRatio", "jet recLepRatio", 10, 0.0, 1.0);
  TH1 *TH_noCuts_met = new TH1F("noCuts_met", "met", 100, 0.0, 1000);

  TH1 *TH_noCuts_cos_lep_MET = new TH1F("noCuts_cos_lep_MET", "cos_lep_MET", 100, -1, 1);
  TH1 *TH_noCuts_W_guess_mass = new TH1F("noCuts_W_guess_mass", "W_guess_mass", 100, 50, 400);
  TH1 *TH_noCuts_DeltaR_lepJet_lepton = new TH1F("noCuts_DeltaR_lepJet_lepton", "DeltaR_lepJet_lepton", 100, 0.0, 1.5);

  TH1 *TH_cutflow = new TH1F("TH_cutflow", "cutflow", 10, -.5, 9.5);




  TCanvas c1;

  int muonCount = 0;
  int totalCount = 0;
  bool verbose = false;
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
    TLorentzVector TL_lepTop_guess;
    TLorentzVector TL_neutrino_guess;
    TLorentzVector TL_W_guess;
    TLorentzVector TL_object_closest_to_MET;


    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(verbose) std::cout << "new event\n" << std::endl;
    if(verbose) std::cout << "\t nJets: " << branchJet->GetEntries() << " n leptons (electron or muon) " << branchMuon->GetEntries() + branchElectron->GetEntries() << std::endl;

    //number of events
    totalCount = totalCount +1;


    //get met
    MET = (MissingET*) branchMet->At(0);
    TL_object_closest_to_MET = -MET->P4();


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


    // select for highest bjet near lepton
     std::cout << "\tone jet combos \n";
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


    // find best three jet combo
    float bestMass = 1000;
    float nextBest = 1000;
    float topMass = 172.5;
    if(branchJet->GetEntries() > 0 and TL_lepton.Pt() > 0)
    {
      nBJets_int_top = 0;
      //one jet
      for(int i=0; i < branchJet->GetEntries(); i++){
        hadSubJet = (Jet*) branchJet->At(i);
        if (hadSubJet->PT < 25) continue;
        if (TL_lepton.DeltaR(hadSubJet->P4()) < 1.5) continue;

        if(verbose) std::cout << "\t \t mass of top " << hadSubJet->P4().M() << std::endl;

        if( abs( hadSubJet->P4().M() - topMass)  <  abs( TL_hadTop.M() - topMass)  )
        {
          nextBest = bestMass;
          TL_hadTop = hadSubJet->P4();
          bestMass = hadSubJet->P4().M();
        }

        nBJets_int_top = nBJets_int_top + 1;
      }
      std::cout << "\ttwo jet combos \n";


      // two jet combinations
      TLorentzVector tempTop;
      for(int i=0; i < branchJet->GetEntries(); i++){
        hadSubJet = (Jet*) branchJet->At(i);
        if (hadSubJet->PT < 25) continue;
        if (TL_lepton.DeltaR(hadSubJet->P4()) < 1.5) continue;


        for(int j=i+1; j < branchJet->GetEntries(); j++){
          hadSubJet2 = (Jet*) branchJet->At(j);
          if (hadSubJet2->PT < 25) continue;
          if (TL_lepton.DeltaR(hadSubJet2->P4()) < 1.5) continue;

          tempTop = (hadSubJet->P4()+hadSubJet2->P4());
          if(verbose) std::cout << "\t \t mass of top " << tempTop.M() << " " << i << " " << j << std::endl;
  
          if( abs( tempTop.M() - topMass)  <  abs( TL_hadTop.M() - topMass)  )
          {
            nextBest = bestMass;
            TL_hadTop = tempTop;
            bestMass = tempTop.M();
          }
        }
      }

      std::cout << "\tthree jet combos \n";
      // three jet combinations
      for(int i=0; i < branchJet->GetEntries(); i++){
        hadSubJet = (Jet*) branchJet->At(i);
        if (hadSubJet->PT < 25) continue;
        if (TL_lepton.DeltaR(hadSubJet->P4()) < 1.5) continue;


        for(int j=i+1; j < branchJet->GetEntries(); j++){
          hadSubJet2 = (Jet*) branchJet->At(j);
          if (hadSubJet2->PT < 25) continue;
          if (TL_lepton.DeltaR(hadSubJet2->P4()) < 1.5) continue;
            for(int j=i+1; j < branchJet->GetEntries(); j++){
              hadSubJet3 = (Jet*) branchJet->At(j);
              if (hadSubJet3->PT < 25) continue;
              if (TL_lepton.DeltaR(hadSubJet3->P4()) < 1.5) continue;

              tempTop = (hadSubJet->P4()+hadSubJet2->P4()+hadSubJet3->P4());
              if(verbose) std::cout << "\t \t mass of top " << tempTop.M() << " " << i << " " << j << std::endl;
      
              if( abs( tempTop.M() - topMass)  <  abs( TL_hadTop.M() - topMass)  )
              {
                nextBest = bestMass;
                TL_hadTop = tempTop;
                bestMass = tempTop.M();
              }
          }
        }
      }
      std::cout << "\n";


    }

     if(verbose) std::cout << "\t mass of top " << TL_hadTop.M() << " next best " << nextBest << std::endl;




    TH_noCuts_lepEta->Fill(TL_leptJet.Eta());
    TH_noCuts_lepPt->Fill(TL_leptJet.Pt());
    TH_noCuts_lepRatio->Fill(TL_leptJet.Pt()/(TL_leptJet.Pt() + TL_lepton.Pt()) );
    TH_noCuts_lepJetPT->Fill(TL_leptJet.Pt());
    TH_noCuts_lepJetEta->Fill(TL_leptJet.Eta());
    TH_noCuts_lep_plus_JetMass->Fill( (TL_leptJet + TL_lepton).M() );
    TH_noCuts_met->Fill(MET->MET);
    TH_noCuts_had_top_mass->Fill(TL_hadTop.M());
    TH_noCuts_had_top_pt->Fill(TL_hadTop.Pt());
    TH_noCuts_had_top_eta->Fill(TL_hadTop.Eta());
    TH_noCuts_lept_top_mass->Fill(TL_lepTop_guess.M());
    TH_noCuts_lept_top_pt->Fill(TL_lepTop_guess.Pt());
    TH_noCuts_lept_top_eta->Fill(TL_lepTop_guess.Eta());
    TH_noCuts_recLepRatio->Fill(TL_leptJet.E()/TL_lepTop_guess.E());
    TH_noCuts_cos_lep_MET->Fill(cos(TL_lepton.DeltaPhi(MET->P4()) ));
    TH_noCuts_W_guess_mass->Fill(TL_W_guess.M());
    TH_noCuts_DeltaR_lepJet_lepton->Fill(TL_leptJet.DeltaR(TL_lepton));



    cutflow_count = 0;
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;

    if (MET->MET < 125) continue; 
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;


    if (TL_lepton.Pt() < 20 || abs(TL_lepton.Eta()) > 2.4) continue; 
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;

    if (TL_leptJet.Pt() < 25 || abs(TL_leptJet.Eta()) > 2.4) continue; 
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;

    if( abs(cos(TL_lepton.DeltaPhi(MET->P4()) ) ) > 0.7) continue;
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;

    if(  abs( TL_hadTop.M() - topMass)  > 15) continue;
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;

    TL_lepTop_guess = -TL_hadTop;
    TL_neutrino_guess = TL_lepTop_guess -TL_leptJet - TL_lepton;
    TL_W_guess = TL_lepton + TL_neutrino_guess;

    float massW = 80.4;
    if(  abs( TL_W_guess.M() - massW)  < 40) continue;
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;


    if(  TL_leptJet.DeltaR(TL_lepton)  > 1.5 ) continue; //pointless cut as it is required already
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;

    if(abs(MET->P4().DeltaPhi(TL_lepton))  < abs(MET->P4().DeltaPhi(TL_object_closest_to_MET)) )
    {
      TL_object_closest_to_MET = TL_lepton;
    }
    if(abs(MET->P4().DeltaPhi(TL_hadTop))  < abs(MET->P4().DeltaPhi(TL_object_closest_to_MET)) )
    {
      TL_object_closest_to_MET = TL_hadTop;
    }
    if(abs(MET->P4().DeltaPhi(TL_leptJet))  < abs(MET->P4().DeltaPhi(TL_object_closest_to_MET)) )
    {
      TL_object_closest_to_MET = TL_leptJet;
    }
    if(verbose) std::cout << " \tDeltaR met and closest object " << abs(MET->P4().DeltaPhi(TL_object_closest_to_MET)) << std::endl;


    /*if(  abs(MET->P4().DeltaR(TL_object_closest_to_MET))  > 1.5  ||  abs(MET->P4().DeltaR(TL_object_closest_to_MET))  < 0.8 ) continue; //pointless cut as it is required already
    TH_cutflow->Fill(cutflow_count);
    cutflow_count = cutflow_count + 1;*/


    std::cout << "fill hists" <<  std::endl;
    TH_lepEta->Fill(TL_leptJet.Eta());
    TH_lepPt->Fill(TL_leptJet.Pt());
    TH_lepRatio->Fill(TL_leptJet.Pt()/(TL_leptJet.Pt() + TL_lepton.Pt()) );
    TH_lepJetPT->Fill(TL_leptJet.Pt());
    TH_lepJetEta->Fill(TL_leptJet.Eta());
    TH_lep_plus_JetMass->Fill( (TL_leptJet + TL_lepton).M() );
    TH_met->Fill(MET->MET);
    TH_had_top_mass->Fill(TL_hadTop.M());
    TH_had_top_pt->Fill(TL_hadTop.Pt());
    TH_had_top_eta->Fill(TL_hadTop.Eta());
    TH_lept_top_mass->Fill(TL_lepTop_guess.M());
    TH_lept_top_pt->Fill(TL_lepTop_guess.Pt());
    TH_lept_top_eta->Fill(TL_lepTop_guess.Eta());
    TH_recLepRatio->Fill(TL_leptJet.E()/TL_lepTop_guess.E());
    TH_cos_lep_MET->Fill(cos(TL_lepton.DeltaPhi(MET->P4()) ));
    TH_W_guess_mass->Fill(TL_W_guess.M());
    TH_DeltaR_lepJet_lepton->Fill(TL_leptJet.DeltaR(TL_lepton));

    
  }


  //end job
  //std::cout << "muon events vs total " << muonCount << " " << totalCount << std::endl;


  TH_cutflow->Draw();
  c1.SaveAs("output/cutflow.png");

  // Show resulting histograms

  TH_lepPt->Draw();
  c1.SaveAs("output/lepPt.png");

  TH_lep_plus_JetMass->Draw();
  c1.SaveAs("output/lep_plus_JetMass.png");

  TH_lepEta->Draw();
  c1.SaveAs("output/lepEta.png");

  TH_lepRatio->Draw();
  c1.SaveAs("output/lepRatio.png");

  TH_lepJetPT->Draw();
  c1.SaveAs("output/lepJetPT.png");

  TH_met->Draw();
  c1.SaveAs("output/met.png");

  f.Write();

}

