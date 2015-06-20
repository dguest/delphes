/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l examples/ExampleTrack.C'("delphes_output.root")'
*/

//------------------------------------------------------------------------------

struct TestPlots
{
  TH1 *fOldTrackD0, *fOldTrackZ0, *fOldTrackPhi, *fOldTrackEta, *fOldTrackPt;
  TH1 *fTrackD0, *fTrackZ0, *fTrackPhi, *fTrackEta, *fTrackPt, *fTrackD0Sig, *fTrackPt2;
  TH1 *fTrackDeltaD0, *fTrackDeltaZ0, *fTrackDeltaPhi, *fTrackDeltaEta, *fTrackDeltaPt;
  
  TH2 *fTrackDeltaD0vsEta, *fTrackDeltaZ0vsEta, *fTrackDeltaPhivsEta, *fTrackDeltaThetavsEta, *fTrackDeltaPtvsEta;
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  plots->fOldTrackD0 = result->AddHist1D("oldd0", "original d0", "d0","number of tracks", 100,-0.1,0.1);
  plots->fOldTrackZ0 = result->AddHist1D("oldz0", "orignal z0", "z0","number of tracks", 100,-10,10);
  plots->fOldTrackPhi = result->AddHist1D("oldphi", "original phi", "phi","number of tracks", 64,-3.2,3.2);
  plots->fOldTrackEta= result->AddHist1D("oldeta", "original eta", "eta","number of tracks", 100,-3,3);
//  plots->fOldTrackPt = result->AddHist1D("oldpt", "original pt", "pt","number of tracks", 100,0,200);
  plots->fOldTrackPt = result->AddHist1D("oldpt", "pt", "pt","number of tracks", 100,0,200.0);

  plots->fTrackD0Sig = result->AddHist1D("d0Sig", "d0/d0Err", "d0/d0Err","number of tracks", 100,-10,10);
  plots->fTrackD0 = result->AddHist1D("d0", "d0", "d0","number of tracks", 100,-0.1,0.1);
  plots->fTrackZ0 = result->AddHist1D("z0", "z0", "z0","number of tracks", 100,-10,10);
  plots->fTrackPhi = result->AddHist1D("phi", "phi", "phi","number of tracks", 64,-3.2,3.2);
  plots->fTrackEta= result->AddHist1D("eta", "eta", "eta","number of tracks", 100,-3,3);
  //plots->fTrackPt = result->AddHist1D("pt", "pt", "pt","number of tracks", 100,0,200);
  plots->fTrackPt = result->AddHist1D("PT" , "pt", "pt","number of tracks", 100,0,200.0);

  plots->fTrackDeltaD0vsEta = result->AddHist2D("Deltad0vsEta","d0vsEta", "eta","d0 - d0_{truth}", 100, -3,3, 100,-0.1,0.1);

  plots->fTrackDeltaThetavsEta = result->AddHist2D("DeltaThetavsEta","DeltaThetavsEta", "eta","theta  - theta_{truth}", 100, -3,3, 100,-0.01,0.01);
  plots->fTrackDeltaZ0vsEta = result->AddHist2D("DeltaZ0vsEta","DeltaZ0vsEta", "eta","z0 - z0_{truth}", 100, -3,3, 100,-1,1);
  plots->fTrackDeltaPhivsEta = result->AddHist2D("DeltaPhivsEta","DeltaPhivsEta", "eta","phi - phi_{truth}", 100, -3,3, 100,-0.01,0.01);
  plots->fTrackDeltaPtvsEta = result->AddHist2D("DeltapTvsEta","DeltaPtvsEta", "eta","(pt - pt_{truth})/pt_{truth}", 100, -3,3, 100,-0.05,0.05);

}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchOriginalTrack = treeReader->UseBranch("OriginalTrack");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;



  TLorentzVector momentum;

  Bool_t skip;

  Long64_t entry;

  Int_t i, j, pdgCode;

  // Loop over all events
  for(entry = 0; entry < 3; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

cout << "Event "<< entry << endl;
    // Loop over all tracks in event
    for(i = 0; i < branchTrack->GetEntriesFast(); ++i)
    {
      Track* track = (Track*) branchTrack->At(i);
      Track* particle = (Track*) track->Particle.GetObject();
      
      float eta= TMath::Log( TMath::Tan(particle->trkPar[3]/2));
      float pt = particle->Charge/(particle->trkPar[4]*cosh(eta));
 
      plots->fOldTrackD0 ->Fill(particle->trkPar[0]);
      plots->fOldTrackZ0 ->Fill(particle->trkPar[1]);
      plots->fOldTrackPhi->Fill(particle->trkPar[2]);
      plots->fOldTrackEta->Fill(eta);
      plots->fOldTrackPt ->Fill(pt);
      
      float d0Err = TMath::Sqrt(track->trkCov[0]);
      plots->fTrackD0Sig ->Fill(track->trkPar[0]/d0Err);
      plots->fTrackD0 ->Fill(track->trkPar[0]);
      plots->fTrackZ0 ->Fill(track->trkPar[1]);
      plots->fTrackPhi->Fill(track->trkPar[2]);
      plots->fTrackEta->Fill(track->P4().Eta());
      plots->fTrackPt ->Fill(track->P4().Pt()); 

      plots->fTrackDeltaD0vsEta->Fill( track->P4().Eta(),  track->trkPar[0]-particle->trkPar[0]);
      plots->fTrackDeltaZ0vsEta->Fill( track->P4().Eta(),  track->trkPar[1]-particle->trkPar[1]);
      plots->fTrackDeltaPhivsEta->Fill( track->P4().Eta(),  track->trkPar[2]-particle->trkPar[2]);
      plots->fTrackDeltaThetavsEta->Fill( track->P4().Eta(),  track->trkPar[3]-particle->trkPar[3]);
      plots->fTrackDeltaPtvsEta->Fill( track->P4().Eta(),  (track->P4().Pt()-particle->P4().Pt())/ particle->P4().Pt());


    }

    // cout << "--  New event -- " << endl;

  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void ExampleTrack(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
