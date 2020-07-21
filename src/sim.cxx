//!sim.cxx    ~~~    Modified by V. Verkest - June 2020
//!Isaac Mooney, WSU - July 2019
//!This file runs the (pp) analysis on the STAR-tuned Pythia6 and Pythia6+Geant for the jet mass project.
//!It takes in the Picos, performs selections, clusters particles, performs selections on the resulting jets,
//!applies the Soft Drop grooming procedure to a copy of the jet population, fills jet trees, fills responses, and writes to files.
//!It also produces histograms for use in the unfolding closure test, and responses for use in the propagation of systematic uncertainty variations through the unfolding.

#include "params.hh"
#include "funcs.hh"
#include "TStarJetPicoDefinitions.h"

using namespace fastjet;
using namespace std;
using namespace Analysis;

//! Command line arguments: ( Defaults
//! Defined for debugging in main )
//! [0]: output directory
//! [1]: output filename
//! [2]: flag: ch/ch+ne jets. Options: "ch" (full = 0), or "full" (full = 1).
//! [3]: flag: require match? Options: "nomatch" (match = 0), or "match" (match = 1).
//! [4]: input data: can be a single .root or a .txt or .list of root files - should always be last argument

int main (int argc, const char ** argv) {
  TH1::SetDefaultSumw2( );    TH2::SetDefaultSumw2( );    TH3::SetDefaultSumw2( );
  
  // Defaults
  std::string outputDir = "out/"; // directory where everything will be saved
  std::string outFileName = "test.root"; //output file
  std::string chainList = "lists/simlist.txt"; // input file: can be .root, .txt, .list
  double radius = 0.4; //jet radius parameter; input value can range from 0.1 to 9.9.
  bool full = 1;  //full = 1 => ch+ne; full = 0 => ch only.
  bool match = 1; //match = 0 => no match between Pythia & Pythia+Geant events.
  
  // Now check to see if we were given modifying arguments
  switch ( argc ) {
  case 1: // Default case
    __OUT("Using Default Settings");
    break;
  case 7: { // Custom case
    __OUT("Using Custom Settings");
    std::vector<std::string> arguments( argv+1, argv+argc );
    
    // Set non-default values
    outputDir         = arguments[0];
    outFileName       = arguments[1];
    radius            = radius_str_to_double (arguments[2]);
    if (arguments[3] == "ch") {full = 0;} else {full = 1;}
    if (arguments[4] == "matched") {match = 1;} else if (arguments[4] == "unmatched") {match = 0;} else {cerr << "Not a valid flag!" << endl; exit(1);}
    chainList         = arguments[5];
    
    //Printing settings:
    cout << "Outputting to: " << (outputDir+outFileName).c_str() << "\nSettings:\nR = " << radius << " " <<  arguments[3] << " jets;\n match pythia and geant? " << match << ";\n input file: " << chainList << "\n";
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
    return -1;
    break;
  }
  }
  
  //  Initialize readers and provide chains
  TStarJetPicoReader* P6Reader = new TStarJetPicoReader();
  TChain* P6Chain = new TChain( "JetTreeMc" ); // PURE PYTHIA (particle)
  TStarJetPicoReader* GEANTReader = new TStarJetPicoReader();
  TChain* GEANTChain = new TChain( "JetTree" ); // CORRESPONDING GEANT (detector)
  
  // Check to see if the input is a .root file or a .txt
  bool inputIsRoot = Analysis::HasEnding( chainList.c_str(), ".root" );
  bool inputIsTxt  = Analysis::HasEnding( chainList.c_str(), ".txt"  );
  bool inputIsList = Analysis::HasEnding( chainList.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot ) { P6Chain->Add( chainList.c_str()); GEANTChain->Add( chainList.c_str());}
  else if ( inputIsTxt || inputIsList)  { P6Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str()); GEANTChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }
  
  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
  fout->cd();
  
  TString pythiaFilename, geantFilename;  //define relevant data structures
  TStarJetPicoEventHeader* p_header; TStarJetPicoEventHeader* g_header;
  TStarJetPicoEvent* p_event; TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector>* p_container; TStarJetVectorContainer<TStarJetVector>* g_container;
  TStarJetVector *p_sv = 0; TStarJetVector *g_sv = 0;
  
  int p_EventID, g_EventID;  //defining local containers to be linked to tree branches
  double p_n_jets, p_wt, g_n_jets, g_wt;
  vector<vector<double> > p_conspT, g_conspT;
  vector<double> p_jetMult, p_Pt, p_Eta, p_Phi, p_M, p_E, p_ch_e_frac, p_zg, p_rg, p_mg, p_ptg, p_mcd;
  vector<double> g_jetMult, g_Pt, g_Eta, g_Phi, g_M, g_E, g_ch_e_frac, g_zg, g_rg, g_mg, g_ptg, g_mcd;
  

  TTree *eventTree = new TTree("event","event");  //tree to hold jet and constituent quantites
  eventTree->Branch("p_n_jets", &p_n_jets);  eventTree->Branch("p_conspT",&p_conspT);  eventTree->Branch("p_jetMult",&p_jetMult);
  eventTree->Branch("p_Pt", &p_Pt); eventTree->Branch("p_Eta",&p_Eta); eventTree->Branch("p_Phi",&p_Phi); eventTree->Branch("p_M",&p_M); eventTree->Branch("p_E",&p_E);
  eventTree->Branch("p_ch_e_frac", &p_ch_e_frac);
  eventTree->Branch("p_zg", &p_zg); eventTree->Branch("p_rg", &p_rg); eventTree->Branch("p_mg", &p_mg); eventTree->Branch("p_ptg",&p_ptg);
  eventTree->Branch("p_mcd",&p_mcd);  eventTree->Branch("p_weight", &p_wt); eventTree->Branch("p_EventID", &p_EventID);
  
  eventTree->Branch("g_n_jets", &g_n_jets);  eventTree->Branch("g_conspT",&g_conspT);  eventTree->Branch("g_jetMult",&g_jetMult);
  eventTree->Branch("g_Pt", &g_Pt); eventTree->Branch("g_Eta",&g_Eta); eventTree->Branch("g_Phi",&g_Phi); eventTree->Branch("g_M",&g_M); eventTree->Branch("g_E",&g_E);
  eventTree->Branch("g_ch_e_frac", &g_ch_e_frac);
  eventTree->Branch("g_zg", &g_zg); eventTree->Branch("g_rg", &g_rg); eventTree->Branch("g_mg", &g_mg); eventTree->Branch("g_ptg",&g_ptg);
  eventTree->Branch("g_mcd",&g_mcd);  eventTree->Branch("g_weight", &g_wt); eventTree->Branch("g_EventID", &g_EventID);
  
  //hist for statistical error correction - for matched events, but unmatched jets
  TH1D* pt_gen_match_plus_miss = new TH1D("pt_gen_match_plus_miss","",15,5,80);
  
  //1D hists for closure test
  TH1D* sampleA_pt_gen = new TH1D("sampleA_pt_gen","",15,5,80);   TH1D* sampleA_pt_det = new TH1D("sampleA_pt_det","",9,15,60);
    
  TH1D* sampleB_pt_gen = new TH1D("sampleB_pt_gen","",15,5,80);   TH1D* sampleB_pt_det = new TH1D("sampleB_pt_det","",9,15,60);

  //~~~counts histograms for bin_drop, later~~~//
    
  //1D hists for closure test
  TH1D* sampleA_pt_gen_counts = new TH1D("sampleA_pt_gen_counts","",15,5,80);
  TH1D* sampleA_pt_det_counts = new TH1D("sampleA_pt_det_counts","",9,15,60);

  TH1D* sampleB_pt_gen_counts = new TH1D("sampleB_pt_gen_counts","",15,5,80);
  TH1D* sampleB_pt_det_counts = new TH1D("sampleB_pt_det_counts","",9,15,60);
  //~~~~~~//
    
  // 1D responses
  RooUnfoldResponse *pt_response = new RooUnfoldResponse(55,4.5,59.5,55,4.5,59.5,"pt_response","");
  //responses for training in the closure test
  RooUnfoldResponse *sampleA_pt_response = new RooUnfoldResponse(55,4.5,59.5,55,4.5,59.5,"sampleA_pt_response","");
  //responses for validation in the closure test (unfolding pseudo-data with response constructed from same sample)
  RooUnfoldResponse *sampleB_pt_response = new RooUnfoldResponse(55,4.5,59.5,55,4.5,59.5,"sampleB_pt_response","");

  TH1D *hMisses = new TH1D( "hMisses","Missses;missing part-level leading jet p_{T} (GeV)",55,4.5,59.5);
  TH1D *hFakes = new TH1D( "hFakes","Fakes;fake det-level leading jet p_{T} (GeV)",55,4.5,59.5);

  //vectors of responses & hists for easy writing to file later
  std::vector<RooUnfoldResponse*> res = {pt_response};
  std::vector<RooUnfoldResponse*> sampleA_res = {sampleA_pt_response};
  std::vector<RooUnfoldResponse*> sampleB_res = {sampleB_pt_response};
  
  std::vector<TH1D*> sampleA_h1Ds = {sampleA_pt_gen,sampleA_pt_det};
  std::vector<TH1D*> sampleB_h1Ds = {sampleB_pt_gen,sampleB_pt_det};    
    
  //defining the algorithm and radius parameter for clustering jets
  JetDefinition jet_def(antikt_algorithm, radius/*R*/);
  
  //SELECTORS
  // Constituent selectors
  Selector select_track_rap = fastjet::SelectorAbsRapMax(max_track_rap);
  Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
  Selector select_loptmax   = fastjet::SelectorPtMax( partMaxPt );
  Selector spart = select_track_rap * select_lopt * select_loptmax;
  
  // Jet candidate selectors
  Selector select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
  Selector select_det_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
  Selector select_gen_jet_pt_min = fastjet::SelectorPtMin( jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector select_det_jet_m_min = fastjet::SelectorMassMin( mass_min );
  Selector select_gen_jet_m_min = fastjet::SelectorMassMin( 0.0 );
  
  Selector sjet_gen = select_jet_rap && select_gen_jet_pt_min && select_jet_pt_max && select_gen_jet_m_min;
  Selector sjet_det = select_jet_rap && select_det_jet_pt_min && select_jet_pt_max && select_det_jet_m_min;
  
  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles, p_JetsInitial, g_JetsInitial;
  int p_n_accepted = 0, g_n_accepted = 0; int p_NJets = 0, g_NJets = 0;
  int counter_debug = 0;
  double mc_weight = -1;
  
  double hc = 0.9999; //to be varied in the systematic uncertainty variation
  const int nSources = /*7*/1; //includes the nominal settings as a "systematic".//ipse fixit
  for (int iSyst = 0; iSyst < nSources; ++ iSyst) {
    if (iSyst == 0) {cout << endl << "RUNNING WITH NOMINAL SETTINGS!" << endl << endl;}
    if (iSyst == 1) {cout << endl << "RUNNING WITH INCREASED TOWER SCALE!" << endl << endl;}
    if (iSyst == 2) {cout << endl << "RUNNING WITH DECREASED TRACKING EFFICIENCY!" << endl << endl;}
    if (iSyst == 3) {cout << endl << "RUNNING WITH 50% HADRONIC CORRECTION!" << endl << endl;}
    if (iSyst == 4) {cout << endl << "RUNNING WITH SMEARED DETECTOR SPECTRUM!" << endl << endl;}
    if (iSyst == 5) {cout << endl << "RUNNING WITH SMEARED GENERATOR SPECTRUM!" << endl << endl;}
    if (iSyst == 6) {cout << endl << "RUNNING WITH SMEARED GENERATOR ~MASS~ SPECTRUM!" << endl << endl;}
    
    p_NJets = 0; g_NJets = 0; p_n_accepted = 0; g_n_accepted = 0; counter_debug = 0;
    //set parameters back to their nominal values after the previous iteration changed them.
    hc = 0.9999;
    
    //change the nominal values
    if (iSyst == 3) { //this means the systematic we're examining is the hadronic correction. Set it to 50%
      hc = 0.5;
    }
   
    //initialize both readers
    InitReader(P6Reader, P6Chain, nEvents, "All", truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, hc, false, sim_badTowers, sim_bad_run_list);
    InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, hc, false, /*sim_badTowers, sim_bad_run_list);*/det_badTowers, dat_bad_run_list);
    
    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
    
    for (int event = 0; event < P6Chain->GetEntries(); ++ event) {
      P6Reader->ReadEvent(event);
      GEANTReader->ReadEvent(event);
      
      g_EventID = GEANTReader->GetNOfCurrentEvent();
      p_EventID = P6Reader->GetNOfCurrentEvent();
      //NOTE: even # events will be used for training (response), odd will be used for validation (pseudo-data)
      
      //clearing vectors; initializing variables to -9999
      mc_weight = -9999;
      p_n_jets = -9999; p_wt = -9999;
      p_conspT.clear();
      p_jetMult.clear();
      p_Pt.clear(); p_Eta.clear(); p_Phi.clear(); p_M.clear(); p_E.clear();
      p_ch_e_frac.clear();
      p_zg.clear(); p_rg.clear(); p_mg.clear(); p_ptg.clear();
      p_mcd.clear();
      
      g_n_jets = -9999; g_wt = -9999;
      g_conspT.clear();
      g_jetMult.clear();
      g_Pt.clear(); g_Eta.clear(); g_Phi.clear(); g_M.clear(); g_E.clear();
      g_ch_e_frac.clear();
      g_zg.clear(); g_rg.clear(); g_mg.clear(); g_ptg.clear();
      g_mcd.clear();
      
      p_Particles.clear(); g_Particles.clear();
      p_JetsInitial.clear(); g_JetsInitial.clear();
      
      //if we require matching, must have an event in both P+G and P
      if ( match && GEANTReader->ReadEvent(p_EventID) != 1 ) {
	// cout << "no corresponding geant event...skipping event " << p_EventID <<endl;
	continue;//goes to the next event
      }
      //sanity check:
      if ( match && (p_EventID != g_EventID) ) { cerr << "ERROR: READING DIFFERENT EVENTS. EXITING." << endl; exit(1);}
      
      if (P6Reader->ReadEvent(p_EventID) == 1) {p_n_accepted++;} // == 1 => event loaded and passed selections
      if (GEANTReader->ReadEvent(g_EventID) == 1) {g_n_accepted++;}
            
      P6Reader->PrintStatus(10); GEANTReader->PrintStatus(10);     // Print out reader status every 10 seconds
      
      //filling the data structures that were defined before the event loop:
      p_event = P6Reader->GetEvent();
      p_header = p_event->GetHeader();
      g_event = GEANTReader->GetEvent();
      g_header = g_event->GetHeader();
      
      if (!(g_header->HasTriggerId(500205) || g_header->HasTriggerId(500215))) {continue;}   //  ONLY SELECT HT TRIGGER EVENTS

      TList *SelectedTowers = GEANTReader->GetListOfSelectedTowers();
      int nTowers = CountTowers( SelectedTowers );

      int trigTowId;
      TStarJetPicoTriggerInfo *trig;
      TStarJetPicoTower *tow, *triggerTower;
      double trigTowEt = 0.0;
      std::vector<int> trigTowers;
      for ( int i=0; i<g_event->GetTrigObjs()->GetEntries(); ++i ) {
	trig = (TStarJetPicoTriggerInfo *)g_event->GetTrigObj(i);
	if ( trig->isBHT2() && UseTriggerTower( trig->GetId()) ) { trigTowers.push_back( trig->GetId() ); }
      }
      sort(trigTowers.begin(), trigTowers.end());

      int nmatched = 0;
      for (int i=0; i<nTowers; ++i){				// loop throught selected towers in event
	tow = (TStarJetPicoTower*)SelectedTowers->At(i);
	if ( !(tow->GetEt()>=5.4 && count(trigTowers.begin(), trigTowers.end(), tow->GetId()))) { // min 5.4 GeV tower and must be in list of HT towers
	  
	}
      }



      
      p_container = P6Reader->GetOutputContainer();
      g_container = GEANTReader->GetOutputContainer();
    
      pythiaFilename =  P6Reader->GetInputChain()->GetCurrentFile()->GetName();
      geantFilename =  GEANTReader->GetInputChain()->GetCurrentFile()->GetName();
      
      if (match && (pythiaFilename != geantFilename)) {std::cerr << "FILES DON'T MATCH! EXITING." << std::endl; exit(1);}
 
      p_wt = LookupRun15Xsec( pythiaFilename );
      g_wt = LookupRun15Xsec( geantFilename );
      if (match && (p_wt != g_wt)) {std::cerr << "WEIGHTS DON'T MATCH! EXITING." << std::endl; exit(1);}
      mc_weight = p_wt; //arbitrarily setting it to pythia's but can be either
      
      if (iSyst == 1) {//varying the gain of the towers
	for (int i = 0; i < g_container->GetEntries(); ++ i) {
	  g_sv = g_container->Get(i);
	  if (!(g_sv->IsCharged())) {
	    //cout << "DEBUG: Pre-change: " << g_sv->E() << " " << g_sv->Eta() << " " << g_sv->Phi() << " " << g_sv->M() << endl;
	    double Enew = 1.038*g_sv->E();
	    g_sv->SetE(Enew);
	    //g_sv->SetPtEtaPhiM(sqrt(Etnew*Etnew - g_sv->M()*g_sv->M()), g_sv->Eta(), g_sv->Phi(), g_sv->M());
	    //cout << "DEBUG: Post-change: " << g_sv->E() << " " << g_sv->Eta() << " " << g_sv->Phi() << " " << g_sv->M() << endl;
	    //cout << "DEBUG: Percent diff: " << ((g_sv->E()/(double)(Enew/(double)1.038)) - 1)*100 << "%" << endl;
	  }
	}
      }
      
      // converts TStarJetVectors to PseudoJets, carefully assigning the proper particle mass in either case
      GatherParticles ( p_container, p_sv, p_Particles, full, 1); //Pythia; full = 0 => charged-only, 1 => ch+ne
      GatherParticles ( g_container, g_sv, g_Particles, full, 0); //GEANT
      
      if (iSyst == 2) {//varying the tracking efficiency randomly by 4%
	double effic_num;
	for (int i = 0; i < g_Particles.size(); ++ i) {
	  if (g_Particles[i].user_index() != 0) {
	    effic_num = gRandom->Uniform(0.0, 1.0);
	    if (effic_num > 0.96) {
	      g_Particles.erase(g_Particles.begin() + i);
	      i --; //need to account for the shrinking of the list.
	    }
	  }
	}
      }
      
      // applying particle-level cuts
      vector<PseudoJet> p_cut_Particles = spart(p_Particles);
      vector<PseudoJet> g_cut_Particles = spart(g_Particles);
      
      //Clustering jets
      ClusterSequence p_Cluster(p_cut_Particles, jet_def);
      ClusterSequence g_Cluster(g_cut_Particles, jet_def);
      p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets()));
      g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));
      
      vector<PseudoJet> p_Jets;
      vector<PseudoJet> g_Jets;
      
      //Implementing a neutral energy fraction cut (of 90% currently) on inclusive det-level jets
      p_Jets = p_JetsInitial; //just passing intermediate -> final vector (no NEF selection on Pythia)
      ApplyNEFSelection(g_JetsInitial, g_Jets);
      
      //at this point, all criteria have been met (except checking if the event is bad in DiscardEvent() below)
      //so we can start filling jet information e.g.:
      p_n_jets = p_Jets.size();
      g_n_jets = g_Jets.size();
      
      if (DiscardEmbedEvent(pythiaFilename, p_Jets, g_Jets)) { counter_debug ++; continue; }

      if ( p_Jets.size()==0 || g_Jets.size()==0 ) { continue; }
      
      //for calculating charged energy fraction of the jets
      vector<double> pch_e_frac, gch_e_frac;
      double ch_e = 0; double tot_e = 0;//names are misnomers here since we use pT, not E.
      vector<PseudoJet> cons = p_Jets[0].constituents();
      if ( p_Jets.size()!=0 ) {
	for (int j = 0; j < cons.size(); ++ j) {
	  if (cons[j].user_index() != 0) {ch_e += cons[j].pt();}
	  tot_e += cons[j].pt();
	}
	pch_e_frac.push_back(ch_e/(double)tot_e);
      }
      ch_e = 0;  tot_e = 0;//names are misnomers here since we use pT, not E.
      cons = g_Jets[0].constituents();
      if ( g_Jets.size()!=0 ) {
	for (int j = 0; j < cons.size(); ++ j) {
	  if (cons[j].user_index() != 0) {ch_e += cons[j].pt();}
	  tot_e += cons[j].pt();
	}
	gch_e_frac.push_back(ch_e/(double)tot_e);
      }

      //We have two vectors to be filled with matched jets. If they aren't, when looping over pythia jets, we have misses. Same goes when looping over geant jets with fakes. And for matches, we just fill with however many entries there are in the matched vectors.
      //MatchJets returns a vector of pairs of indices (i,j). The first entry is the position of the jet to match, the second its match's position, the third the position of the next jet to match, the fourth its match's position, etc.
      //FakesandMisses returns a vector of indices (i) corresponding to the indices of misses or fakes from the original candidate vector.
      if (match) {
	//before the actual matching, fill the spectrum that will be used for statistical error correction later
	if ( p_Jets.size()!=0 ) { pt_gen_match_plus_miss->Fill(p_Jets[0].pt(), mc_weight); }
	
	//before the actual matching, fill the spectra for validation in the closure test
	//(we only require there to be a matching Geant & Pythia event so our event reconstruction efficiency isn't folded unnecessarily into the closure)
	if (p_EventID % 2 == 0) { //even events are sampleA
	  if ( p_Jets.size()!=0 ) { sampleA_pt_gen->Fill(p_Jets[0].pt(),mc_weight); sampleA_pt_gen_counts->Fill(p_Jets[0].pt()); }
	  
	  if ( g_Jets.size()!=0 ) { sampleA_pt_det->Fill(g_Jets[0].pt(),mc_weight); sampleA_pt_det_counts->Fill(g_Jets[0].pt()); }
	}//end sampleA
	if (p_EventID % 2 != 0) { //odd events are sampleB
	  if ( p_Jets.size()!=0 ) { sampleB_pt_gen->Fill(p_Jets[0].pt(),mc_weight); sampleB_pt_gen_counts->Fill(p_Jets[0].pt()); }
	  
	  if ( g_Jets.size()!=0 ) { sampleB_pt_det->Fill(g_Jets[0].pt(),mc_weight); sampleB_pt_det_counts->Fill(g_Jets[0].pt()); }
	}//end sampleB pseudo-data filling
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MATCHING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	std::vector<fastjet::PseudoJet> g_matches, p_matches, g_matches_for_fakes, p_matches_for_fakes, fakes,
	  misses, g_sd_matches, p_sd_matches, sd_fakes, sd_misses;
	std::vector<int> match_indices, miss_indices, fake_indices;

	//matches & misses
	if (p_Jets.size() != 0) {
	  g_matches.clear(); p_matches.clear(); misses.clear();
	  match_indices.clear(); miss_indices.clear();
                
	  match_indices = MatchJets(g_Jets, p_Jets, g_matches, p_matches); //find matches
                
	  if (g_matches.size() != p_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}

	  if (g_matches.size() == 0 || p_matches.size() ==0) {continue;}//require a matched jet!
	  
	  if (g_matches.size() < p_Jets.size()) { //then we have misses
	    miss_indices = FakesandMisses(p_matches, p_Jets, misses); //find misses
	  }
	}
        
	//fakes
	if (g_Jets.size() != 0) {
	  //clear the vectors to be used for determination of fakes (jets we find in Geant that don't have a match in Pythia)
	  fakes.clear(); fake_indices.clear();
	  g_matches_for_fakes.clear(); p_matches_for_fakes.clear();
		
	  MatchJets(p_Jets, g_Jets, p_matches_for_fakes, g_matches_for_fakes);
                
	  if (g_matches_for_fakes.size() != p_matches_for_fakes.size()) {std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}
                
	  if (p_matches_for_fakes.size() < g_Jets.size()) { //then we have fakes
	    fake_indices = FakesandMisses(g_matches_for_fakes, g_Jets, fakes);
	  }
	}
        
	double prior_adjust = 0, prior_adjust_g = 0;
	if ( misses.size()!=0 ) {
	  //determine on a per-jet basis the pT and M smearing for the systematic prior variation
	  
	  pt_response->Miss(misses[0].pt(), mc_weight);
	  hMisses->Fill(misses[0].pt(), mc_weight);

	  if (match && p_EventID % 2 == 0) { //throughout this section, checking if (match) is unnecessary because of the overall conditional
	    sampleA_pt_response->Miss(misses[0].pt(), mc_weight);
	  }
	  //closure sampleB responses
	  if (match && p_EventID % 2 != 0) {
	    sampleB_pt_response->Miss(misses[0].pt(), mc_weight);
	  }
	}
	
	pt_response->Fill(p_matches[0].pt(), g_matches[0].pt(), mc_weight);
	//closure sampleA responses
	if (match && p_EventID % 2 == 0) {
	  sampleA_pt_response->Fill(g_matches[0].pt(), p_matches[0].pt(), mc_weight);
	}
	//closure sampleB responses
	if (match && p_EventID % 2 != 0) {
	  sampleB_pt_response->Fill(g_matches[0].pt(), p_matches[0].pt(), mc_weight);
	}
	      
	//(matched) TREES:
	//ungroomed
	p_Pt.push_back(p_matches[0].pt());
	p_M.push_back(p_matches[0].m());
	p_Eta.push_back(p_matches[0].eta());
	p_Phi.push_back(p_matches[0].phi());
	p_E.push_back(p_matches[0].e());
	//ungroomed
	g_Pt.push_back(g_matches[0].pt());
	g_M.push_back(g_matches[0].m());
	g_Eta.push_back(g_matches[0].eta());
	g_Phi.push_back(g_matches[0].phi());
	g_E.push_back(g_matches[0].e());
	//assigning vectors of values we calculated earlier to the tree
	p_ch_e_frac.push_back(pch_e_frac[match_indices[1]]);
	g_ch_e_frac.push_back(gch_e_frac[match_indices[1]]);
	// }//for loop over matches
            
	if ( fakes.size()!=0 ) {
	  hFakes->Fill(fakes[0].pt(), mc_weight);
	  pt_response->Fake(fakes[0].pt(), mc_weight);
	  //closure sampleA responses
	  if (match && p_EventID % 2 == 0) {
	    sampleA_pt_response->Fake(fakes[0].pt(), mc_weight);
	  }
	  //closure sampleB responses
	  if (match && p_EventID % 2 != 0) {
	    sampleB_pt_response->Fake(fakes[0].pt(), mc_weight);
	  }
	}	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      }//matching-required conditional


	
      if (p_Jets.size() != 0 && iSyst == 0) { eventTree->Fill(); } //when !match, will fill sometimes with empty geant vectors in the geant branches
        
      p_NJets += p_Jets.size(); g_NJets += g_Jets.size(); // add jets to total
    }//event loop
    
    //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
    
    fout->cd();
    
    cout << endl << endl << "For systematic variation " << iSyst << endl
	 << "Of " << p_n_accepted << " pythia events and " << g_n_accepted << " geant events" << endl
	 << p_NJets << " gen jets have been found" << endl
	 << g_NJets << " det jets have been found" << endl << endl
	 << "Discarded " << counter_debug << " events on grounds of the found jets being too much higher than the pT-hat range" << endl;
    
    if (iSyst == 0) {//nominal case, so write closure and regular responses
      //trees
      eventTree->Write();
      
      if (match) {
	//hists
	pt_gen_match_plus_miss->Write();
	hFakes->Write();
	hMisses->Write();
	
	for (int i = 0; i < sampleA_h1Ds.size(); ++ i) { sampleA_h1Ds[i]->Write(); }
	for (int i = 0; i < sampleB_h1Ds.size(); ++ i) { sampleB_h1Ds[i]->Write(); }

	//responses
	for (int i = 0; i < res.size(); ++ i) { res[i]->Write(); }
	for (int i = 0; i < sampleA_res.size(); ++ i) { sampleA_res[i]->Write(); }
	for (int i = 0; i < sampleB_res.size(); ++ i) { sampleB_res[i]->Write(); }
      }//match conditional
    }//iSyst==0 conditional
    
  }//systematic variation loop

  cout << endl << "Writing to:  " << fout->GetName() << endl;
  fout->Close();
    
  return 0;
}//main
