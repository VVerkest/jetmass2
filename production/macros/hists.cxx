//Isaac Mooney, WSU, June 2019
//This file takes in a root file, pulls the tree(s) out, fills histograms with the entries, and writes to an output root file

#include <ctime>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TRandom.h>

using namespace std;

//self explanatory: takes branches in a tree, "treestr", in file f, and fills histograms with them
void TreetoHist (TFile *f, string treestr, vector<string> branchnames, vector<TH1D*> h1Ds, vector<TH2D*> h2Ds, vector<TH3D*> h3Ds, bool hadronic) {  
  cout << "STARTING TREETOHIST" << endl;
  
  //initializing the variables that will be filled by the values in the branches later
  vector<double> *Pt = 0; vector<double> *Eta = 0; vector<double> *Phi = 0;
  vector<double> *M = 0; vector<double> *zg = 0; vector<double> *rg = 0; vector<double> *mg = 0;
  vector<int> *qvg = 0;
  vector<vector<double> > *PID = 0; vector<vector<double> > *consPt = 0;
  
  double weight = -1;
  
  cout << "INITIALIZED BRANCH VARIABLES!" << endl;
  
  //getting the tree and linking the branches with the variables
  TTree *t = (TTree*) f->Get(treestr.c_str());
  t->SetBranchAddress(branchnames[0].c_str(),&Pt);
  t->SetBranchAddress(branchnames[1].c_str(),&Eta);
  t->SetBranchAddress(branchnames[2].c_str(),&Phi);
  t->SetBranchAddress(branchnames[3].c_str(),&M);
  t->SetBranchAddress(branchnames[4].c_str(),&mg);
  t->SetBranchAddress(branchnames[5].c_str(),&zg);
  t->SetBranchAddress(branchnames[6].c_str(),&rg);

  if (hadronic) {
    t->SetBranchAddress("qvg",&qvg);
    t->SetBranchAddress("consPID",&PID);
    t->SetBranchAddress("conspT",&consPt);
  }
  
  t->SetBranchAddress("mcweight", &weight);

  cout << "ASSIGNED VARIABLES TO BRANCHES!" << endl;
  
  vector<double> JEF_pi(15),JEF_K(15),JEF_p(15),JEF_other(15),JEF_tot(15);
  
  cout << "BOOKED SPACE FOR JEF VECTORS" << endl;
  
  cout << /*(*/"RUNNING OVER TREE ResultTree"/*+treestr+*/ << "! Entries: "/*).c_str()*/ << t->GetEntries() << endl;
  const clock_t begin_time = clock(); //timing - for debugging and for fun
  for (unsigned i = 0; i < t->GetEntries(); ++ i) { //"event" loop
    //cout << "IN EVENT LOOP" << endl;
    if (i % 1000 == 0 && i != 0) { //can change this to a frequency of your preference (for real data I use 1e5 or 1e6)
      cout << "Still chuggin. On event " << i << endl;
      cout << "Total time passed: " << fixed << setprecision(5) << double(clock() - begin_time) /(double) CLOCKS_PER_SEC << " secs" << endl;
    }
    //    cout << "GETTING ENTRY " << i << endl;
    t->GetEntry(i);
    cout << "GOT ENTRY " << i << endl;
    
    //filling "event" observables
    //
    //looping over jets and filling histograms
    for (unsigned j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size
      /*      //cout << "IN JET LOOP" << endl;
      if (hadronic) {//don't have the following branches coded in the parton-tree
	for (unsigned k = 0; k < PID->at(j).size(); ++ k) { //constituent loop
	  //	  cout << "IN CONSTITUENT LOOP" << endl;
	  //~~~this block will handle the pi/K/p jet energy fraction calculation~~~//
	  if (fabs(PID->at(j).at(k)) == (double) 211 || PID->at(j).at(k) == (double) 111) {
	    //cout << "DEBUG: PIONS: PID = " << PID->at(j).at(k) << " should be |211| or 111" << endl;
	    //for each jet pT we have a running count of the JEF of pions. We therefore index the vector of pions depending on the jet's pT. E.g. for 5 GeV bins starting at 5 GeV, a 20 GeV jet is in the 3rd bin (indexing from 0) or: (20 - 5)/5. 
	    JEF_pi[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	    //cout << "DEBUG: PIONS: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	  }
	  else if (fabs(PID->at(j).at(k)) == (double) 321 || PID->at(j).at(k) == (double) 311 || PID->at(j).at(k) == (double) 310 || PID->at(j).at(k) == (double) 130) {
	    //cout << "DEBUG: KAONS: PID = " << PID->at(j).at(k) << " should be |321|, 311, 310, or 130" << endl;
	    JEF_K[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	    //cout << "DEBUG: KAONS: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	  }
	  else if (fabs(PID->at(j).at(k)) == (double) 2212) {
	    //cout << "DEBUG: PROTONS: PID = " << PID->at(j).at(k) << " should be |2212|" << endl;
	    JEF_p[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	    //cout << "DEBUG: PROTONS: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	  }
	  else {
	    //cout << "DEBUG: OTHERS: PID = " << PID->at(j).at(k) << " should be anything but |211|, 111, |321|, 311, 310, 130, or |2212|" << endl;
	    JEF_other[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	    //cout << "DEBUG: OTHERS: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	  }
	  JEF_tot[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	  //cout << "DEBUG: TOTAL: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	  //~~~          ~~~//
	}//for loop over constituents
      }//hadronic conditional
      *///2Ds!                                                                                                                                           
      h2Ds[0]->Fill(M->at(j), Pt->at(j), weight);
      //quark v. gluon jets
      if (hadronic && qvg->at(j) == 0) {//q jet!                                                                                                                           
        h2Ds[1]->Fill(M->at(j),Pt->at(j),weight);
      }
      else if (hadronic && qvg->at(j) == 1) {//g jet!                                                                                                                      
        h2Ds[2]->Fill(M->at(j),Pt->at(j),weight);
      }
      else {//neither (or no match to hard parton)!
        h2Ds[3]->Fill(M->at(j),Pt->at(j),weight);
      }

      h2Ds[4]->Fill(Eta->at(j), Pt->at(j), weight);
      h2Ds[5]->Fill(Phi->at(j), Pt->at(j), weight);
      if (zg->at(j) >= 0.1) { //otherwise, we tag and drop the groomed jet
	h2Ds[6]->Fill(mg->at(j), Pt->at(j), weight);
	h2Ds[7]->Fill(zg->at(j), Pt->at(j), weight);
	h2Ds[8]->Fill(rg->at(j), Pt->at(j), weight);
      
	//3Ds!
	h3Ds[0]->Fill(M->at(j),mg->at(j),Pt->at(j),weight);
	h3Ds[1]->Fill(M->at(j),zg->at(j),Pt->at(j),weight);
	h3Ds[2]->Fill(M->at(j),rg->at(j),Pt->at(j),weight);
	h3Ds[3]->Fill(mg->at(j),zg->at(j),Pt->at(j),weight);
	h3Ds[4]->Fill(mg->at(j),rg->at(j),Pt->at(j),weight);
	h3Ds[5]->Fill(zg->at(j),rg->at(j),Pt->at(j),weight);
      } 
    }//!jet loop
    //    cout << "DONE WITH JET LOOP" << endl;
  }//!event loop
  cout << "DONE WITH EVENT LOOP" << endl;
  /*
  if (hadronic) {//normalizing the JEFs
    for (unsigned i = 0; i < JEF_tot.size(); ++ i) {
      cout << "DEBUG: " << JEF_pi[i] << " " << JEF_K[i] << " " << JEF_p[i] << " " << JEF_other[i] << " " << JEF_tot[i] << endl;
      if (JEF_tot[i] != 0) {
	JEF_pi[i] /= (double) JEF_tot[i];
	JEF_K[i] /= (double) JEF_tot[i];
	JEF_p[i] /= (double) JEF_tot[i];
	JEF_other[i] /= (double) JEF_tot[i];
      }
      else { //to avoid dividing by zero when there are no jets in the bin
	JEF_pi[i] = 0;
	JEF_K[i] = 0;
	JEF_p[i] = 0;
	JEF_other[i] = 0;
      }
      h1Ds[0]->Fill(5*i+5,JEF_pi[i]);//re-converts the bin number to a jet pT value for filling the hist
      h1Ds[1]->Fill(5*i+5,JEF_K[i]);
      h1Ds[2]->Fill(5*i+5,JEF_p[i]);
      h1Ds[3]->Fill(5*i+5,JEF_other[i]);
    }
  }
*/
  cout << "HeLP!" << endl;
  //!needs to be outside the event loop; not sure exactly what it does
  t->ResetBranchAddresses();
  cout << "HeLP2!" << endl;
  return;
}

//self explanatory: takes branches in a tree, "treestr", in file f, and fills histograms with them
//this time, requires matches between hadron and parton level so we can make 2D plots with both
void MatchedTreetoHist (TFile *f, string treestr, vector<TH2D*> h2Ds) {  
  //initializing the variables that will be filled by the values in the branches later
 
  vector<double> *PL_pt = 0;
  vector<double> *PL_M = 0;
  vector<double> *PL_mg = 0;
  vector<double> *PL_zg = 0;
  vector<double> *HL_pt = 0;
  vector<double> *HL_M = 0;
  vector<double> *HL_mg = 0;
  vector<double> *HL_zg = 0;
  
  double weight = -1;

  //getting the tree and linking the branches with the variables                                                                                                           
  TTree *t = (TTree*) f->Get(treestr.c_str());
  //parton-level                                                                                                                                
  t->SetBranchAddress("PL_pt_match",&PL_pt);
  t->SetBranchAddress("PL_m_match",&PL_M);
  t->SetBranchAddress("PL_mg_match",&PL_mg);
  t->SetBranchAddress("PL_zg_match",&PL_zg);
  //hadron-level
  t->SetBranchAddress("HL_pt_match",&HL_pt);
  t->SetBranchAddress("HL_m_match",&HL_M);
  t->SetBranchAddress("HL_mg_match",&HL_mg);
  t->SetBranchAddress("HL_zg_match",&HL_zg);
  //simulation needs to be x-section weighted                                                                                                                             
  t->SetBranchAddress("mcweight", &weight);

  cout << ("RUNNING OVER TREE "+treestr+"! Entries: ").c_str() << t->GetEntries() << endl;
  const clock_t begin_time = clock(); //timing - for debugging and for fun
  for (unsigned i = 0; i < t->GetEntries(); ++ i) { //"event" loop
    if (i % 1000 == 0 && i != 0) { //can change this to a frequency of your preference (for real data I use 1e5 or 1e6)
      cout << "Still chuggin. On event " << i << endl;
      cout << "Total time passed: " << fixed << setprecision(5) << double(clock() - begin_time) /(double) CLOCKS_PER_SEC << " secs" << endl;
    }
    t->GetEntry(i);
    
    //filling "event" observables
    //
    //looping over "tracks" and filling histograms
    for (unsigned j = 0; j < PL_pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size 
      //2Ds!                                                                                                                                           
      h2Ds[0]->Fill(PL_M->at(j), HL_pt->at(j), weight);
      //      if (PL_zg->at(j) >= 0.1) { //otherwise, we tag and drop the groomed jet //could also choose the requirement to be at hadron-level or both or either
	h2Ds[1]->Fill(PL_mg->at(j), HL_pt->at(j), weight);
	//      }//CHANGE BACK LATER!!!!!!!!! 
    }//!jet loop
  }//!event loop
  cout << "HELP!" << endl;
  //!needs to be outside the event loop; not sure exactly what it does
  t->ResetBranchAddresses();
  cout << "HELP2!" << endl;
  return;
}


int main (int argc, const char ** argv) {
  //intro
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //basic argument checking.
  if (argc != 4) {
    cerr << "Should be three arguments: output location, output name, input name. Received "
	 << argc << ". Exiting." << endl;
    exit(1);
  }

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  //opening file containing some example trees
  //argv[3] should be the name of the input file
  //string fin_name = (string) argv[3];
  cout << "WHAT ARE YOUUUUUUUUU? " << (string) argv[3] << endl;
  string fin_name; fin_name.assign((string) argv[3]);//strcpy(fin_name,(string) argv[3]);
  TFile *fin = new TFile(fin_name.c_str(),"READ");
  cout << "DEBUG: input file name is " << fin->GetName() << endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hists~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //event observables:
  //
  
  //pT range for particle-level is different from detector-level so define it here for ease of change later if necessary:
  const int ptbins = 15;
  const double ptlow = 5; const double pthigh = 80;
  
  //for variable bin size pT spectrum                                                                                                                       
  const int nBins_pt = 8;
  double edges[nBins_pt + 1] = {5,10,15,20,25,30,40,60,100};
  
  TH1D* dummy1D = new TH1D("dummy1D","",1,0,1);
  TH2D* dummy = new TH2D("dummy","",1,0,1,1,0,1);
  
  //~~~HADRON LEVEL~~~
  //2D correlations with pT
  TH2D* mvpt = new TH2D("mvpt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);
  
  TH2D* mvpt_q = new TH2D("mvpt_q",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);//quark jets
  TH2D* mvpt_g = new TH2D("mvpt_g",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);//gluon jets
  TH2D* mvpt_n = new TH2D("mvpt_n",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);//neither

  
  TH2D* mgvpt = new TH2D("mgvpt",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);
  TH2D* etavpt = new TH2D("etavpt",";#eta;p_{T} [GeV/c]",50,-1,1,ptbins,ptlow,pthigh);
  TH2D* phivpt = new TH2D("phivpt",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,ptbins,ptlow,pthigh);
  TH2D* zgvpt = new TH2D("zgvpt",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,ptbins,ptlow,pthigh);
  TH2D* rgvpt = new TH2D("rgvpt",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,ptbins,ptlow,pthigh);

  //jet energy fraction by species, and delta_R & z as function of species
  TH1D* ptvJEF_pi = new TH1D("ptvJEF_pi","",15,5,80);
  TH1D* ptvJEF_K = new TH1D("ptvJEF_K","",15,5,80);
  TH1D* ptvJEF_p = new TH1D("ptvJEF_p","",15,5,80);
  TH1D* ptvJEF_other = new TH1D("ptvJEF_other","",15,5,80);
  TH2D* deltaRvpt_pi = new TH2D("deltaRvpt_pi","",10,0,0.5,15,5,80);
  TH2D* deltaRvpt_K = new TH2D("deltaRvpt_K","",10,0,0.5,15,5,80);
  TH2D* deltaRvpt_p = new TH2D("deltaRvpt_p","",10,0,0.5,15,5,80);
  TH2D* deltaRvpt_other = new TH2D("deltaRvpt_other","",10,0,0.5,15,5,80);
  TH2D* zvpt_pi = new TH2D("zvpt_pi","",10,0,1,15,5,80);
  TH2D* zvpt_K = new TH2D("zvpt_K","",10,0,1,15,5,80);
  TH2D* zvpt_p = new TH2D("zvpt_p","",10,0,1,15,5,80);
  TH2D* zvpt_other = new TH2D("zvpt_other","",10,0,1,15,5,80);
  
  //need 4 1D hists of jet pT weighted by the JEF. How is JEF calculated? Add up all the pTs of the constituents and divide by the jet pT. Then weight that number with the MC weight. Keep a running tally of these numbers. Add the next jet's weighted number to this first number, and also to the running total. At the end, divide each of the four weighted numbers by the weighted total to get the JEF for each of the species. Stack them here. Or create the stack in another plotting macro if ownership is weird w/r/t prettification later.
  //do the same thing at the constituent level with z and delta_R?
  
  //3D correlations
  TH3D* mvmgvpt = new TH3D("mvmgvpt",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,ptbins,ptlow,pthigh);
  TH3D* mvzgvpt = new TH3D("mvzgvpt",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* mvrgvpt = new TH3D("mvrgvpt",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* mgvzgvpt = new TH3D("mgvzgvpt",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* mgvrgvpt = new TH3D("mgvrgvpt",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* zgvrgvpt = new TH3D("zgvrgvpt",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,ptbins,ptlow,pthigh);

  //~~~PARTON LEVEL~~~
  //2D correlations with pT
  TH2D* PLmvpt = new TH2D("PLmvpt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);
  TH2D* PLmgvpt = new TH2D("PLmgvpt",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);
  TH2D* PLetavpt = new TH2D("PLetavpt",";#eta;p_{T} [GeV/c]",50,-1,1,ptbins,ptlow,pthigh);
  TH2D* PLphivpt = new TH2D("PLphivpt",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,ptbins,ptlow,pthigh);
  TH2D* PLzgvpt = new TH2D("PLzgvpt",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,ptbins,ptlow,pthigh);
  TH2D* PLrgvpt = new TH2D("PLrgvpt",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,ptbins,ptlow,pthigh);
  
  //3D correlations
  TH3D* PLmvmgvpt = new TH3D("PLmvmgvpt",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,ptbins,ptlow,pthigh);
  TH3D* PLmvzgvpt = new TH3D("PLmvzgvpt",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* PLmvrgvpt = new TH3D("PLmvrgvpt",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* PLmgvzgvpt = new TH3D("PLmgvzgvpt",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* PLmgvrgvpt = new TH3D("PLmgvrgvpt",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* PLzgvrgvpt = new TH3D("PLzgvrgvpt",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,ptbins,ptlow,pthigh);
  
  //~~~MATCHED PARTON- TO HADRON-LEVEL~~~//
  TH2D* PLmvHLpt = new TH2D("PLmvHLpt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);
  TH2D* PLmgvHLpt = new TH2D("PLmgvHLpt",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);


  //putting them in a vector to more easily shuttle them back and forth in the function. Drawback: have to know their order.
  vector<TH1D*> h1Ds = {ptvJEF_pi, ptvJEF_K, ptvJEF_p, ptvJEF_other};
  vector<TH2D*> h2Ds = {mvpt,mvpt_q,mvpt_g,mvpt_n,etavpt,phivpt,mgvpt,zgvpt,rgvpt,deltaRvpt_pi,deltaRvpt_K,deltaRvpt_p,deltaRvpt_other,zvpt_pi,zvpt_K,zvpt_p,zvpt_other};
  vector<TH3D*> h3Ds = {mvmgvpt,mvzgvpt,mvrgvpt,mgvzgvpt,mgvrgvpt,zgvrgvpt}; 
  vector<TH1D*> PLh1Ds = {dummy1D,dummy1D,dummy1D,dummy1D};
  vector<TH2D*> PLh2Ds = {PLmvpt,dummy,dummy,dummy,PLetavpt,PLphivpt,PLmgvpt,PLzgvpt,PLrgvpt,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy};
  vector<TH3D*> PLh3Ds = {PLmvmgvpt,PLmvzgvpt,PLmvrgvpt,PLmgvzgvpt,PLmgvrgvpt,PLzgvrgvpt};
  vector<TH2D*> matchh2Ds = {PLmvHLpt,PLmgvHLpt};
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  cout << "HISTS DONE!" << endl;
  
  //inelegant way to reuse the same function twice: hardcode the branch names
  //order of appearance: pt, eta, phi, m, mg, zg, rg
  vector<string> HLbranches = {"jetpT","jeteta","jetphi","jetM","sdjetM","zg","rg"};
  vector<string> PLbranches = {"PLpt","PLeta","PLphi","PLm","PLmg","PLzg","PLrg"};
  
  //  cout << "VECTORS OF (STRING) BRANCHES DONE!" << endl;
  
  //calling analysis function(s)! "event" here is the internal name of the tree in "fin"  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TreetoHist (fin, "ResultTree", HLbranches, h1Ds, h2Ds, h3Ds, 1); //1 = hadronic
  cout << "HeLP3!" << endl;
  TreetoHist (fin, "PartonTree", PLbranches, PLh1Ds, PLh2Ds, PLh3Ds, 0); //0 = partonic
  //MatchedTreetoHist (fin, "MatchTree", matchh2Ds);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  cout << "HeLP4!" << endl;
  
  //outro
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //creating output file in which to deposit histograms
  //argv[1] should be the desired location of the output file, argv[2] should be the desired name
  TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
  cout << "DEBUG: output file name is " << fout->GetName() << endl;
  
  //writing hists to file
  for (unsigned i = 0; i < h1Ds.size(); ++ i) {
    h1Ds[i]->Write();
    PLh1Ds[i]->Write();
  }
  for (unsigned i = 0; i < h2Ds.size(); ++ i) {
    h2Ds[i]->Write();
    PLh2Ds[i]->Write();
  }
  for (unsigned i = 0; i < h3Ds.size(); ++ i) {
    h3Ds[i]->Write();
    PLh2Ds[i]->Write();
  }
  for (unsigned i = 0; i < matchh2Ds.size(); ++ i) {
    //matchh2Ds[i]->Write();
  }
 
  cout << "Wrote to " << fout->GetName() << endl;
  
  //closing file
  fout->Close();
  cout << "Closed " << fout->GetName() << endl;
  fin->Close();
  cout << "Closed " << fin->GetName() << endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  return 0;
}
