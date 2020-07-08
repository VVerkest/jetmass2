{

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };
  
  TFile *inFile = new TFile("out/sim/Cleanpp12Pico_full_matched_R04.root","READ");

  vector<double> *p_Pt, *g_Pt;

  TTree *eTree = (TTree*)inFile->Get("event");
  eTree->SetBranchAddress("p_Pt", &p_Pt);
  eTree->SetBranchAddress("g_Pt", &g_Pt);

  TH2D *hPtResponse = new TH2D("hPtResponse",";part. jet p_{T} (GeV);det. jet p_{T} (GeV)", 55,4.5,59.5, 55,4.5,59.5);
  
  for (int i=0; i<eTree->GetEntries(); ++i) {
    eTree->GetEvent(i);
    if ( p_Pt->size()==0 || g_Pt->size()==0) { continue; }
    int maxval = min( p_Pt->size(), g_Pt->size() );
    for (int j=0; j<maxval; ++j) {
      double ppt = (double) p_Pt->at(j);
      double gpt = (double) g_Pt->at(j);
      hPtResponse->Fill( ppt , gpt );
    }
  }

  hPtResponse->Scale(1./hPtResponse->GetEntries());

  TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );              // CANVAS 0
  can->SetLogz();
  hPtResponse->Draw("COLZ");
  can->SaveAs("pTresponse.pdf","PDF");

  hPtResponse->Scale(hPtResponse->GetEntries());
  
  TH1D *hDet[nPtBins]; TH1D *hPart[nPtBins];
  
  for (int p=0; p<nPtBins; ++p) {
    hPtResponse->GetXaxis()->SetRangeUser(ptLo[p],ptHi[p]);
    TString name = "hPtResponse"; name += ptBinName[p];
    hDet[p] = (TH1D*) hPtResponse->ProjectionY(name);
    hDet[p]->DrawNormalized();
    name += ".pdf";
    can->SaveAs( name, "PDF" );
    hDet[p]->Scale(1./hDet[p]->GetEntries());
  }

  TFile *UEfile = new TFile("../pAu_analysis/out/UE/pAuHTjetUE_UE3Dplot.root","READ");

  TH3D *hChgUE3D = (TH3D*)UEfile->Get("hChgUE");  // leadPt:UEpt:UEeta
  TH1D *hleadPt = (TH1D*)UEfile->Get("hLeadPt");  // leadPt

  TH2D *hChgUE2D = (TH2D*)hChgUE3D->Project3D("YX");
  //hChgUE2D->Scale(1./hleadPt->GetEntries());
  hChgUE2D->GetXaxis()->SetRangeUser(0.0,60.0);
  hChgUE2D->GetYaxis()->SetRangeUser(0.0,15.0);
  hChgUE2D->GetYaxis()->SetTitleOffset(1.25);
  hChgUE2D->Draw("COLZ");
  can->SaveAs("hChgUE2D.pdf","PDF");
  hChgUE2D->GetXaxis()->SetRangeUser(1,-1);
  hChgUE2D->GetYaxis()->SetRangeUser(1,-1);
  //hChgUE2D->Scale(hleadPt->GetEntries());
  
  TH1D *hUEpt[55];

  double /*nCh[55], ptMean[55],*/ detPt[nPtBins], detChPt[nPtBins], detNch[nPtBins];

  TH1D *nCh[nPtBins];
  TH1D *ptMean[nPtBins];

  TH1D *nCh_wt[nPtBins];
  TH1D *ptMean_wt[nPtBins];

  for (int p=0; p<nPtBins; ++p) {
    TString name = "nCh"; name += ptBinName[p];
    nCh[p] = new TH1D(name,";det-level #LT N_{ch}#GT", 60,0.0,6.0);
    name = "ptMean"; name += ptBinName[p];
    ptMean[p] = new TH1D(name,";det-level #LT p_{T}^{ch}#GT", 30,0.63,0.72);
    name = "nCh_wt"; name += ptBinName[p];
    nCh_wt[p] = new TH1D(name,";part-level #LT N_{ch}#GT", 60,0.0,6.0);
    name = "ptMean_wt"; name += ptBinName[p];
    ptMean_wt[p] = new TH1D(name,";part-level #LT p_{T}^{ch}#GT", 30,0.63,0.72);
  }

  can->SetLogy();
  for (int i=0; i<55; ++i) {  // det-level fractional pT contribution to part-level pT --> fc(pT_det)
    int ptVal = i + 5;
    int binno = i + 1;
    TString name = "hUEpt_"; name += ptVal; name += "GeV";
    hUEpt[i] = (TH1D*)hChgUE2D->ProjectionY(name,binno,binno);
    if ( (int)hleadPt->GetBinCenter(binno) != ptVal ) { cerr<<"failed pT matching"<<endl; }
    int nJets = hleadPt->GetBinContent( binno );
    if (nJets>0){
      hUEpt[i]->Scale(1./nJets);
      hUEpt[i]->Draw();
      name += ".pdf";
      can->SaveAs(name,"PDF");
      
      for (int p=0; p<nPtBins; ++p) {
	if (ptVal > ptLo[p] && ptVal < ptHi[p] ) {
	  ptMean[p]->Fill(hUEpt[i]->GetMean(1));
	  nCh[p]->Fill(hUEpt[i]->Integral());
	}
	  ptMean_wt[p]->Fill(hUEpt[i]->GetMean(1),hDet[p]->GetBinContent(binno));
	  nCh_wt[p]->Fill(hUEpt[i]->Integral(),hDet[p]->GetBinContent(binno));
      }
	// ptMean_p->SetBinContent(i,hUEpt[i]->GetMean(1));
	// ptMean_p->SetBinError(i,hUEpt[i]->GetMeanError(1));
	// nCh_p->SetBinContent(i,hUEpt[i]->Integral());
	// ptMean[i] = hUEpt[i]->GetMean(1);
	// nCh[i] = hUEpt[i]->Integral();
    }
    //else { ptMean[i] = 0.0;  nCh[i] = 0.0; }
  }

  can->SetLogy();
  for (int p=0; p<nPtBins; ++p) {
    ptMean_wt[p]->Scale(1./ptMean_wt[p]->Integral("width"));
    ptMean_wt[p]->Draw();
    TString name = "ptMean_wt"; name += ptBinName[p]; name += ".pdf";
    can->SaveAs(name,"PDF");

    nCh_wt[p]->Scale(1./nCh_wt[p]->Integral("width"));
    nCh_wt[p]->Draw();
    name = "nCh_wt"; name += ptBinName[p]; name += ".pdf";
    can->SaveAs(name,"PDF");

    ptMean[p]->Scale(1./ptMean[p]->Integral("width"));
    ptMean[p]->Draw();
    name = "ptMean"; name += ptBinName[p]; name += ".pdf";
    can->SaveAs(name,"PDF");

    nCh[p]->Scale(1./nCh[p]->Integral("width"));
    nCh[p]->Draw();
    name = "nCh"; name += ptBinName[p]; name += ".pdf";
    can->SaveAs(name,"PDF");
  }
  can->SetLogy(0);
  
  // for (int p=0; p<nPtBins; ++p) {

  //   for (int i=0; i<55; ++i) {
      
  //     int ptVal = i + 5;
  //     int binno = i + 1;
  //     double content, error;
	
  //     if ( (int)hleadPt->GetBinCenter(binno) != ptVal ) { cerr<<"failed pT matching"<<endl; continue; }
  //     if ( isnan(hDet[p]->GetBinContent(binno)) ) { continue; }

  //     if ( ptVal >= ptLo[p] && ptVal < ptHi[p] ) {

  //     	content = hUEpt[i]->GetMean(1)*hDet[p]->GetBinContent(binno);
  //     	ptMean_p[p]->SetBinContent(i,content);
  //     	error = hUEpt[i]->GetMean(1)*(hDet[p]->GetBinError(binno)*hDet[p]->GetBinContent(binno));
  //     	ptMean_p[p]->SetBinError(i,error);

  //     	content = hUEpt[i]->Integral()*hDet[p]->GetBinContent(binno);
  //     	nCh_p[p]->SetBinContent(i,content);
  // 	// 	// error = hUEpt[i]->Integral()*(hDet[p]->GetBinError(binno)*hDet[p]->GetBinContent(binno));
  // 	// 	// nCh_p[p]->SetBinError(i,error);

	
  //     }
  //   }
    
  // }
  
  // TH1D *hDetPt[nPtBins], *hDetNch[nPtBins];
  // double value, error;
  // TString val;
  
  // for (int p=0; p<nPtBins; ++p) {

  //   TString name = "hPartPt"; name += ptBinName[p];
  //   hDetPt[p] = new TH1D(name,";det. lead p_{T} (GeV)",55,4.5,59.5);
  //   name = "hPartNch"; name += ptBinName[p];
  //   hDetNch[p] = new TH1D(name,";det. lead p_{T} (GeV)",55,4.5,59.5);

  //   detPt[p] = 0;  detChPt[p] = 0;  detNch[p] = 0;
    
  //     for (int i=0; i<55; ++i) {

  // 	int ptVal = i + 5;
  // 	int binno = i + 1;

  // 	if ( (int)hleadPt->GetBinCenter(binno) != ptVal ) { cerr<<"failed pT matching"<<endl; }

  // 	if ( !isnan(hDet[p]->GetBinContent(binno)) ) {
  // 	  detPt[p] += ptVal * hDet[p]->GetBinContent(binno);

  // 	  if ( ptVal >= ptLo[p] && ptVal < ptHi[p] ) {
  // 	    detChPt[p] += ptMean[i];
  // 	    detNch[p] += nCh[i];
  // 	  }

  // 	  if (hDet[p]->GetBinContent(binno)>0.0) {
  // 	    value = ptMean[i] * hDet[p]->GetBinContent(binno);
  // 	    hDetPt[p]->SetBinContent( binno, value );
  // 	    // error = ptMean[i] * ( hDet[p]->GetBinError(binno)/hDet[p]->GetBinContent(binno) );
  // 	    // hDetPt[p]->SetBinError( binno, error );
  // 	  }
  // 	  else { hDetPt[p]->SetBinContent( binno, 0.0 ); }
	  
  // 	  if ( hDet[p]->GetBinContent(binno)>0.0) {
  // 	    value = nCh[i] * hDet[p]->GetBinContent(binno);
  // 	    hDetNch[p]->SetBinContent( binno, value );
  // 	    // error = nCh[i] * ( hDet[p]->GetBinError(binno)/hDet[p]->GetBinContent(binno) );
  // 	    // hDetNch[p]->SetBinError( binno, error );
  // 	  }
  // 	  else { hDetNch[p]->SetBinContent( binno, 0.0 ); }
  // 	}
  //     }

  //     detChPt[p] /= (ptHi[p]-ptLo[p]);
  //     detNch[p] /= (ptHi[p]-ptLo[p]);

  //     value = detChPt[p];
  //     val = "det. #LT p_{T}^{ch}#GT = "; val += value; val=val(0,31);
  //     TLatex *tex2 = new TLatex(0.6,0.5,val);
  //     tex2->SetTextFont(63);      tex2->SetTextSize(16);      tex2->SetTextColor(kBlack);      tex2->SetLineWidth(1);      tex2->SetNDC();
  //     hDetPt[p]->Draw();
  //     tex2->Draw();
      
  //     value = hDetPt[p]->Integral();
  //     val = "part. #LT p_{T}^{ch}#GT = "; val += value; val=val(0,32);
  //     tex2->DrawLatex(0.6,0.4,val);
      
  //     name = hDetPt[p]->GetName(); name += ".pdf";
  //     can->SaveAs(name,"PDF");

      
  //     value = detNch[p];
  //     val = "det. #LT N_{ch}#GT = "; /*"part. #LT #frac{dN_{ch}}{d#eta d#phi}#GT = ";*/ val += value; val=val(0,26);
  //     TLatex *tex1 = new TLatex(0.6,0.5,val);
  //     tex1->SetTextFont(63);      tex1->SetTextSize(16);      tex1->SetTextColor(kBlack);      tex1->SetLineWidth(1);      tex1->SetNDC();
  //     hDetNch[p]->Draw();
  //     tex1->Draw();
  //     value = hDetNch[p]->Integral();
  //     val = "part. #LT N_{ch}#GT = "; /*"part. #LT #frac{dN_{ch}}{d#eta d#phi}#GT = ";*/ val += value; val=val(0,27);
  //     tex1->DrawLatex(0.6,0.4,val);

  //     name = hDetNch[p]->GetName(); name += ".pdf";
  //     can->SaveAs(name,"PDF");
  // }
  
}
