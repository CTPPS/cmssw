void PlotTimingEfficiencies(TString infile = "timingHistograms.root")
{
  TFile *f = TFile::Open(infile);

  TH2F *den45 = (TH2F *)f->Get("diamondHistograms/deffden45");
  TH2F *den56 = (TH2F *)f->Get("diamondHistograms/deffden56");

  TH2F *num45 = (TH2F *)f->Get("diamondHistograms/deffnum45");
  TH2F *num56 = (TH2F *)f->Get("diamondHistograms/deffnum56");

  TH2F *boxnum45 = (TH2F *)f->Get("diamondHistograms/dboxeffnum45");
  TH2F *boxnum56 = (TH2F *)f->Get("diamondHistograms/dboxeffnum56");

  TH2F *pl0st145 = (TH2F *)f->Get("diamondHistograms/deffnum45plane0");
  TH2F *pl1st145 = (TH2F *)f->Get("diamondHistograms/deffnum45plane1");
  TH2F *pl2st145 = (TH2F *)f->Get("diamondHistograms/deffnum45plane2");
  TH2F *pl3st145 = (TH2F *)f->Get("diamondHistograms/deffnum45plane3");

  TH2F *pl0st156 = (TH2F *)f->Get("diamondHistograms/deffnum56plane0");
  TH2F *pl1st156 = (TH2F *)f->Get("diamondHistograms/deffnum56plane1");
  TH2F *pl2st156 = (TH2F *)f->Get("diamondHistograms/deffnum56plane2");
  TH2F *pl3st156 = (TH2F *)f->Get("diamondHistograms/deffnum56plane3");

  TH2F *pl0st245box = (TH2F *)f->Get("diamondHistograms/deffnum45boxplane0");
  TH2F *pl1st245box = (TH2F *)f->Get("diamondHistograms/deffnum45boxplane1");
  TH2F *pl2st245box = (TH2F *)f->Get("diamondHistograms/deffnum45boxplane2");
  TH2F *pl3st245box = (TH2F *)f->Get("diamondHistograms/deffnum45boxplane3");

  TH2F *pl0st256box = (TH2F *)f->Get("diamondHistograms/deffnum56boxplane0");
  TH2F *pl1st256box = (TH2F *)f->Get("diamondHistograms/deffnum56boxplane1");
  TH2F *pl2st256box = (TH2F *)f->Get("diamondHistograms/deffnum56boxplane2");
  TH2F *pl3st256box = (TH2F *)f->Get("diamondHistograms/deffnum56boxplane3");


  pl0st145->Divide(den45);
  pl1st145->Divide(den45);
  pl2st145->Divide(den45);
  pl3st145->Divide(den45);
  pl0st156->Divide(den56);
  pl1st156->Divide(den56);
  pl2st156->Divide(den56);
  pl3st156->Divide(den56);

  pl0st245box->Divide(den45);
  pl1st245box->Divide(den45);
  pl2st245box->Divide(den45);
  pl3st245box->Divide(den45);
  pl0st256box->Divide(den56);
  pl1st256box->Divide(den56);
  pl2st256box->Divide(den56);
  pl3st256box->Divide(den56);

  num45->Divide(den45);
  num56->Divide(den56);

  boxnum45->Divide(den45);
  boxnum56->Divide(den56);

  TCanvas *c1 = new TCanvas("c1","c1",800,400);
  c1->Divide(4,2);
  c1->cd(1); pl0st145->SetMaximum(1); pl0st145->Draw("colz"); pl0st145->SetTitle("45 cyl plane0"); pl0st145->SetStats(0);
  c1->cd(2); pl1st145->SetMaximum(1); pl1st145->Draw("colz"); pl1st145->SetTitle("45 cyl plane1"); pl1st145->SetStats(0);
  c1->cd(5); pl2st145->SetMaximum(1); pl2st145->Draw("colz"); pl2st145->SetTitle("45 cyl plane2"); pl2st145->SetStats(0);
  c1->cd(6); pl3st145->SetMaximum(1); pl3st145->Draw("colz"); pl3st145->SetTitle("45 cyl plane3"); pl3st145->SetStats(0);

  c1->cd(3); pl0st156->SetMaximum(1); pl0st156->Draw("colz"); pl0st156->SetTitle("56 cyl plane0"); pl0st156->SetStats(0);
  c1->cd(4); pl1st156->SetMaximum(1); pl1st156->Draw("colz"); pl1st156->SetTitle("56 cyl plane1"); pl1st156->SetStats(0);
  c1->cd(7); pl2st156->SetMaximum(1); pl2st156->Draw("colz"); pl3st156->SetTitle("56 cyl plane2"); pl2st156->SetStats(0);
  c1->cd(8); pl3st156->SetMaximum(1); pl3st156->Draw("colz"); pl3st156->SetTitle("56 cyl plane3"); pl3st156->SetStats(0);

  TCanvas *c2 = new TCanvas("c2","c2",800,400);
  c2->Divide(4,2);
  c2->cd(1); pl0st245box->SetMaximum(1); pl0st245box->Draw("colz"); pl0st245box->SetTitle("45 box plane0"); pl0st245box->SetStats(0);
  c2->cd(2); pl1st245box->SetMaximum(1); pl1st245box->Draw("colz"); pl1st245box->SetTitle("45 box plane1"); pl1st245box->SetStats(0);
  c2->cd(5); pl2st245box->SetMaximum(1); pl2st245box->Draw("colz"); pl2st245box->SetTitle("45 box plane2"); pl2st245box->SetStats(0);
  c2->cd(6); pl3st245box->SetMaximum(1); pl3st245box->Draw("colz"); pl3st245box->SetTitle("45 box plane3"); pl3st245box->SetStats(0);

  c2->cd(3); pl0st256box->SetMaximum(1); pl0st256box->Draw("colz"); pl0st256box->SetTitle("56 box plane0"); pl0st256box->SetStats(0);
  c2->cd(4); pl1st256box->SetMaximum(1); pl1st256box->Draw("colz"); pl1st256box->SetTitle("56 box plane1"); pl1st256box->SetStats(0);
  c2->cd(7); pl2st256box->SetMaximum(1); pl2st256box->Draw("colz"); pl2st256box->SetTitle("56 box plane2"); pl2st256box->SetStats(0);
  c2->cd(8); pl3st256box->SetMaximum(1); pl3st256box->Draw("colz"); pl3st256box->SetTitle("56 box plane3"); pl3st256box->SetStats(0);


  TCanvas *c3 = new TCanvas("c3","c3");
  c3->Divide(2,1);
  c3->cd(1); num45->Draw("colz"); num45->SetTitle("45 time-tracks"); num45->SetStats(0); num45->SetMaximum(1);
  c3->cd(2); num56->Draw("colz"); num56->SetTitle("56 time-tracks"); num56->SetStats(0); num56->SetMaximum(1);

  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(2,1);
  c4->cd(1); boxnum45->Draw("colz"); boxnum45->SetTitle("45 time-tracks"); boxnum45->SetStats(0); boxnum45->SetMaximum(1);
  c4->cd(2); boxnum56->Draw("colz"); boxnum56->SetTitle("56 time-tracks"); boxnum56->SetStats(0); boxnum56->SetMaximum(1);

  TObjArray *tx = infile.Tokenize(".");
  TString outplot = ((TObjString *)(tx->At(0)))->String();
  c1->SaveAs("CylindricalPlanes_" + outplot + ".pdf");
  c2->SaveAs("BoxPlanes_" + outplot + ".pdf");
  c3->SaveAs("CylindricalTimeTracks_" + outplot + ".pdf");
  c4->SaveAs("BoxTimeTracks_" + outplot + ".pdf");

}


void PlotToTs(TString infile = "timingHistograms.root")
{
  gStyle->SetPalette(107);

  TFile *f = TFile::Open(infile);

  TH2F *pl0st145 = (TH2F *)f->Get("diamondHistograms/deffnum45plane0");
  TH2F *pl1st145 = (TH2F *)f->Get("diamondHistograms/deffnum45plane1");
  TH2F *pl2st145 = (TH2F *)f->Get("diamondHistograms/deffnum45plane2");
  TH2F *pl3st145 = (TH2F *)f->Get("diamondHistograms/deffnum45plane3");

  TH2F *pl0st156 = (TH2F *)f->Get("diamondHistograms/deffnum56plane0");
  TH2F *pl1st156 = (TH2F *)f->Get("diamondHistograms/deffnum56plane1");
  TH2F *pl2st156 = (TH2F *)f->Get("diamondHistograms/deffnum56plane2");
  TH2F *pl3st156 = (TH2F *)f->Get("diamondHistograms/deffnum56plane3");

  TH2F *pl0st245box = (TH2F *)f->Get("diamondHistograms/deffnum45boxplane0");
  TH2F *pl1st245box = (TH2F *)f->Get("diamondHistograms/deffnum45boxplane1");
  TH2F *pl2st245box = (TH2F *)f->Get("diamondHistograms/deffnum45boxplane2");
  TH2F *pl3st245box = (TH2F *)f->Get("diamondHistograms/deffnum45boxplane3");

  TH2F *pl0st256box = (TH2F *)f->Get("diamondHistograms/deffnum56boxplane0");
  TH2F *pl1st256box = (TH2F *)f->Get("diamondHistograms/deffnum56boxplane1");
  TH2F *pl2st256box = (TH2F *)f->Get("diamondHistograms/deffnum56boxplane2");
  TH2F *pl3st256box = (TH2F *)f->Get("diamondHistograms/deffnum56boxplane3");

  // 
  TH2F *totpl0st145 = (TH2F *)f->Get("diamondHistograms/tot45plane0");
  TH2F *totpl1st145 = (TH2F *)f->Get("diamondHistograms/tot45plane1");
  TH2F *totpl2st145 = (TH2F *)f->Get("diamondHistograms/tot45plane2");
  TH2F *totpl3st145 = (TH2F *)f->Get("diamondHistograms/tot45plane3");

  TH2F *totpl0st156 = (TH2F *)f->Get("diamondHistograms/tot56plane0");
  TH2F *totpl1st156 = (TH2F *)f->Get("diamondHistograms/tot56plane1");
  TH2F *totpl2st156 = (TH2F *)f->Get("diamondHistograms/tot56plane2");
  TH2F *totpl3st156 = (TH2F *)f->Get("diamondHistograms/tot56plane3");

  TH2F *totpl0st245box = (TH2F *)f->Get("diamondHistograms/tot45boxplane0");
  TH2F *totpl1st245box = (TH2F *)f->Get("diamondHistograms/tot45boxplane1");
  TH2F *totpl2st245box = (TH2F *)f->Get("diamondHistograms/tot45boxplane2");
  TH2F *totpl3st245box = (TH2F *)f->Get("diamondHistograms/tot45boxplane3");

  TH2F *totpl0st256box = (TH2F *)f->Get("diamondHistograms/tot56boxplane0");
  TH2F *totpl1st256box = (TH2F *)f->Get("diamondHistograms/tot56boxplane1");
  TH2F *totpl2st256box = (TH2F *)f->Get("diamondHistograms/tot56boxplane2");
  TH2F *totpl3st256box = (TH2F *)f->Get("diamondHistograms/tot56boxplane3");

  totpl0st145->Divide(pl0st145);
  totpl1st145->Divide(pl1st145);
  totpl2st145->Divide(pl2st145);
  totpl3st145->Divide(pl3st145);

  totpl0st156->Divide(pl0st156);
  totpl1st156->Divide(pl1st156);
  totpl2st156->Divide(pl2st156);
  totpl3st156->Divide(pl3st156);

  totpl0st245box->Divide(pl0st245box);
  totpl1st245box->Divide(pl1st245box);
  totpl2st245box->Divide(pl2st245box);

  totpl0st256box->Divide(pl0st256box);
  totpl1st256box->Divide(pl1st256box);
  totpl2st256box->Divide(pl2st256box);

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(4,2);
  c1->cd(1); totpl0st145->SetMaximum(15); totpl0st145->SetMinimum(7); totpl0st145->Draw("colz"); totpl0st145->SetTitle("45 cyl plane0"); totpl0st145->SetStats(0);
  c1->cd(2); totpl1st145->SetMaximum(15); totpl1st145->SetMinimum(7); totpl1st145->Draw("colz"); totpl1st145->SetTitle("45 cyl plane1"); totpl1st145->SetStats(0);
  c1->cd(5); totpl2st145->SetMaximum(15); totpl2st145->SetMinimum(7); totpl2st145->Draw("colz"); totpl2st145->SetTitle("45 cyl plane2"); totpl2st145->SetStats(0);
  c1->cd(6); totpl3st145->SetMaximum(15); totpl3st145->SetMinimum(7); totpl3st145->Draw("colz"); totpl3st145->SetTitle("45 cyl plane3"); totpl3st145->SetStats(0);

  c1->cd(3); totpl0st156->SetMaximum(15); totpl0st156->SetMinimum(7); totpl0st156->Draw("colz"); totpl0st156->SetTitle("56 cyl plane0"); totpl0st156->SetStats(0);
  c1->cd(4); totpl1st156->SetMaximum(15); totpl1st156->SetMinimum(7); totpl1st156->Draw("colz"); totpl1st156->SetTitle("56 cyl plane1"); totpl1st156->SetStats(0);
  c1->cd(7); totpl2st156->SetMaximum(15); totpl2st156->SetMinimum(7); totpl2st156->Draw("colz"); totpl3st156->SetTitle("56 cyl plane2"); totpl2st156->SetStats(0);
  c1->cd(8); totpl3st156->SetMaximum(15); totpl3st156->SetMinimum(7); totpl3st156->Draw("colz"); totpl3st156->SetTitle("56 cyl plane3"); totpl3st156->SetStats(0);

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(4,2);
  c2->cd(1); totpl0st245box->SetMaximum(15); totpl0st245box->SetMinimum(7); totpl0st245box->Draw("colz"); totpl0st245box->SetTitle("45 box plane0"); totpl0st245box->SetStats(0);
  c2->cd(2); totpl1st245box->SetMaximum(15); totpl1st245box->SetMinimum(7); totpl1st245box->Draw("colz"); totpl1st245box->SetTitle("45 box plane1"); totpl1st245box->SetStats(0);
  c2->cd(5); totpl2st245box->SetMaximum(15); totpl2st245box->SetMinimum(7); totpl2st245box->Draw("colz"); totpl2st245box->SetTitle("45 box plane2"); totpl2st245box->SetStats(0);
  c2->cd(6); totpl3st245box->SetMaximum(15); totpl3st245box->SetMinimum(7); totpl3st245box->Draw("colz"); totpl3st245box->SetTitle("45 box plane3"); totpl3st245box->SetStats(0);

  c2->cd(3); totpl0st256box->SetMaximum(15); totpl0st256box->SetMinimum(7); totpl0st256box->Draw("colz"); totpl0st256box->SetTitle("56 box plane0"); totpl0st256box->SetStats(0);
  c2->cd(4); totpl1st256box->SetMaximum(15); totpl1st256box->SetMinimum(7); totpl1st256box->Draw("colz"); totpl1st256box->SetTitle("56 box plane1"); totpl1st256box->SetStats(0);
  c2->cd(7); totpl2st256box->SetMaximum(15); totpl2st256box->SetMinimum(7); totpl2st256box->Draw("colz"); totpl3st256box->SetTitle("56 box plane2"); totpl2st256box->SetStats(0);
  c2->cd(8); totpl3st256box->SetMaximum(15); totpl3st256box->SetMinimum(7); totpl3st256box->Draw("colz"); totpl3st256box->SetTitle("56 box plane3"); totpl3st256box->SetStats(0);

  TObjArray *tx = infile.Tokenize(".");
  TString outplot = ((TObjString *)(tx->At(0)))->String();
  c1->SaveAs("ToT_CylindricalPlanes_" + outplot + ".pdf");
  c2->SaveAs("ToT_BoxPlanes_" + outplot + ".pdf");
}

void PlotRadiographies()
{
  TFile *f = TFile::Open("timingHistograms_Run378867_LS1to1000.root");

  TH2F *rad45 = (TH2F *)f->Get("diamondHistograms/drad45");
  TH2F *rad56 = (TH2F *)f->Get("diamondHistograms/drad56");

  TH2F *rad45box = (TH2F *)f->Get("diamondHistograms/dboxrad45");
  TH2F *rad56box = (TH2F *)f->Get("diamondHistograms/dboxrad56");


  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(2,1);
  c1->cd(1); rad45->Draw("colz"); rad45->SetStats(0); rad45->SetTitle("45 cylindrical"); 
  c1->cd(2); rad56->Draw("colz"); rad56->SetStats(0); rad56->SetTitle("56 cylindrical");

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(2,1);
  c2->cd(1); rad45box->Draw("colz"); rad45box->SetStats(0); rad45box->SetTitle("45 box");
  c2->cd(2); rad56box->Draw("colz"); rad56box->SetStats(0); rad56box->SetTitle("56 box");

}
