// Parametrise the Pb-Pb 5.02 TeV 30-50% model prediction of Pz,s2 vs pT (black curve in PzsVsPt_PbPbComp),
// rescale it so that its spectrum-weighted average matches the OO 0-50% pT-integrated Pzs2,
// and overlay the rescaled curve on the PzsVsPt_PbPbComp plot.
//
// Scale factor: k = Pzs2(OO, 0-50%) / <f>_w, with <f>_w = sum_i w_i f(pT_i) / sum_i w_i over pT > ptMin.
// Weights: Pb-Pb 30-50% Lambda + AntiLambda corrected spectra (isOOWeights = 0, default), as done for the original
// normalisation, or the OO 0-50% raw Lambda yield vs pT (isOOWeights = 1), which is how the OO pT-integrated value is built
// when no efficiency weights are applied.
//
// The OO 0-50% pT-integrated value and its stat. error are read from fileNameOOInt (output of FitV2orPol_CB.C for Cent0-50):
// histoPzs2PtInt if isFromFit = 1, histoPzs2PtIntNoFit / histoPurityPtInt if isFromFit = 0 (as in PzsVsCentrality.C).
// The OO points are read from fileNameOOPt (output of PzsVsPt.C): histoPzs2 if isFromFit = 1, histoPzs2NoFit otherwise.
//
// Usage: root -l -b -q ScaledModelVsPt.C
//        root -l -b -q 'ScaledModelVsPt.C(1)'
//        root -l -b -q 'ScaledModelVsPt.C(0, 0.5, 0, "OutputAnalysis/<Cent0-50 fit file>.root", "Pzs2VsPt/<pT file, NoFit>.root", "_NoFit")'

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TDirectoryFile.h"
#include <fstream>
#include <iostream>

#include "CommonVarPub.h" // pTPalermo, Pzs2VsPtPalermo

TGraph *ReadModel(TString path)
{
  std::ifstream in(path.Data());
  if (!in)
  {
    std::cout << "Cannot open " << path << std::endl;
    return nullptr;
  }
  TGraph *g = new TGraph();
  double x, y;
  while (in >> x >> y)
    g->AddPoint(x, y / 100.); // file is in %
  return g;
}

TH1 *GetHisto(TFile *f, TString name)
{
  TH1 *h = f ? (TH1 *)f->Get(name) : nullptr;
  if (!h)
    std::cout << "Missing histogram " << name << (f ? TString(" in ") + f->GetName() : TString("")) << std::endl;
  return h;
}

void ScaledModelVsPt(Bool_t isOOWeights = 0,
                     Double_t ptMin = 0.5,
                     Bool_t isFromFit = 1,
                     TString fileNameOOInt = "OutputAnalysis/FitPzs2_LHC25_OO_pass2_Train742311_Lambda_DSCB_BkgCheb_Cent0-50_CentWeighted_EffW_Eta08_TightMassCut2.1_ReducedPtBins_AcceptanceInMacro_ResoOnTheFly_isBkgPol0_Rebin.root",
                     TString fileNameOOPt = "Pzs2VsPt/Pzs2_LHC25_OO_pass2_Train742311_Lambda_DSCB_BkgCheb_Pzs2_CentWeighted_Eta08_TightMassCut2.1_ReducedPtBins_AcceptanceInMacro_ResoOnTheFly_isBkgPol0_EffW_Rebin.root",
                     TString outSuffix = "", // appended to the output file names
                     TString dir = "/Users/mbp-cdm-01/Desktop/ResearchTeachingActivity/MyRun3Analyses/CascadesFlow/")
{
  gStyle->SetOptStat(0);

  // 0. OO 0-50% pT-integrated value
  TFile *fileOOInt = TFile::Open(dir + fileNameOOInt);
  if (!fileOOInt || fileOOInt->IsZombie())
  {
    std::cout << "Cannot open " << fileNameOOInt << std::endl;
    return;
  }
  Double_t PzsOO = 0, PzsOOErr = 0;
  if (isFromFit)
  {
    TH1 *hInt = GetHisto(fileOOInt, "histoPzs2PtInt");
    if (!hInt)
      return;
    PzsOO = hInt->GetBinContent(1);
    PzsOOErr = hInt->GetBinError(1);
  }
  else
  {
    TH1 *hInt = GetHisto(fileOOInt, "histoPzs2PtIntNoFit");
    TH1 *hPur = GetHisto(fileOOInt, "histoPurityPtInt");
    if (!hInt || !hPur)
      return;
    PzsOO = hInt->GetBinContent(1) / hPur->GetBinContent(1);
    PzsOOErr = hInt->GetBinError(1) / hPur->GetBinContent(1);
  }
  std::cout << "\nOO 0-50% pT-integrated Pzs2 (" << (isFromFit ? "fit" : "NoFit / purity") << ") from " << fileNameOOInt << std::endl;
  std::cout << Form("Pzs2 = %.6f +- %.6f (stat)", PzsOO, PzsOOErr) << std::endl;

  // 1. Parametrisation of the black curve
  TGraph *gModel502 = ReadModel(dir + "Fig_5_Prediction_1-1-2_model_PbPb_5020_30-50_LocalSpinPol_z_s2_vs_pT.dat");
  TGraph *gModel536 = ReadModel(dir + "Fig_5_Prediction_1-1-2_model_PbPb_5360_30-50_LocalSpinPol_z_s2_vs_pT.dat");
  if (!gModel502 || !gModel536)
    return;
  TGraphErrors *gFit = new TGraphErrors(gModel502->GetN());
  for (Int_t i = 0; i < gModel502->GetN(); i++)
  {
    gFit->SetPoint(i, gModel502->GetX()[i], gModel502->GetY()[i]);
    gFit->SetPointError(i, 0, 0.01 * gModel502->GetY()[i]); // equal relative weight to every point
  }
  TF1 *fModel = new TF1("fModel", "[0]*x*x/pow(x*x+[1],[2])*exp(-[3]*x)", 0, 10);
  fModel->SetParameters(1.5e-3, 1., 1.2, 0.05);
  for (Int_t it = 0; it < 3; it++)
    gFit->Fit(fModel, "QN0M", "", 0.1, 5.5);
  Double_t maxDev = 0;
  for (Int_t i = 0; i < gModel502->GetN(); i++)
    if (gModel502->GetX()[i] >= 0.2)
      maxDev = TMath::Max(maxDev, TMath::Abs(fModel->Eval(gModel502->GetX()[i]) / gModel502->GetY()[i] - 1));
  std::cout << "\nParametrisation f(pT) = A pT^2 / (pT^2 + B)^C exp(-D pT)" << std::endl;
  std::cout << Form("A = %.6g, B = %.6g, C = %.6g, D = %.6g", fModel->GetParameter(0), fModel->GetParameter(1), fModel->GetParameter(2), fModel->GetParameter(3)) << std::endl;
  std::cout << Form("max |f/model - 1| for pT >= 0.2 GeV/c: %.2f%%", 100 * maxDev) << std::endl;

  // 2. Weights
  TH1 *hW = nullptr;
  if (!isOOWeights)
  {
    TFile *fL = TFile::Open(dir + "CorrectedSpectra_Lambda_withEvtLoss_withFeeddown.root");
    TFile *fA = TFile::Open(dir + "CorrectedSpectra_Alambda_withEvtLoss_withFeeddown.root");
    TH1 *h1 = GetHisto(fL, "hLambdaCorrected1_CENT_30_40"), *h2 = GetHisto(fL, "hLambdaCorrected1_CENT_40_50");
    TH1 *h3 = GetHisto(fA, "hAlambdaCorrected1_CENT_30_40"), *h4 = GetHisto(fA, "hAlambdaCorrected1_CENT_40_50");
    if (!h1 || !h2 || !h3 || !h4)
      return;
    hW = (TH1 *)h1->Clone("hWeights");
    hW->Add(h2);
    hW->Add(h3);
    hW->Add(h4);
  }
  else
  {
    TH1 *hY = GetHisto(fileOOInt, "histoYield"); // raw yield per event per unit pT
    if (!hY)
      return;
    hW = (TH1 *)hY->Clone("hWeights");
  }
  Double_t num = 0, den = 0;
  for (Int_t b = 1; b <= hW->GetNbinsX(); b++)
  {
    Double_t lo = hW->GetBinLowEdge(b), hi = hW->GetBinLowEdge(b + 1);
    if (hi <= ptMin)
      continue;
    lo = TMath::Max(lo, ptMin);
    Double_t w = hW->GetBinContent(b) * (hi - lo);
    num += w * fModel->Integral(lo, hi) / (hi - lo); // average of f over the bin
    den += w;
  }
  Double_t avgModel = num / den;
  Double_t k = PzsOO / avgModel;
  std::cout << "\nWeights: " << (isOOWeights ? "OO 0-50% raw Lambda yield" : "Pb-Pb 30-50% Lambda+AntiLambda spectra") << Form(", pT > %.1f GeV/c", ptMin) << std::endl;
  std::cout << Form("<f>_w = %.6f  ->  scale factor k = %.6f / %.6f = %.3f (+- %.3f from the OO stat. error)", avgModel, PzsOO, avgModel, k, k * PzsOOErr / PzsOO) << std::endl;

  TF1 *fScaled = new TF1("fScaled", Form("%.10g*fModel", k), 0, 6);
  fScaled->SetLineColor(kMagenta + 1);
  fScaled->SetLineStyle(2);
  fScaled->SetLineWidth(3);
  fScaled->SetNpx(500);

  // 3. Re-draw PzsVsPt_PbPbComp with the scaled curve
  TFile *fileOO = TFile::Open(dir + fileNameOOPt);
  TH1 *hOO = GetHisto(fileOO, isFromFit ? "histoPzs2" : "histoPzs2NoFit");
  TH1 *hOOSyst = GetHisto(fileOO, "fHistPzsSist");
  TFile *fileRun2 = TFile::Open(dir + "Run2Results/HEPData-ins1891389-PzVsPt.root");
  TDirectoryFile *dirRun2 = fileRun2 ? (TDirectoryFile *)fileRun2->Get("P_z vs. pT (5.02TeV, 30-50%)") : nullptr;
  TFile *fileRun3 = TFile::Open(dir + "LambdaJunlee/PolIntgOut.root");
  if (!hOO || !hOOSyst || !dirRun2 || !fileRun3)
    return;
  TH1 *hRun2 = (TH1 *)dirRun2->Get("Hist1D_y1");
  TH1 *hRun2Stat = (TH1 *)dirRun2->Get("Hist1D_y1_e1");
  TH1 *hRun2Syst = (TH1 *)dirRun2->Get("Hist1D_y1_e2");
  TGraphErrors *gRun3 = (TGraphErrors *)fileRun3->Get("PtPolRun3_stat");
  TGraphErrors *gRun3Syst = (TGraphErrors *)fileRun3->Get("PtPolRun3_syst");
  if (!hRun2 || !hRun2Stat || !hRun2Syst || !gRun3 || !gRun3Syst)
  {
    std::cout << "Missing Pb-Pb inputs" << std::endl;
    return;
  }
  TH1 *hRun2Sist = (TH1 *)hRun2->Clone("hRun2Sist");
  for (Int_t i = 1; i <= hRun2->GetNbinsX(); i++)
  {
    hRun2->SetBinError(i, hRun2Stat->GetBinContent(i));
    hRun2Sist->SetBinError(i, hRun2Syst->GetBinContent(i));
  }

  TCanvas *c = new TCanvas("cScaled", "cScaled", 900, 700);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.13);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.05);
  TH1F *hFrame = new TH1F("hFrame", ";#it{p}_{T} (GeV/#it{c});#it{P}_{z,s2}", 100, 0, 6);
  hFrame->GetYaxis()->SetRangeUser(-0.001, 0.011);
  hFrame->GetYaxis()->SetTitleOffset(1.6);
  hFrame->GetXaxis()->SetTitleSize(0.05);
  hFrame->GetYaxis()->SetTitleSize(0.05);
  hFrame->GetXaxis()->SetLabelSize(0.045);
  hFrame->GetYaxis()->SetLabelSize(0.045);
  hFrame->Draw();

  for (TH1 *h : {hRun2, hRun2Sist})
  {
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.5);
    h->SetMarkerColor(kBlue);
    h->SetLineColor(kBlue);
  }
  hRun2Sist->SetFillStyle(0);
  gRun3->SetMarkerStyle(20);
  gRun3->SetMarkerSize(1.5);
  gRun3->SetMarkerColor(kGreen + 1);
  gRun3->SetLineColor(kGreen + 1);
  gRun3Syst->SetLineColor(kGreen + 1);
  gRun3Syst->SetFillStyle(0);
  for (TH1 *h : {hOO, hOOSyst})
  {
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.5);
    h->SetMarkerColor(kMagenta);
    h->SetLineColor(kMagenta);
  }
  hOOSyst->SetFillStyle(0);
  gModel502->SetLineColor(kBlack);
  gModel502->SetLineWidth(2);
  gModel536->SetLineColor(kRed);
  gModel536->SetLineWidth(2);

  hRun2->Draw("same ex0");
  hRun2Sist->Draw("same e2");
  gRun3->Draw("same p");
  gRun3Syst->Draw("same e2");
  hOO->Draw("same ex0");
  hOOSyst->Draw("same e2");
  TGraph *gPalermo = new TGraph(20, pTPalermo, Pzs2VsPtPalermo);
  gPalermo->SetLineColor(kBlue);
  gPalermo->SetLineWidth(2);
  gPalermo->Draw("same l");
  gModel536->Draw("same l");
  gModel502->Draw("same l");
  fScaled->Draw("same");

  TLegend *leg = new TLegend(0.2, 0.66, 0.7, 0.93);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(hRun2, "Pb-Pb 5.02, 30-50%, |y| < 0.5", "pl");
  leg->AddEntry(gRun3, "Pb-Pb 5.36, 30-50%, |y| < 0.5", "pl");
  leg->AddEntry(hOO, "OO 5.36, 0-50%, |#eta| < 0.8", "pl");
  leg->AddEntry(gModel502, "Model, Pb-Pb 5.02, 30-50%", "l");
  leg->AddEntry(fScaled, Form("Model #times %.2f (scaled to OO 0-50%%)", k), "l");
  leg->Draw();

  TString sout = dir + (isOOWeights ? "PzsVsPt_PbPbComp_ScaledModel_OOWeights" : "PzsVsPt_PbPbComp_ScaledModel") + outSuffix;
  c->SaveAs(sout + ".pdf");
  c->SaveAs(sout + ".png");

  TFile *fout = new TFile(sout + ".root", "RECREATE");
  fModel->Write();
  fScaled->Write();
  gModel502->Write("gModel502");
  c->Write();
  fout->Close();
  std::cout << "\nI have created " << sout << ".{pdf,png,root}" << std::endl;
}
