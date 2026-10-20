// Pz,s2 vs pT for Pb-Pb 5.02 TeV 30-50% (Run 2), Pb-Pb 5.36 TeV 30-50% (Run 3) and OO 5.36 TeV 0-50% (this analysis),
// with the Pb-Pb model predictions (top panel), and the ratio OO / Pb-Pb for each Pb-Pb dataset (bottom panel).
// Each ratio is fitted with a pol0, using stat. uncertainties only. Bins where the Pb-Pb value is compatible with zero within
// minSigmaDen stat. sigma are drawn with open markers and excluded from the fit, since the ratio is ill-defined there.
// As a cross-check, a scale factor k is also fitted directly, minimising sum (OO_i - k PbPb_i)^2 / (sOO_i^2 + k^2 sPbPb_i^2)
// over all bins, which avoids dividing by values close to zero.
// If isPbPbOverOO = 1, the inverse ratio Pb-Pb / OO is shown instead (the denominator compatibility check and the scale-factor
// fit, PbPb = k OO, are then done with OO as denominator).
//
// Usage: root -l -b -q PzsVsPtRatioToPbPb.C
//        root -l -b -q 'PzsVsPtRatioToPbPb.C(1)'   // fit with stat. and syst. uncertainties added in quadrature
//        root -l -b -q 'PzsVsPtRatioToPbPb.C(0, 1., 0, "OutputAnalysis/<Cent0-50 fit file>.root", "Pzs2VsPt/<pT file, NoFit>.root", 1, "_NoFit")'
//        root -l -b -q 'PzsVsPtRatioToPbPb.C(0, 1., 1, "OutputAnalysis/<Cent0-50 fit file>.root", "Pzs2VsPt/<pT file>.root", 1, "", 1)'   // Pb-Pb / OO
//
// The OO points are read from fileNameOOPt (output of PzsVsPt.C): histoPzs2 if isFromFit = 1, histoPzs2NoFit otherwise.
// The OO 0-50% pT-integrated value is read from fileNameOOInt (output of FitV2orPol_CB.C for Cent0-50):
// histoPzs2PtInt if isFromFit = 1, histoPzs2PtIntNoFit / histoPurityPtInt if isFromFit = 0 (as in PzsVsCentrality.C).
// If isDrawModelOO = 1, the model prediction for OO 0-50% is also drawn: the Pb-Pb 5.02 TeV 30-50% prediction, parametrised as
// f(pT) = A pT^2 / (pT^2 + B)^C exp(-D pT), scaled by k = PzsOOInt / <f>_w, with <f>_w the average of f weighted by the
// Pb-Pb 30-50% Lambda + AntiLambda spectrum for pT > 0.5 GeV/c (same procedure as in ScaledModelVsPt.C).

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TDirectoryFile.h"
#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"
#include <fstream>
#include <iostream>
#include <vector>

struct PtSet
{
  TString label;
  Int_t color;
  Int_t marker;
  std::vector<double> x, ex, y, stat, syst;
};

void PzsVsPtRatioToPbPb(Bool_t isFitWithSyst = 0,
                        Double_t minSigmaDen = 1.,
                        Bool_t isFromFit = 1,
                        TString fileNameOOInt = "OutputAnalysis/FitPzs2_LHC25_OO_pass2_Train742311_Lambda_DSCB_BkgCheb_Cent0-50_CentWeighted_EffW_Eta08_TightMassCut2.1_ReducedPtBins_AcceptanceInMacro_ResoOnTheFly_Rebin.root",
                        TString fileNameOOPt = "Pzs2VsPt/Pzs2_LHC25_OO_pass2_Train742311_Lambda_DSCB_BkgCheb_Pzs2_CentWeighted_Eta08_TightMassCut2.1_ReducedPtBins_AcceptanceInMacro_ResoOnTheFly_EffW_Rebin.root",
                        Bool_t isDrawModelOO = 1, // draw the model scaled to the OO 0-50% pT-integrated value
                        TString outSuffix = "",   // appended to the output file names
                        Bool_t isPbPbOverOO = 0,  // 0: ratio OO / Pb-Pb, 1: ratio Pb-Pb / OO
                        TString dir = "/Users/mbp-cdm-01/Desktop/ResearchTeachingActivity/MyRun3Analyses/CascadesFlow/")
{
  gStyle->SetOptStat(0);
  const Int_t nPt = 7; // 0.5-6 GeV/c, common binning of the three datasets
  const Double_t ptEdges[nPt + 1] = {0.5, 1, 1.5, 2, 2.5, 3, 4, 6};

  TFile *fRun2 = TFile::Open(dir + "Run2Results/HEPData-ins1891389-PzVsPt.root");
  TDirectoryFile *dRun2 = fRun2 ? (TDirectoryFile *)fRun2->Get("P_z vs. pT (5.02TeV, 30-50%)") : nullptr;
  TFile *fJ = TFile::Open(dir + "LambdaJunlee/PolIntgOut.root");
  TString histNameOO = isFromFit ? "histoPzs2" : "histoPzs2NoFit";
  TFile *fOO = TFile::Open(dir + fileNameOOPt);
  std::cout << "OO input: " << fileNameOOPt << " (" << histNameOO << ")" << std::endl;
  if (!dRun2 || !fJ || !fOO)
  {
    std::cout << "Missing input files" << std::endl;
    return;
  }
  TH1 *hR2 = (TH1 *)dRun2->Get("Hist1D_y1"), *hR2Stat = (TH1 *)dRun2->Get("Hist1D_y1_e1"), *hR2Syst = (TH1 *)dRun2->Get("Hist1D_y1_e2");
  TGraphErrors *gR3 = (TGraphErrors *)fJ->Get("PtPolRun3_stat"), *gR3Syst = (TGraphErrors *)fJ->Get("PtPolRun3_syst");
  TH1 *hOO = (TH1 *)fOO->Get(histNameOO), *hOOSyst = (TH1 *)fOO->Get("fHistPzsSist");
  if (!hR2 || !hR2Stat || !hR2Syst || !gR3 || !gR3Syst || !hOO || !hOOSyst)
  {
    std::cout << "Missing input histograms" << std::endl;
    return;
  }

  // 0: Pb-Pb 5.02, 1: Pb-Pb 5.36, 2: OO
  std::vector<PtSet> d(3);
  d[0] = {"Pb-Pb 5.02 TeV, 30-50%, |#it{y}| < 0.5", kBlue, 20};
  d[1] = {"Pb-Pb 5.36 TeV, 30-50%, |#it{y}| < 0.5", kGreen + 2, 21};
  d[2] = {"OO 5.36 TeV, 0-50%, |#it{#eta}| < 0.8", kMagenta, 20};
  for (Int_t b = 1; b <= nPt; b++)
  {
    d[0].x.push_back(hR2->GetBinCenter(b));
    d[0].ex.push_back(hR2->GetBinWidth(b) / 2);
    d[0].y.push_back(hR2->GetBinContent(b));
    d[0].stat.push_back(hR2Stat->GetBinContent(b));
    d[0].syst.push_back(hR2Syst->GetBinContent(b));
    d[2].x.push_back(hOO->GetBinCenter(b));
    d[2].ex.push_back(hOO->GetBinWidth(b) / 2);
    d[2].y.push_back(hOO->GetBinContent(b));
    d[2].stat.push_back(hOO->GetBinError(b));
    d[2].syst.push_back(hOOSyst->GetBinError(b));
  }
  for (Int_t i = 0; i < gR3->GetN(); i++)
  {
    d[1].x.push_back(gR3->GetX()[i]);
    d[1].ex.push_back(gR3Syst->GetEX()[i]);
    d[1].y.push_back(gR3->GetY()[i]);
    d[1].stat.push_back(gR3->GetEY()[i]);
    d[1].syst.push_back(gR3Syst->GetEY()[i]);
  }
  for (auto &s : d)
  {
    if ((Int_t)s.x.size() != nPt)
    {
      std::cout << s.label << ": expected " << nPt << " pT bins, found " << s.x.size() << std::endl;
      return;
    }
    for (Int_t i = 0; i < nPt; i++)
      if (TMath::Abs(s.x[i] - s.ex[i] - ptEdges[i]) > 1e-3 || TMath::Abs(s.x[i] + s.ex[i] - ptEdges[i + 1]) > 1e-3)
      {
        std::cout << s.label << ": pT bin " << i << " does not match the common binning" << std::endl;
        return;
      }
  }

  // Model predictions
  const char *modelFile[2] = {"Fig_5_Prediction_1-1-2_model_PbPb_5020_30-50_LocalSpinPol_z_s2_vs_pT.dat",
                              "Fig_5_Prediction_1-1-2_model_PbPb_5360_30-50_LocalSpinPol_z_s2_vs_pT.dat"};
  const char *modelLabel[2] = {"Model, Pb-Pb 5.02 TeV, 30-50%", "Model, Pb-Pb 5.36 TeV, 30-50%"};
  const Int_t modelColor[2] = {kBlue, kGreen + 2};
  TGraph *gModel[2];
  for (Int_t m = 0; m < 2; m++)
  {
    std::ifstream in((dir + modelFile[m]).Data());
    if (!in)
    {
      std::cout << "Cannot open " << modelFile[m] << std::endl;
      return;
    }
    gModel[m] = new TGraph();
    Double_t x, y;
    while (in >> x >> y)
      gModel[m]->AddPoint(x, y / 100.); // file is in %
    gModel[m]->SetLineColor(modelColor[m]);
    gModel[m]->SetLineWidth(2);
  }

  // Data graphs (small x offsets so the points do not overlap)
  const Double_t shift[3] = {-0.06, 0., 0.06};
  TGraphErrors *gD[3], *gDSyst[3];
  for (Int_t k = 0; k < 3; k++)
  {
    gD[k] = new TGraphErrors(nPt);
    gDSyst[k] = new TGraphErrors(nPt);
    for (Int_t i = 0; i < nPt; i++)
    {
      gD[k]->SetPoint(i, d[k].x[i] + shift[k], d[k].y[i]);
      gD[k]->SetPointError(i, 0, d[k].stat[i]);
      gDSyst[k]->SetPoint(i, d[k].x[i] + shift[k], d[k].y[i]);
      gDSyst[k]->SetPointError(i, 0.08, d[k].syst[i]);
    }
    gD[k]->SetMarkerStyle(d[k].marker);
    gD[k]->SetMarkerSize(1.4);
    gD[k]->SetMarkerColor(d[k].color);
    gD[k]->SetLineColor(d[k].color);
    gDSyst[k]->SetLineColor(d[k].color);
    gDSyst[k]->SetFillStyle(0);
  }

  // Ratios OO / Pb-Pb (or Pb-Pb / OO if isPbPbOverOO) and pol0 fits
  TGraphErrors *gRFit[2], *gRExcl[2], *gRSyst[2];
  TF1 *fPol0[2];
  for (Int_t p = 0; p < 2; p++)
  {
    gRFit[p] = new TGraphErrors();
    gRExcl[p] = new TGraphErrors();
    gRSyst[p] = new TGraphErrors();
    const PtSet &sN = isPbPbOverOO ? d[p] : d[2]; // numerator
    const PtSet &sD = isPbPbOverOO ? d[2] : d[p]; // denominator
    std::cout << "\n" << sN.label << " / " << sD.label << std::endl;
    for (Int_t i = 0; i < nPt; i++)
    {
      Double_t num = sN.y[i], den = sD.y[i];
      Double_t r = num / den;
      Double_t eStat = TMath::Abs(r) * sqrt(pow(sN.stat[i] / num, 2) + pow(sD.stat[i] / den, 2));
      Double_t eSyst = TMath::Abs(r) * sqrt(pow(sN.syst[i] / num, 2) + pow(sD.syst[i] / den, 2));
      Double_t eFit = isFitWithSyst ? sqrt(eStat * eStat + eSyst * eSyst) : eStat;
      Bool_t isExcluded = TMath::Abs(den) < minSigmaDen * sD.stat[i];
      Double_t xr = d[p].x[i] + (p == 0 ? -0.04 : 0.04);
      TGraphErrors *g = isExcluded ? gRExcl[p] : gRFit[p];
      g->AddPoint(xr, r);
      g->SetPointError(g->GetN() - 1, 0, eFit);
      gRSyst[p]->AddPoint(xr, r);
      gRSyst[p]->SetPointError(gRSyst[p]->GetN() - 1, 0.08, eSyst);
      std::cout << Form("  pT %4.2f: %7.3f +- %.3f (stat) +- %.3f (syst)%s", d[p].x[i], r, eStat, eSyst,
                        isExcluded ? Form("   excluded: denominator is %.1f sigma from zero", TMath::Abs(den) / sD.stat[i]) : "")
                << std::endl;
    }
    fPol0[p] = new TF1(Form("fPol0_%d", p), "pol0", 0.5, 6);
    fPol0[p]->SetLineColor(d[p].color);
    fPol0[p]->SetLineStyle(1); // solid, the dashed line is the pT-integrated ratio
    fPol0[p]->SetLineWidth(2);
    gRFit[p]->Fit(fPol0[p], "Q0");
    Double_t chi2 = fPol0[p]->GetChisquare();
    Int_t ndf = fPol0[p]->GetNDF();
    std::cout << Form("  pol0 fit (%s): %.3f +- %.3f, chi2/ndf = %.2f/%d, p-value = %.3f", isFitWithSyst ? "stat+syst" : "stat",
                      fPol0[p]->GetParameter(0), fPol0[p]->GetParError(0), chi2, ndf, TMath::Prob(chi2, ndf))
              << std::endl;

    // Cross-check: direct fit of the scale factor numerator = k denominator on all bins
    auto chi2k = [&](Double_t k)
    {
      Double_t c = 0;
      for (Int_t i = 0; i < nPt; i++)
      {
        Double_t sNum = sN.stat[i], sDen = sD.stat[i];
        if (isFitWithSyst)
        {
          sNum = sqrt(sNum * sNum + pow(sN.syst[i], 2));
          sDen = sqrt(sDen * sDen + pow(sD.syst[i], 2));
        }
        c += pow(sN.y[i] - k * sD.y[i], 2) / (sNum * sNum + k * k * sDen * sDen);
      }
      return c;
    };
    ROOT::Math::Functor1D func(chi2k);
    ROOT::Math::BrentMinimizer1D minim;
    minim.SetFunction(func, 0., 20.);
    minim.Minimize(1000, 1e-8, 1e-8);
    Double_t kBest = minim.XMinimum(), chi2Min = minim.FValMinimum();
    // 1-sigma interval from chi2 = chi2min + 1
    Double_t kLo = kBest, kHi = kBest;
    while (kLo > 0 && chi2k(kLo) < chi2Min + 1)
      kLo -= 1e-4;
    while (kHi < 20 && chi2k(kHi) < chi2Min + 1)
      kHi += 1e-4;
    std::cout << Form("  scale-factor fit on all %d bins: k = %.3f +%.3f -%.3f, chi2/ndf = %.2f/%d, p-value = %.3f", nPt, kBest,
                      kHi - kBest, kBest - kLo, chi2Min, nPt - 1, TMath::Prob(chi2Min, nPt - 1))
              << std::endl;
  }

  // Drawing
  TCanvas *c = new TCanvas("cRatio", "cRatio", 800, 900);
  TPad *pTop = new TPad("pTop", "", 0, 0.35, 1, 1);
  TPad *pBot = new TPad("pBot", "", 0, 0, 1, 0.35);
  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.15);
  pTop->SetRightMargin(0.03);
  pTop->SetTopMargin(0.05);
  pBot->SetTopMargin(0.03);
  pBot->SetBottomMargin(0.3);
  pBot->SetLeftMargin(0.15);
  pBot->SetRightMargin(0.03);
  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  TH1F *hTop = new TH1F("hTop", ";;#it{P}_{z,s2}", 100, 0, 6.2);
  hTop->GetYaxis()->SetRangeUser(-0.001, 0.0105);
  hTop->GetYaxis()->SetTitleSize(0.065);
  hTop->GetYaxis()->SetTitleOffset(1.05);
  hTop->GetYaxis()->SetLabelSize(0.05);
  hTop->GetYaxis()->SetMaxDigits(3);
  hTop->GetXaxis()->SetLabelSize(0);
  hTop->SetLineWidth(0);
  hTop->Draw("axis");
  TLine *l0 = new TLine(0, 0, 6.2, 0);
  l0->SetLineStyle(2);
  l0->Draw();
  for (Int_t m = 0; m < 2; m++)
  {
    if (m == 0)
      continue;
    gModel[m]->Draw("l same");
  }

  // Model prediction for OO 0-50%: Pb-Pb 5.02 TeV prediction scaled to the OO pT-integrated value
  TF1 *fModelOO = nullptr;
  Double_t kOO = 0;
  Double_t PzsOOInt = 0;
  if (isDrawModelOO)
  {
    TFile *fOOInt = TFile::Open(dir + fileNameOOInt);
    TH1 *hOOInt = fOOInt ? (TH1 *)fOOInt->Get(isFromFit ? "histoPzs2PtInt" : "histoPzs2PtIntNoFit") : nullptr;
    TH1 *hOOPur = fOOInt ? (TH1 *)fOOInt->Get("histoPurityPtInt") : nullptr;
    if (!hOOInt || (!isFromFit && !hOOPur))
    {
      std::cout << "Missing OO pT-integrated histograms in " << fileNameOOInt << std::endl;
      return;
    }
    PzsOOInt = hOOInt->GetBinContent(1) / (isFromFit ? 1. : hOOPur->GetBinContent(1));
    Double_t PzsOOIntErr = hOOInt->GetBinError(1) / (isFromFit ? 1. : hOOPur->GetBinContent(1));
    std::cout << "\nOO pT-integrated input: " << fileNameOOInt << " (" << (isFromFit ? "fit" : "NoFit / purity") << ")" << std::endl;
    std::cout << Form("Pzs2(OO, 0-50%%) = %.6f +- %.6f (stat)", PzsOOInt, PzsOOIntErr) << std::endl;

    TGraphErrors *gFit = new TGraphErrors(gModel[0]->GetN());
    for (Int_t i = 0; i < gModel[0]->GetN(); i++)
    {
      gFit->SetPoint(i, gModel[0]->GetX()[i], gModel[0]->GetY()[i]);
      gFit->SetPointError(i, 0, 0.01 * gModel[0]->GetY()[i]); // equal relative weight to every point
    }
    TF1 *fModel = new TF1("fModel", "[0]*x*x/pow(x*x+[1],[2])*exp(-[3]*x)", 0, 10);
    fModel->SetParameters(1.5e-3, 1., 1.2, 0.05);
    for (Int_t it = 0; it < 3; it++)
      gFit->Fit(fModel, "QN0M", "", 0.1, 5.5);
    TFile *fL = TFile::Open(dir + "CorrectedSpectra_Lambda_withEvtLoss_withFeeddown.root");
    TFile *fA = TFile::Open(dir + "CorrectedSpectra_Alambda_withEvtLoss_withFeeddown.root");
    TH1 *hW = nullptr;
    for (TString cl : {"30_40", "40_50"})
      for (TH1 *h : {fL ? (TH1 *)fL->Get("hLambdaCorrected1_CENT_" + cl) : nullptr, fA ? (TH1 *)fA->Get("hAlambdaCorrected1_CENT_" + cl) : nullptr})
      {
        if (!h)
        {
          std::cout << "Missing Pb-Pb spectrum " << cl << std::endl;
          return;
        }
        if (!hW)
          hW = (TH1 *)h->Clone("hWeights");
        else
          hW->Add(h);
      }
    Double_t num = 0, den = 0;
    for (Int_t b = 1; b <= hW->GetNbinsX(); b++)
    {
      Double_t lo = TMath::Max(hW->GetBinLowEdge(b), 0.5), hi = hW->GetBinLowEdge(b + 1);
      if (hi <= lo)
        continue;
      Double_t w = hW->GetBinContent(b) * (hi - lo);
      num += w * fModel->Integral(lo, hi) / (hi - lo);
      den += w;
    }
    kOO = PzsOOInt / (num / den);
    std::cout << Form("\nModel for OO: <f>_w = %.6f, k = %.6f / %.6f = %.3f", num / den, PzsOOInt, num / den, kOO) << std::endl;
    fModelOO = new TF1("fModelOO", Form("%.10g*fModel", kOO), 0, 6);
    fModelOO->SetLineColor(kMagenta + 1);
    fModelOO->SetLineStyle(2);
    fModelOO->SetLineWidth(3);
    fModelOO->SetNpx(500);
    // fModelOO->Draw("same");
  }
  TLegend *leg = new TLegend(0.19, 0.52, 0.65, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.042);
  for (Int_t k = 0; k < 3; k++)
  {
    gDSyst[k]->Draw("2 same");
    gD[k]->Draw("p same");
    leg->AddEntry(gD[k], d[k].label, "pe");
  }
  for (Int_t m = 0; m < 2; m++)
  {
    if (m == 0)
      continue;
    leg->AddEntry(gModel[m], modelLabel[m], "l");
  }
  // if (fModelOO)
  //   leg->AddEntry(fModelOO, Form("Model #times %.2f, scaled to OO 0-50%%", kOO), "l");
  leg->Draw();

  pBot->cd();
  TH1F *hBot = new TH1F("hBot", isPbPbOverOO ? ";#it{p}_{T} (GeV/#it{c});Pb-Pb / OO" : ";#it{p}_{T} (GeV/#it{c});OO / Pb-Pb", 100, 0, 6.2);
  if (isPbPbOverOO)
    hBot->GetYaxis()->SetRangeUser(-0.5, 2.5);
  else
    hBot->GetYaxis()->SetRangeUser(-2, 11);
  hBot->GetYaxis()->SetNdivisions(505);
  hBot->GetYaxis()->SetTitleSize(0.1);
  hBot->GetYaxis()->SetTitleOffset(0.62);
  hBot->GetYaxis()->SetLabelSize(0.09);
  hBot->GetXaxis()->SetTitleSize(0.11);
  hBot->GetXaxis()->SetLabelSize(0.09);
  hBot->SetLineWidth(0);
  hBot->Draw("axis");

  const char *denLabel[2] = {"Pb-Pb 5.02", "Pb-Pb 5.36"};
  TString ratioLabel[2];
  for (Int_t p = 0; p < 2; p++)
    ratioLabel[p] = isPbPbOverOO ? Form("%s / OO", denLabel[p]) : Form("OO / %s", denLabel[p]);

  // pT-integrated ratios (Pb-Pb 30-50% pT-integrated values)
  const Double_t PzsPbPbInt[2] = {0.000882, 0.000683}; // 5.02 TeV, 5.36 TeV
  TF1 *fRatioInt[2] = {nullptr, nullptr};
  if (PzsOOInt != 0) // only available if isDrawModelOO
    for (Int_t p = 0; p < 2; p++)
    {
      fRatioInt[p] = new TF1(Form("fRatioInt_%d", p), "pol0", 0, 6);
      fRatioInt[p]->SetParameter(0, isPbPbOverOO ? PzsPbPbInt[p] / PzsOOInt : PzsOOInt / PzsPbPbInt[p]);
      fRatioInt[p]->SetLineColor(d[p].color);
      fRatioInt[p]->SetLineStyle(2);
      fRatioInt[p]->SetLineWidth(3);
      fRatioInt[p]->Draw("same");
    }

  TLegend *legR = new TLegend(0.18, 0.7, 0.9, 0.96);
  legR->SetBorderSize(0);
  legR->SetFillStyle(0);
  legR->SetTextSize(0.065);
  for (Int_t p = 0; p < 2; p++)
    if (fRatioInt[p])
      legR->AddEntry(fRatioInt[p], Form("%s, #it{p}_{T}-integrated: %.2f", ratioLabel[p].Data(), fRatioInt[p]->GetParameter(0)), "l");

  for (Int_t p = 0; p < 2; p++)
  {
    gRSyst[p]->SetLineColor(d[p].color);
    gRSyst[p]->SetFillStyle(0);
    gRSyst[p]->Draw("2 same");
    gRFit[p]->SetMarkerStyle(d[p].marker);
    gRFit[p]->SetMarkerSize(1.3);
    gRFit[p]->SetMarkerColor(d[p].color);
    gRFit[p]->SetLineColor(d[p].color);
    gRFit[p]->Draw("p same");
    gRExcl[p]->SetMarkerStyle(d[p].marker == 20 ? 24 : 25);
    gRExcl[p]->SetMarkerSize(1.3);
    gRExcl[p]->SetMarkerColor(d[p].color);
    gRExcl[p]->SetLineColor(d[p].color);
    // if (gRExcl[p]->GetN() > 0)
    gRExcl[p]->Draw("p same");
  }
  legR->Draw();
  c->Modified();
  c->Update();

  TString sout = dir + (isPbPbOverOO ? "PzsVsPt_RatioPbPbToOO" : "PzsVsPt_RatioOOToPbPb") + (isFitWithSyst ? "_FitStatSyst" : "") + outSuffix;
  c->SaveAs(sout + ".pdf");
  c->SaveAs(sout + ".png");
  std::cout << "\nI have created " << sout << ".{pdf,png}" << std::endl;
}
