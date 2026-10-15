// Pz,s2 vs pT divided by the pT-integrated Pz,s2 of the same centrality class, for
//   Pb-Pb 5.02 TeV 30-50% (Run 2, HEPData ins1891389), Pb-Pb 5.36 TeV 30-50% (Run 3, LambdaJunlee/PolIntgOut.root)
//   and OO 5.36 TeV 0-50% (this analysis).
// Top panel: P(pT) / P_int. Bottom panel: ratio of the normalised shapes to the reference dataset (isRefOO = 1: OO, 0: Pb-Pb 5.36).
//
// The Pb-Pb 30-50% integrated values are obtained from the 30-40% and 40-50% values weighted by the Lambda + AntiLambda dN/dy
// (Pb-Pb corrected spectra, pT > 0.6 GeV/c). Stat. and syst. uncertainties of the pT points are divided by P_int; the
// uncertainty of P_int itself (normalisation) is printed and not propagated to the points.
// In the ratio panel, stat. and syst. uncertainties of the two datasets are added in quadrature.
//
// The OO points are read from fileNameOOPt (output of PzsVsPt.C): histoPzs2 if isFromFit = 1, histoPzs2NoFit otherwise.
// The OO 0-50% pT-integrated value and its stat. error are read from fileNameOOInt (output of FitV2orPol_CB.C for Cent0-50):
// histoPzs2PtInt if isFromFit = 1, histoPzs2PtIntNoFit / histoPurityPtInt if isFromFit = 0 (as in PzsVsCentrality.C).
//
// Usage: root -l -b -q NormalizedPzsVsPt.C
//        root -l -b -q 'NormalizedPzsVsPt.C(0)'
//        root -l -b -q 'NormalizedPzsVsPt.C(1, 0, "OutputAnalysis/<Cent0-50 fit file>.root", "Pzs2VsPt/<pT file, NoFit>.root", 0, "_NoFit")'

#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TDirectoryFile.h"
#include <fstream>
#include <iostream>
#include <vector>

struct PtData
{
  TString label;
  Int_t color;
  Int_t marker;
  std::vector<double> x, ex, y, stat, syst;
  double pInt, pIntStat, pIntSyst;
};

void NormalizedPzsVsPt(Bool_t isRefOO = 1,
                       Bool_t isFromFit = 1,
                       TString fileNameOOInt = "OutputAnalysis/FitPzs2_LHC25_OO_pass2_Train742311_Lambda_DSCB_BkgCheb_Cent0-50_CentWeighted_EffW_Eta08_TightMassCut2.1_ReducedPtBins_AcceptanceInMacro_ResoOnTheFly_isBkgPol0_Rebin.root",
                       TString fileNameOOPt = "Pzs2VsPt/Pzs2_LHC25_OO_pass2_Train742311_Lambda_DSCB_BkgCheb_Pzs2_CentWeighted_Eta08_TightMassCut2.1_ReducedPtBins_AcceptanceInMacro_ResoOnTheFly_isBkgPol0_EffW_Rebin.root",
                       Double_t PzsOOIntSyst = 0, // syst. uncertainty of the OO 0-50% pT-integrated value (printed only)
                       TString outSuffix = "",    // appended to the output file names
                       TString dir = "/Users/mbp-cdm-01/Desktop/ResearchTeachingActivity/MyRun3Analyses/CascadesFlow/")
{
  gStyle->SetOptStat(0);
  const Int_t nPt = 7; // 0.5-6 GeV/c, common binning of the three datasets
  const Double_t ptEdges[nPt + 1] = {0.5, 1, 1.5, 2, 2.5, 3, 4, 6};

  // Weights to combine 30-40% and 40-50% into 30-50% (Pb-Pb Lambda + AntiLambda dN/dy)
  TFile *fL = TFile::Open(dir + "CorrectedSpectra_Lambda_withEvtLoss_withFeeddown.root");
  TFile *fA = TFile::Open(dir + "CorrectedSpectra_Alambda_withEvtLoss_withFeeddown.root");
  if (!fL || !fA)
    return;
  Double_t yield[2];
  const char *cl[2] = {"30_40", "40_50"};
  TH1 *hSpectrumPbPb = nullptr; // Lambda + AntiLambda, 30-50%
  for (Int_t i = 0; i < 2; i++)
  {
    TH1 *hl = (TH1 *)fL->Get(Form("hLambdaCorrected1_CENT_%s", cl[i]));
    TH1 *ha = (TH1 *)fA->Get(Form("hAlambdaCorrected1_CENT_%s", cl[i]));
    if (!hl || !ha)
    {
      std::cout << "Missing Pb-Pb spectrum " << cl[i] << std::endl;
      return;
    }
    yield[i] = hl->Integral("width") + ha->Integral("width");
    if (!hSpectrumPbPb)
      hSpectrumPbPb = (TH1 *)hl->Clone("hSpectrumPbPb");
    else
      hSpectrumPbPb->Add(hl);
    hSpectrumPbPb->Add(ha);
  }
  Double_t w1 = yield[0] / (yield[0] + yield[1]), w2 = 1 - w1;
  std::cout << Form("Weights 30-40%% / 40-50%%: %.3f / %.3f", w1, w2) << std::endl;

  std::vector<PtData> data(3);

  // Pb-Pb 5.02 TeV, Run 2
  TFile *fRun2 = TFile::Open(dir + "Run2Results/HEPData-ins1891389-PzVsPt.root");
  TDirectoryFile *dRun2 = fRun2 ? (TDirectoryFile *)fRun2->Get("P_z vs. pT (5.02TeV, 30-50%)") : nullptr;
  TFile *fJ = TFile::Open(dir + "LambdaJunlee/PolIntgOut.root");
  if (!dRun2 || !fJ)
  {
    std::cout << "Missing Pb-Pb input files" << std::endl;
    return;
  }
  TH1 *hR2 = (TH1 *)dRun2->Get("Hist1D_y1"), *hR2Stat = (TH1 *)dRun2->Get("Hist1D_y1_e1"), *hR2Syst = (TH1 *)dRun2->Get("Hist1D_y1_e2");
  TH1 *hR2Cent = (TH1 *)fJ->Get("CentPolRun2_stat"), *hR2CentSyst = (TH1 *)fJ->Get("CentPolRun2_syst");
  TGraphErrors *gR3 = (TGraphErrors *)fJ->Get("PtPolRun3_stat"), *gR3Syst = (TGraphErrors *)fJ->Get("PtPolRun3_syst");
  TGraphErrors *gR3Cent = (TGraphErrors *)fJ->Get("CentPolRun3_stat"), *gR3CentSyst = (TGraphErrors *)fJ->Get("CentPolRun3_syst");
  if (!hR2 || !hR2Stat || !hR2Syst || !hR2Cent || !hR2CentSyst || !gR3 || !gR3Syst || !gR3Cent || !gR3CentSyst)
  {
    std::cout << "Missing Pb-Pb histograms" << std::endl;
    return;
  }

  data[0].label = "Pb-Pb 5.02 TeV, 30-50%";
  data[0].color = kBlue;
  data[0].marker = 20;
  for (Int_t b = 1; b <= nPt; b++)
  {
    data[0].x.push_back(hR2->GetBinCenter(b));
    data[0].ex.push_back(hR2->GetBinWidth(b) / 2);
    data[0].y.push_back(hR2->GetBinContent(b));
    data[0].stat.push_back(hR2Stat->GetBinContent(b));
    data[0].syst.push_back(hR2Syst->GetBinContent(b));
  }
  // Run 2 vs centrality: bins 4 (30-40) and 5 (40-50); syst. combined linearly (fully correlated)
  data[0].pInt = w1 * hR2Cent->GetBinContent(4) + w2 * hR2Cent->GetBinContent(5);
  data[0].pIntStat = sqrt(pow(w1 * hR2Cent->GetBinError(4), 2) + pow(w2 * hR2Cent->GetBinError(5), 2));
  data[0].pIntSyst = w1 * hR2CentSyst->GetBinError(4) + w2 * hR2CentSyst->GetBinError(5);

  // Pb-Pb 5.36 TeV, Run 3
  data[1].label = "Pb-Pb 5.36 TeV, 30-50%";
  data[1].color = kGreen + 2;
  data[1].marker = 21;
  for (Int_t i = 0; i < gR3->GetN(); i++)
  {
    data[1].x.push_back(gR3->GetX()[i]);
    data[1].ex.push_back(gR3Syst->GetEX()[i]);
    data[1].y.push_back(gR3->GetY()[i]);
    data[1].stat.push_back(gR3->GetEY()[i]);
    data[1].syst.push_back(gR3Syst->GetEY()[i]);
  }
  // Run 3 vs centrality: points at 35 (index 3) and 45 (index 4)
  data[1].pInt = w1 * gR3Cent->GetY()[3] + w2 * gR3Cent->GetY()[4];
  data[1].pIntStat = sqrt(pow(w1 * gR3Cent->GetEY()[3], 2) + pow(w2 * gR3Cent->GetEY()[4], 2));
  data[1].pIntSyst = w1 * gR3CentSyst->GetEY()[3] + w2 * gR3CentSyst->GetEY()[4];

  // OO 5.36 TeV, this analysis
  TString histNameOO = isFromFit ? "histoPzs2" : "histoPzs2NoFit";
  TFile *fOO = TFile::Open(dir + fileNameOOPt);
  std::cout << "OO input: " << fileNameOOPt << " (" << histNameOO << ")" << std::endl;
  TH1 *hOO = fOO ? (TH1 *)fOO->Get(histNameOO) : nullptr;
  TH1 *hOOSyst = fOO ? (TH1 *)fOO->Get("fHistPzsSist") : nullptr;
  if (!hOO || !hOOSyst)
  {
    std::cout << "Missing OO pT-differential histograms" << std::endl;
    return;
  }
  data[2].label = "OO 5.36 TeV, 0-50%";
  data[2].color = kMagenta;
  data[2].marker = 20;
  for (Int_t b = 1; b <= nPt; b++)
  {
    data[2].x.push_back(hOO->GetBinCenter(b));
    data[2].ex.push_back(hOO->GetBinWidth(b) / 2);
    data[2].y.push_back(hOO->GetBinContent(b));
    data[2].stat.push_back(hOO->GetBinError(b));
    data[2].syst.push_back(hOOSyst->GetBinError(b));
  }
  // OO 0-50% pT-integrated value
  TFile *fOOInt = TFile::Open(dir + fileNameOOInt);
  TH1 *hOOInt = fOOInt ? (TH1 *)fOOInt->Get(isFromFit ? "histoPzs2PtInt" : "histoPzs2PtIntNoFit") : nullptr;
  TH1 *hOOPur = fOOInt ? (TH1 *)fOOInt->Get("histoPurityPtInt") : nullptr;
  if (!hOOInt || (!isFromFit && !hOOPur))
  {
    std::cout << "Missing OO pT-integrated histograms in " << fileNameOOInt << std::endl;
    return;
  }
  Double_t purOO = isFromFit ? 1. : hOOPur->GetBinContent(1);
  std::cout << "OO pT-integrated input: " << fileNameOOInt << " (" << (isFromFit ? "fit" : "NoFit / purity") << ")" << std::endl;
  data[2].pInt = hOOInt->GetBinContent(1) / purOO;
  data[2].pIntStat = hOOInt->GetBinError(1) / purOO;
  data[2].pIntSyst = PzsOOIntSyst;

  // Check binning and print the normalisations
  for (auto &d : data)
  {
    if ((Int_t)d.x.size() != nPt)
    {
      std::cout << d.label << ": expected " << nPt << " pT bins, found " << d.x.size() << std::endl;
      return;
    }
    for (Int_t i = 0; i < nPt; i++)
      if (TMath::Abs(d.x[i] - d.ex[i] - ptEdges[i]) > 1e-3 || TMath::Abs(d.x[i] + d.ex[i] - ptEdges[i + 1]) > 1e-3)
      {
        std::cout << d.label << ": pT bin " << i << " does not match the common binning" << std::endl;
        return;
      }
    std::cout << Form("%-26s P_int = %.6f +- %.6f (stat) +- %.6f (syst)  -> normalisation uncertainty %.1f%% (stat)",
                      d.label.Data(), d.pInt, d.pIntStat, d.pIntSyst, 100 * d.pIntStat / d.pInt)
              << std::endl;
  }

  // Model predictions for Pb-Pb 30-50%, normalised to their own pT-integrated value
  // (average weighted by the Pb-Pb 30-50% Lambda + AntiLambda spectrum, pT > 0.5 GeV/c, within the pT range of the model)
  const char *modelFile[2] = {"Fig_5_Prediction_1-1-2_model_PbPb_5020_30-50_LocalSpinPol_z_s2_vs_pT.dat",
                              "Fig_5_Prediction_1-1-2_model_PbPb_5360_30-50_LocalSpinPol_z_s2_vs_pT.dat"};
  const char *modelLabel[2] = {"Model, Pb-Pb 5.02 TeV, 30-50%", "Model, Pb-Pb 5.36 TeV, 30-50%"};
  const Int_t modelColor[2] = {kBlack, kRed};
  TGraph *gModelN[2];
  for (Int_t m = 0; m < 2; m++)
  {
    std::ifstream in((dir + modelFile[m]).Data());
    if (!in)
    {
      std::cout << "Cannot open " << modelFile[m] << std::endl;
      return;
    }
    TGraph *g = new TGraph();
    Double_t x, y;
    while (in >> x >> y)
      g->AddPoint(x, y / 100.); // file is in %
    Double_t num = 0, den = 0;
    for (Int_t b = 1; b <= hSpectrumPbPb->GetNbinsX(); b++)
    {
      Double_t c = hSpectrumPbPb->GetBinCenter(b);
      if (c < 0.5 || c > g->GetX()[g->GetN() - 1])
        continue;
      Double_t w = hSpectrumPbPb->GetBinContent(b) * hSpectrumPbPb->GetBinWidth(b);
      num += w * g->Eval(c);
      den += w;
    }
    Double_t pIntModel = num / den;
    std::cout << Form("%-30s P_int = %.6f", modelLabel[m], pIntModel) << std::endl;
    gModelN[m] = new TGraph(g->GetN());
    for (Int_t i = 0; i < g->GetN(); i++)
      gModelN[m]->SetPoint(i, g->GetX()[i], g->GetY()[i] / pIntModel);
    gModelN[m]->SetLineColor(modelColor[m]);
    gModelN[m]->SetLineWidth(m == 0 ? 3 : 2);
    gModelN[m]->SetLineStyle(m == 0 ? 1 : 2); // the two normalised predictions almost coincide
  }

  // Normalised graphs (small x offsets so the points do not overlap)
  const Double_t shift[3] = {-0.06, 0., 0.06};
  TGraphErrors *gN[3], *gNSyst[3];
  for (Int_t k = 0; k < 3; k++)
  {
    gN[k] = new TGraphErrors(nPt);
    gNSyst[k] = new TGraphErrors(nPt);
    for (Int_t i = 0; i < nPt; i++)
    {
      Double_t r = data[k].y[i] / data[k].pInt;
      gN[k]->SetPoint(i, data[k].x[i] + shift[k], r);
      gN[k]->SetPointError(i, 0, data[k].stat[i] / data[k].pInt);
      gNSyst[k]->SetPoint(i, data[k].x[i] + shift[k], r);
      gNSyst[k]->SetPointError(i, 0.08, data[k].syst[i] / data[k].pInt);
    }
  }

  // Ratios to the reference
  Int_t iRef = isRefOO ? 2 : 1;
  TGraphErrors *gR[3] = {nullptr}, *gRSyst[3] = {nullptr};
  std::cout << "\nRatio of normalised shapes to " << data[iRef].label << std::endl;
  for (Int_t k = 0; k < 3; k++)
  {
    if (k == iRef)
      continue;
    gR[k] = new TGraphErrors(nPt);
    gRSyst[k] = new TGraphErrors(nPt);
    for (Int_t i = 0; i < nPt; i++)
    {
      Double_t a = data[k].y[i] / data[k].pInt, b = data[iRef].y[i] / data[iRef].pInt;
      Double_t ra = a / b;
      Double_t sa = data[k].stat[i] / data[k].pInt, sb = data[iRef].stat[i] / data[iRef].pInt;
      Double_t ya = data[k].syst[i] / data[k].pInt, yb = data[iRef].syst[i] / data[iRef].pInt;
      Double_t eStat = TMath::Abs(ra) * sqrt(pow(sa / a, 2) + pow(sb / b, 2));
      Double_t eSyst = TMath::Abs(ra) * sqrt(pow(ya / a, 2) + pow(yb / b, 2));
      gR[k]->SetPoint(i, data[k].x[i] + shift[k], ra);
      gR[k]->SetPointError(i, 0, eStat);
      gRSyst[k]->SetPoint(i, data[k].x[i] + shift[k], ra);
      gRSyst[k]->SetPointError(i, 0.08, eSyst);
      std::cout << Form("  %-26s pT %.2f: %7.3f +- %.3f (stat) +- %.3f (syst)", data[k].label.Data(), data[k].x[i], ra, eStat, eSyst) << std::endl;
    }
  }

  // Drawing
  TCanvas *c = new TCanvas("cNorm", "cNorm", 800, 900);
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
  TH1F *hTop = new TH1F("hTop", ";;#it{P}_{z,s2}(#it{p}_{T}) / #it{P}_{z,s2}^{#it{p}_{T}-int.}", 100, 0, 6.2);
  hTop->GetYaxis()->SetRangeUser(-0.8, 4.8);
  hTop->GetYaxis()->SetTitleSize(0.06);
  hTop->GetYaxis()->SetTitleOffset(1.1);
  hTop->GetYaxis()->SetLabelSize(0.05);
  hTop->GetXaxis()->SetLabelSize(0);
  hTop->SetLineWidth(0);
  hTop->Draw("axis");
  TLine *l1 = new TLine(0, 0, 6.2, 0);
  l1->SetLineStyle(2);
  l1->Draw();
  for (Int_t m = 0; m < 2; m++)
    gModelN[m]->Draw("l same");
  TLegend *leg = new TLegend(0.19, 0.56, 0.6, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  for (Int_t k = 0; k < 3; k++)
  {
    gNSyst[k]->SetLineColor(data[k].color);
    gNSyst[k]->SetFillStyle(0);
    gN[k]->SetLineColor(data[k].color);
    gN[k]->SetMarkerColor(data[k].color);
    gN[k]->SetMarkerStyle(data[k].marker);
    gN[k]->SetMarkerSize(1.4);
    gNSyst[k]->Draw("2 same");
    gN[k]->Draw("p same");
    leg->AddEntry(gN[k], data[k].label, "pe");
  }
  for (Int_t m = 0; m < 2; m++)
    leg->AddEntry(gModelN[m], modelLabel[m], "l");
  leg->Draw();

  pBot->cd();
  TH1F *hBot = new TH1F("hBot", Form(";#it{p}_{T} (GeV/#it{c});Ratio to %s", isRefOO ? "OO" : "Pb-Pb 5.36"), 100, 0, 6.2);
  hBot->GetYaxis()->SetRangeUser(isRefOO ? -1.5 : -3, isRefOO ? 3.5 : 6);
  hBot->GetYaxis()->SetNdivisions(505);
  hBot->GetYaxis()->SetTitleSize(0.1);
  hBot->GetYaxis()->SetTitleOffset(0.62);
  hBot->GetYaxis()->SetLabelSize(0.09);
  hBot->GetXaxis()->SetTitleSize(0.11);
  hBot->GetXaxis()->SetLabelSize(0.09);
  hBot->SetLineWidth(0);
  hBot->Draw("axis");
  TLine *l2 = new TLine(0, 1, 6.2, 1);
  l2->SetLineStyle(2);
  l2->Draw();
  for (Int_t k = 0; k < 3; k++)
  {
    if (!gR[k])
      continue;
    gRSyst[k]->SetLineColor(data[k].color);
    gRSyst[k]->SetFillStyle(0);
    gR[k]->SetLineColor(data[k].color);
    gR[k]->SetMarkerColor(data[k].color);
    gR[k]->SetMarkerStyle(data[k].marker);
    gR[k]->SetMarkerSize(1.4);
    gRSyst[k]->Draw("2 same");
    gR[k]->Draw("p same");
  }
  c->Modified();
  c->Update();

  TString sout = dir + (isRefOO ? "PzsVsPt_Normalised_RatioToOO" : "PzsVsPt_Normalised_RatioToPbPb536") + outSuffix;
  c->SaveAs(sout + ".pdf");
  c->SaveAs(sout + ".png");
  std::cout << "\nI have created " << sout << ".{pdf,png}" << std::endl;
}
