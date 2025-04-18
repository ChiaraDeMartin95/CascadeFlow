// This macro was originally written by:
// chiara.de.martin@cern.ch
#include "TStyle.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TPad.h"
#include "CommonVar.h"
#include "StyleFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

Double_t SetEfficiencyError(Int_t k, Int_t n)
{
  return sqrt(((Double_t)k + 1) * ((Double_t)k + 2) / (n + 2) / (n + 3) - pow((Double_t)(k + 1), 2) / pow(n + 2, 2));
}

void ComputeEff(Bool_t isMidRapidity = 0, // 0 for |eta| < 0.8, 1 for |y| < 0.5
                Int_t indexMultTrial = 0,
                Int_t ChosenPart = ChosenParticle,
                Int_t RebinFactor = 1,
                TString inputFileNameEff = SinputFileNameEff,
                Int_t EtaSysChoice = ExtrEtaSysChoice,
                Bool_t isSysMultTrial = ExtrisSysMultTrial)
{

  Int_t part = 0;
  if (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5)
  {
    part = 1;
  }

  // if (isSysMultTrial)
  // inputFileNameEff = SinputFileNameSystEff;
  Float_t BDTscoreCut = DefaultBDTscoreCut;
  if (indexMultTrial > trialsBDT)
    return;
  if (isSysMultTrial)
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trialsBDT * indexMultTrial;
  TString SBDT = "";
  if (BDTscoreCut != DefaultBDTscoreCut)
    SBDT = Form("_BDT%.3f", BDTscoreCut);
  
  TString SinputFileReco = "OutputAnalysis/Output_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (!useCommonBDTValue)
    SinputFileReco += "_BDTCentDep";
  if (isRun2Binning)
    SinputFileReco += "_Run2Binning";
  // SinputFileReco += "_EffTrain305015";
  SinputFileReco += ".root";
  cout << "Input file: " << SinputFileReco << endl;

  TFile *inputFileReco = new TFile(SinputFileReco, "READ");
  if (!inputFileReco)
  {
    cout << "File " << SinputFileReco << " not found" << endl;
    return;
  }

  TH2F *histoRecoBeforeBDT = (TH2F *)inputFileReco->Get("PtvsCent_BefBDT");
  //before applying BDT in macro (but after applying BDT in task)
  if (isMidRapidity)
    histoRecoBeforeBDT = (TH2F *)inputFileReco->Get("PtvsCent_Y05_BefBDT");
  TH2F *histoRecoAfterBDT = (TH2F *)inputFileReco->Get("PtvsCent_AftBDT");
  //after applying BDT in macro
  if (isMidRapidity)
    histoRecoAfterBDT = (TH2F *)inputFileReco->Get("PtvsCent_Y05_AftBDT");
  TH1F *histoRecoPtAfterBDT[numCent + 1];
  TH1F *histoRecoPtBeforeBDT[numCent + 1];
  TH1F *histoRecoPt[numCent + 1];
  TH1F *histoPtEff[numCent + 1];      // reco after BDT / gen
  TH1F *histoPtBDTEff[numCent + 1];   // reco after BDT / reco before any BDT selection
  TH1F *histoPtBDTMacroEff[numCent + 1]; // reco after BDT / reco before macro BDT
  TH1F *histoPtEffwoBDT[numCent + 1]; // reco before macro BDT / gen
  TH1F *histoPtRecoEff[numCent + 1];  // reco before any BDT selection / gen

  TString SinputFileGen = "FileForEfficiency/AnalysisResults_" + inputFileNameEff + ".root";
  TFile *inputFileGen = new TFile(SinputFileGen, "READ");
  if (!inputFileGen)
  {
    cout << "File " << SinputFileGen << " not found" << endl;
    return;
  }
  cout << "Reading file: " << SinputFileGen << endl;
  TDirectoryFile *dirMCGen = (TDirectoryFile *)inputFileGen->Get("lf-cascade-flow/histosMCGen");
  if (!dirMCGen)
  {
    cout << "Directory histosMCGen not found" << endl;
    return;
  }
  TH2D *histoGen = (TH2D *)dirMCGen->Get("h2DGen" + ParticleName[part] + RapidityCoverage[isMidRapidity]);
  if (!histoGen)
  {
    cout << "Histogram h2DGen" << ParticleName[part] << " not found" << endl;
    return;
  }
  TH1F *histoGenPt[numCent + 1];
  TH1F *histoGenPtPerEvent[numCent + 1];
  TDirectoryFile *dirReco = (TDirectoryFile *)inputFileGen->Get("lf-cascade-flow/histos");
  if (!dirReco)
  {
    cout << "Directory histos not found" << endl;
    return;
  }
  TH1F *histoRecovsCent = (TH1F *)dirReco->Get("hEventCentrality");
  if (!histoRecovsCent)
  {
    cout << "Histogram hEventCentrality not found" << endl;
    return;
  }
  Float_t NEvents[numCent + 1] = {0};
  //TH2F *histoReco = (TH2F *)dirReco->Get("hXiPtvsCent"); //all reco Xi from derived data, no BDT applied
  TH2F *histoReco = (TH2F *)dirReco->Get("hXiPtvsCent" + RapidityCoverage[isMidRapidity]);
  if (part == 1)
  {
    histoReco = (TH2F *)dirReco->Get("hOmegaPtvsCent" + RapidityCoverage[isMidRapidity]);
  }
  if (!histoReco)
  {
    cout << "Histogram hXiPtvsCent (or hOmegaPtvsCent) not found" << endl;
    return;
  }

  // project spectra in multiplicity classes
  TCanvas *cGen = new TCanvas("cGen", "cGen", 800, 600);
  StyleCanvas(cGen, 0.12, 0.05, 0.02, 0.13);
  TCanvas *cGenPerEvent = new TCanvas("cGenPerEvent", "cGenPerEvent", 800, 600);
  StyleCanvas(cGenPerEvent, 0.12, 0.05, 0.02, 0.13);
  TCanvas *cReco = new TCanvas("cReco", "cReco", 800, 600);
  StyleCanvas(cReco, 0.12, 0.05, 0.02, 0.13);
  TCanvas *cEffReco = new TCanvas("cEffReco", "reco before any BDT selection / gen", 800, 600);
  StyleCanvas(cEffReco, 0.12, 0.05, 0.02, 0.13);
  TCanvas *cEff = new TCanvas("cEff", "reco after BDT / gen", 800, 600);
  StyleCanvas(cEff, 0.12, 0.05, 0.02, 0.13);
  TCanvas *cEffwoBDT = new TCanvas("cEffwoBDT", "reco before any BDT selection / gen", 800, 600);
  StyleCanvas(cEffwoBDT, 0.12, 0.05, 0.02, 0.13);
  TCanvas *cEffBDT = new TCanvas("cEffBDT", "reco after BDT / reco before any BDT", 800, 600);
  StyleCanvas(cEffBDT, 0.12, 0.05, 0.02, 0.13);
  TCanvas *cEffMacroBDT = new TCanvas("cEffMacroBDT", "reco after BDT / reco before macro BDT", 800, 600);
  StyleCanvas(cEffMacroBDT, 0.12, 0.05, 0.02, 0.13);
  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  TLegend *legendMult = new TLegend(0.15, 0.7, 0.5, 0.9);
  legendMult->SetHeader("FT0C Centrality Percentile");
  legendMult->SetNColumns(3);
  legendMult->SetFillStyle(0);
  legendMult->SetTextSize(0.03);
  legendMult->SetBorderSize(0);

  for (Int_t m = numCent; m >= 0; m--)
  {
    if (m == numCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = 80;
    }
    else
    {
      CentFT0CMin = CentFT0C[m];
      CentFT0CMax = CentFT0C[m + 1];
    }
    // cout << "\nCentrality: " << CentFT0CMin << " - " << CentFT0CMax << endl;
    histoRecovsCent->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    NEvents[m] = histoRecovsCent->Integral(histoRecovsCent->FindBin(CentFT0CMin + 0.001), histoRecovsCent->FindBin(CentFT0CMax - 0.001));
    histoGen->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoGenPt[m] = (TH1F *)histoGen->ProjectionY(Form("histoGen" + RapidityCoverage[isMidRapidity] + "Pt_%i-%i", CentFT0CMin, CentFT0CMax));
    histoGenPt[m] = (TH1F *)histoGenPt[m]->Rebin(numPtBinsEff, "", PtBinsEff);
    // histoGenPt[m] = (TH1F *)histoGenPt[m]->Rebin(RebinFactor);
    StyleHistoYield(histoGenPt[m], 0, 1.1 * histoGenPt[numCent]->GetBinContent(histoGenPt[numCent]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1.5, 1.15, 1.6);
    histoGenPtPerEvent[m] = (TH1F *)histoGenPt[m]->Clone(Form("histoGenPtPerEvent_%i-%i", CentFT0CMin, CentFT0CMax));
    histoGenPtPerEvent[m]->Scale(1. / NEvents[m]);
    legendMult->AddEntry(histoGenPt[m], Form("%i-%i %%", CentFT0CMin, CentFT0CMax), "pl");
    histoGenPt[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
    cGen->cd();
    histoGenPt[m]->Draw("same e");
    cGenPerEvent->cd();
    for (Int_t b = 1; b <= histoGenPtPerEvent[m]->GetNbinsX(); b++)
    {
      histoGenPtPerEvent[m]->SetBinContent(b, histoGenPtPerEvent[m]->GetBinContent(b) / histoGenPtPerEvent[m]->GetBinWidth(b));
      histoGenPtPerEvent[m]->SetBinError(b, histoGenPtPerEvent[m]->GetBinError(b) / histoGenPtPerEvent[m]->GetBinWidth(b));
    }
    histoGenPtPerEvent[m]->GetYaxis()->SetRangeUser(0, 0.6);
    histoGenPtPerEvent[m]->Draw("same e");

    histoReco->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoRecoPt[m] = (TH1F *)histoReco->ProjectionY(Form("histoRecoPt_%i-%i", CentFT0CMin, CentFT0CMax));
    histoRecoPt[m] = (TH1F *)histoRecoPt[m]->Rebin(numPtBinsEff, "", PtBinsEff);
    // histoRecoPt[m] = (TH1F *)histoRecoPt[m]->Rebin(RebinFactor);
    StyleHistoYield(histoRecoPt[m], 0, 1.1 * histoRecoPt[numCent]->GetBinContent(histoRecoPt[numCent]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1.5, 1.15, 1.6);
    histoRecoPt[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);

    histoRecoBeforeBDT->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoRecoPtBeforeBDT[m] = (TH1F *)histoRecoBeforeBDT->ProjectionY(Form("histoRecoPtBeforeBDT_%i-%i", CentFT0CMin, CentFT0CMax));
    histoRecoPtBeforeBDT[m] = (TH1F *)histoRecoPtBeforeBDT[m]->Rebin(numPtBinsEff, "", PtBinsEff);
    // histoRecoPtBeforeBDT[m] = (TH1F *)histoRecoPtBeforeBDT[m]->Rebin(RebinFactor);
    StyleHistoYield(histoRecoPtBeforeBDT[m], 0, 1.1 * histoRecoPtBeforeBDT[numCent]->GetBinContent(histoRecoPtBeforeBDT[numCent]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1.5, 1.15, 1.6);
    histoRecoPtBeforeBDT[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
    cReco->cd();
    histoRecoPtBeforeBDT[m]->Draw("same e");

    histoRecoAfterBDT->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoRecoPtAfterBDT[m] = (TH1F *)histoRecoAfterBDT->ProjectionY(Form("histoRecoPtAfterBDT_%i-%i", CentFT0CMin, CentFT0CMax));
    histoRecoPtAfterBDT[m] = (TH1F *)histoRecoPtAfterBDT[m]->Rebin(numPtBinsEff, "", PtBinsEff);
    // histoRecoPtAfterBDT[m] = (TH1F *)histoRecoPtAfterBDT[m]->Rebin(RebinFactor);
    StyleHistoYield(histoRecoPtAfterBDT[m], 0, 1.1 * histoRecoPtAfterBDT[numCent]->GetBinContent(histoRecoPtAfterBDT[numCent]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1.5, 1.15, 1.6);
    histoRecoPtAfterBDT[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
    // histoRecoPtAfterBDT[m]->Draw("same e");

    histoPtEff[m] = (TH1F *)histoRecoPtAfterBDT[m]->Clone(Form("histoPtEff_%i-%i", CentFT0CMin, CentFT0CMax));
    histoPtEff[m]->Divide(histoGenPt[m]);
    for (Int_t i = 1; i <= histoPtEff[m]->GetNbinsX(); i++)
    {
      histoPtEff[m]->SetBinError(i, SetEfficiencyError(histoRecoPtAfterBDT[m]->GetBinContent(i), histoGenPt[m]->GetBinContent(i)));
    }
    StyleHistoYield(histoPtEff[m], 0, 0.2, ColorMult[m], MarkerMult[m], TitleXPt, "Efficiency", "", 1, 1.15, 1.2);
    histoPtEff[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
    cEff->cd();
    histoPtEff[m]->Draw("same e");

    histoPtBDTEff[m] = (TH1F *)histoRecoPtAfterBDT[m]->Clone(Form("histoPtBDTEff_%i-%i", CentFT0CMin, CentFT0CMax));
    histoPtBDTEff[m]->Divide(histoRecoPt[m]);
    for (Int_t i = 1; i <= histoPtBDTEff[m]->GetNbinsX(); i++)
    {
      histoPtBDTEff[m]->SetBinError(i, SetEfficiencyError(histoRecoPtAfterBDT[m]->GetBinContent(i), histoRecoPt[m]->GetBinContent(i)));
    }
    StyleHistoYield(histoPtBDTEff[m], 0, 1.3, ColorMult[m], MarkerMult[m], TitleXPt, "Efficiency of BDT cut", "", 1, 1.15, 1.2);
    histoPtBDTEff[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
    cEffBDT->cd();
    histoPtBDTEff[m]->Draw("same e");

    histoPtBDTMacroEff[m] = (TH1F *)histoRecoPtAfterBDT[m]->Clone(Form("histoPtBDTMacroEff_%i-%i", CentFT0CMin, CentFT0CMax));
    histoPtBDTMacroEff[m]->Divide(histoRecoPtBeforeBDT[m]);
    for (Int_t i = 1; i <= histoPtBDTMacroEff[m]->GetNbinsX(); i++)
    {
      histoPtBDTMacroEff[m]->SetBinError(i, SetEfficiencyError(histoRecoPtAfterBDT[m]->GetBinContent(i), histoRecoPtBeforeBDT[m]->GetBinContent(i)));
    }
    StyleHistoYield(histoPtBDTMacroEff[m], 0, 1.3, ColorMult[m], MarkerMult[m], TitleXPt, "Efficiency of BDT cut", "", 1, 1.15, 1.2);
    histoPtBDTMacroEff[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
    cEffMacroBDT->cd();
    histoPtBDTMacroEff[m]->Draw("same e");

    histoPtEffwoBDT[m] = (TH1F *)histoRecoPtBeforeBDT[m]->Clone(Form("histoPtEffwoBDT_%i-%i", CentFT0CMin, CentFT0CMax));
    histoPtEffwoBDT[m]->Divide(histoGenPt[m]);
    for (Int_t i = 1; i <= histoPtEffwoBDT[m]->GetNbinsX(); i++)
    {
      histoPtEffwoBDT[m]->SetBinError(i, SetEfficiencyError(histoRecoPtBeforeBDT[m]->GetBinContent(i), histoGenPt[m]->GetBinContent(i)));
    }
    StyleHistoYield(histoPtEffwoBDT[m], 0, 0.2, ColorMult[m], MarkerMult[m], TitleXPt, "Efficiency before BDT task cut", "", 1, 1.15, 1.2);
    histoPtEffwoBDT[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
    cEffwoBDT->cd();
    histoPtEffwoBDT[m]->Draw("same e");

    histoPtRecoEff[m] = (TH1F *)histoRecoPt[m]->Clone(Form("histoPtRecoEff_%i-%i", CentFT0CMin, CentFT0CMax));
    histoPtRecoEff[m]->Divide(histoGenPt[m]);
    for (Int_t i = 1; i <= histoPtRecoEff[m]->GetNbinsX(); i++)
    {
      histoPtRecoEff[m]->SetBinError(i, SetEfficiencyError(histoRecoPt[m]->GetBinContent(i), histoGenPt[m]->GetBinContent(i)));
    }
    StyleHistoYield(histoPtRecoEff[m], 0, 0.2, ColorMult[m], MarkerMult[m], TitleXPt, "Efficiency before any BDT cut", "", 1, 1.15, 1.2);
    histoPtRecoEff[m]->GetXaxis()->SetRangeUser(MinPt[part], MaxPt[part]);
    cEffReco->cd();
    histoPtRecoEff[m]->Draw("same e");
  } // end loop on centrality

  TCanvas *cEffvsMult = new TCanvas("cEffvsMult", "cEffvsMult", 1200, 800);
  StyleCanvas(cEffvsMult, 0.12, 0.05, 0.02, 0.13);
  TLegend *legendPt = new TLegend(0.45, 0.6, 0.9, 0.9);
  legendPt->SetHeader("p_{T} bins");
  legendPt->SetNColumns(2);
  legendPt->SetFillStyle(0);
  legendPt->SetTextSize(0.03);
  legendPt->SetBorderSize(0);
  TH1F *histoEffVsMult[numPtBinsEff + 1];
  TGraphAsymmErrors *gEffVsNch[numPtBinsEff + 1];
  TH1F *histoEffVsNch[numPtBinsEff + 1];
  Double_t Efficiency[numPtBinsEff + 1][numCent + 1] = {{0}};
  Double_t EffError[numPtBinsEff + 1][numCent + 1] = {{0}};
  for (Int_t pt = 0; pt < numPtBinsEff; pt++)
  {
    // cout << "Pt bin: " << PtBinsEff[pt] << " - " << PtBinsEff[pt + 1] << endl;
    histoEffVsMult[pt] = new TH1F(Form("histoEffVsMult_%i", pt), Form("Efficiency vs Centrality for %2.1f < p_{T} < %2.1f", PtBinsEff[pt], PtBinsEff[pt + 1]), numCent, 0, numCent);
    for (Int_t m = (numCent - 1); m >= 0; m--)
    {
      if (m == numCent)
      { // 0-80%
        CentFT0CMin = 0;
        CentFT0CMax = 80;
      }
      else
      {
        CentFT0CMin = CentFT0C[m];
        CentFT0CMax = CentFT0C[m + 1];
      }
      // cout << "Centrality: " << CentFT0C[m] << " - " << CentFT0C[m + 1] << endl;

      histoEffVsMult[pt]->GetXaxis()->SetBinLabel(m + 1, Form("%i-%i %%", CentFT0C[m], CentFT0C[m + 1]));
      histoEffVsMult[pt]->SetBinContent(m + 1, histoPtRecoEff[m]->GetBinContent(pt + 1));
      histoEffVsMult[pt]->SetBinError(m + 1, histoPtRecoEff[m]->GetBinError(pt + 1));
      Efficiency[pt][m] = histoPtRecoEff[m]->GetBinContent(pt + 1);
      EffError[pt][m] = histoPtRecoEff[m]->GetBinError(pt + 1);
    }
    StyleHistoYield(histoEffVsMult[pt], 0, 0.3, ColorMult[pt], MarkerMult[pt], TitleXCent, "Efficiency", "", 1, 1.15, 1.2);
    legendPt->AddEntry(histoEffVsMult[pt], Form("%2.1f < p_{T} < %2.1f GeV/c", PtBinsEff[pt], PtBinsEff[pt + 1]), "pl");
    cEffvsMult->cd();
    histoEffVsMult[pt]->Draw("same e");
    gEffVsNch[pt] = new TGraphAsymmErrors(numCent, dNdEtaAbhi, Efficiency[pt], dNdEtaAbhiErr, dNdEtaAbhiErr, EffError[pt], EffError[pt]);
    gEffVsNch[pt]->SetLineColor(ColorMult[pt]);
    gEffVsNch[pt]->SetMarkerColor(ColorMult[pt]);
    gEffVsNch[pt]->SetMarkerStyle(MarkerMult[pt]);
    // StyleTGraphErrors(gEffVsNch[pt], ColorMult[pt], MarkerMult[pt], 1.5, 1);
  }
  legendPt->Draw();

  TCanvas *cEffVsNCh = new TCanvas("cEffVsNCh", "cEffVsNCh", 1200, 800);
  StyleCanvas(cEffVsNCh, 0.12, 0.05, 0.02, 0.13);
  TF1 *fitEffVsNch[numPtBinsEff + 1];
  TH1F *hdummy = new TH1F("hdummy", "Efficiency vs Nch", 10, 0, 2500);
  StyleHistoYield(hdummy, 0, 0.2, 1, 1, "N_{ch}", "Efficiency", "", 1, 1.15, 1.2);
  hdummy->Draw();
  for (Int_t pt = 0; pt < numPtBinsEff; pt++)
  {
    fitEffVsNch[pt] = new TF1(Form("fitEffVsNch_%i", pt), "expo", 0, 2500);
    fitEffVsNch[pt]->SetLineColor(ColorMult[pt]);
    gEffVsNch[pt]->Fit(fitEffVsNch[pt], "R+");
    gEffVsNch[pt]->Draw("same p");
  }
  legendPt->Draw();

  TH1F *hpar0_Eff = new TH1F("hpar0_Eff", "hpar0_Eff", numPtBinsEff, PtBinsEff);
  TH1F *hpar1_Eff = new TH1F("hpar1_Eff", "hpar1_Eff", numPtBinsEff, PtBinsEff);
  for (Int_t pt = 0; pt < numPtBinsEff; pt++)
  {
    hpar0_Eff->SetBinContent(pt + 1, fitEffVsNch[pt]->GetParameter(0));
    hpar0_Eff->SetBinError(pt + 1, fitEffVsNch[pt]->GetParError(0));
    hpar1_Eff->SetBinContent(pt + 1, fitEffVsNch[pt]->GetParameter(1));
    hpar1_Eff->SetBinError(pt + 1, fitEffVsNch[pt]->GetParError(1));
  }

  TCanvas *Chi2VsPt = new TCanvas("Chi2VsPt", "Chi2VsPt", 1200, 800);
  StyleCanvas(Chi2VsPt, 0.12, 0.05, 0.02, 0.13);
  TH1F *histoChi2VsPt = (TH1F *)histoRecoPt[0]->Clone("histoChi2VsPt");
  for (Int_t pt = 0; pt < numPtBinsEff; pt++)
  {
    histoChi2VsPt->SetBinContent(pt + 1, fitEffVsNch[pt]->GetChisquare() / fitEffVsNch[pt]->GetNDF());
    histoChi2VsPt->SetBinError(pt + 1, 0);
  }
  StyleHistoYield(histoChi2VsPt, 0, 3, 1, 1, "p_{T}", "#chi^{2}/NDF", "", 1, 1.15, 1.2);
  histoChi2VsPt->Draw();

  cGen->cd();
  legendMult->Draw();
  cGen->SaveAs("Efficiency/GenSpectra_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".pdf");
  cGen->SaveAs("Efficiency/GenSpectra_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".png");
  cGenPerEvent->cd();
  legendMult->Draw();
  cGenPerEvent->SaveAs("Efficiency/GenSpectraPerEvent_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".pdf");
  cGenPerEvent->SaveAs("Efficiency/GenSpectraPerEvent_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".png");
  cReco->cd();
  legendMult->Draw();
  cReco->SaveAs("Efficiency/RecoSpectra_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".pdf");
  cReco->SaveAs("Efficiency/RecoSpectra_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".png");
  cEff->cd();
  legendMult->Draw();
  cEff->SaveAs("Efficiency/Efficiency_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".pdf");
  cEff->SaveAs("Efficiency/Efficiency_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".png");
  cEffBDT->cd();
  legendMult->Draw();
  cEffBDT->SaveAs("Efficiency/BDTEfficiency_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".pdf");
  cEffBDT->SaveAs("Efficiency/BDTEfficiency_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".png");
  cEffwoBDT->cd();
  legendMult->Draw();
  cEffwoBDT->SaveAs("Efficiency/EfficiencywoBDT_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".pdf");
  cEffwoBDT->SaveAs("Efficiency/EfficiencywoBDT_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".png");
  cEffReco->cd();
  legendMult->Draw();
  cEffReco->SaveAs("Efficiency/EfficiencyReco_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".pdf");
  cEffReco->SaveAs("Efficiency/EfficiencyReco_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + ".png");

  cEffvsMult->SaveAs("Efficiency/EfficiencyVsCentrality_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".pdf");
  cEffvsMult->SaveAs("Efficiency/EfficiencyVsCentrality_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".png");
  cEffVsNCh->SaveAs("Efficiency/EfficiencyVsNch_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".pdf");
  cEffVsNCh->SaveAs("Efficiency/EfficiencyVsNch_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".png");
  Chi2VsPt->SaveAs("Efficiency/Chi2VsPt_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".pdf");
  Chi2VsPt->SaveAs("Efficiency/Chi2VsPt_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".png");

  TString SOutputFile = "Efficiency/Efficiency_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (!useCommonBDTValue)
    SOutputFile += "_BDTCentDep";
  if (isRun2Binning)
    SOutputFile += "_Run2Binning";
  SOutputFile += "_" + RapidityCoverage[isMidRapidity];
  SOutputFile += ".root";
  TFile *file = new TFile(SOutputFile, "RECREATE");
  for (Int_t pt = 0; pt < numPtBinsEff; pt++)
  {
    fitEffVsNch[pt]->Write();
  }
  hpar0_Eff->Write();
  hpar1_Eff->Write();
  for (Int_t m = 0; m < numCent; m++)
  {
    histoPtBDTEff[m]->Write();
  }
  file->Close();

  TCanvas *par0Eff = new TCanvas("par0Eff", "par0Eff", 1200, 800);
  StyleCanvas(par0Eff, 0.12, 0.05, 0.02, 0.13);
  StyleHistoYield(hpar0_Eff, -4, -1, ColorMult[0], MarkerMult[0], TitleXPt, "par0", "", 1, 1.15, 1.2);
  hpar0_Eff->Draw();
  TCanvas *par1Eff = new TCanvas("par1Eff", "par1Eff", 1200, 800);
  StyleCanvas(par1Eff, 0.12, 0.05, 0.02, 0.13);
  StyleHistoYield(hpar1_Eff, -0.0015, -0.0002, ColorMult[0], MarkerMult[0], TitleXPt, "par1", "", 1, 1.15, 1.2);
  hpar1_Eff->Draw();
  par0Eff->SaveAs("Efficiency/par0Eff_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".pdf");
  par0Eff->SaveAs("Efficiency/par0Eff_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".png");
  par1Eff->SaveAs("Efficiency/par1Eff_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".pdf");
  par1Eff->SaveAs("Efficiency/par1Eff_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT + ".png");

  cout << "I created the file " << file->GetName() << endl;
}
