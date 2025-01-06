// This macro was originally written by:
// francesca.ercolessi@cern.ch
// and later modified by:
// chiara.de.martin@cern.ch

#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TPad.h"
#include "StyleFile.h"
#include "CommonVar.h"

void EffMultiTrial(
    Int_t mul = 0,
    TString SisSyst = "BDT",
    Int_t ChosenPart = ChosenParticle,
    Bool_t isMidRapidity = 0, // 0 for |eta| < 0.8, 1 for |y| < 0.5
    Int_t EtaSysChoice = ExtrEtaSysChoice,
    Int_t BkgType = ExtrBkgType,
    TString inputFileName = SinputFileNameEffSyst,
    Bool_t UseTwoGauss = ExtrUseTwoGauss)
{
  Int_t part = 0;
  if (ChosenPart == 1 || ChosenPart == 4 || ChosenPart == 5)
  {
    part = 1;
  }

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
  if (mul == numCent)
  { // 0-80%
    CentFT0CMin = 0;
    CentFT0CMax = 80;
  }
  else
  {
    CentFT0CMin = CentFT0C[mul];
    CentFT0CMax = CentFT0C[mul + 1];
  }

  TString histoName = Form("histoPtBDTEff_%i-%i", CentFT0CMin, CentFT0CMax);

  TString Suffix = inputFileName + Form("_%i-%i_", CentFT0C[mul], CentFT0C[mul + 1]) + ParticleName[ChosenPart] + "_" + SisSyst;

  Int_t trials = 0;
  if (SisSyst == "BDT")
    trials = trialsBDT;
  else if (SisSyst == "eta")
    trials = 2; // eta > 0 and eta < 0
  else if (SisSyst == "IR")
    trials = 5; // different interaction rates
  TString SdefBase = "Efficiency/Efficiency_" + inputFileName + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice];
  TString Sdef = SdefBase;
  if (!useCommonBDTValue)
    Sdef += "_BDTCentDep";
  if (isRun2Binning)
    Sdef += "_Run2Binning";
  Sdef += "_" + RapidityCoverage[isMidRapidity];

  TString Svaried = "";

  TFile *fdef = TFile::Open(Sdef + ".root");
  if (!fdef)
  {
    cout << "File not found: " << Sdef << endl;
    return;
  }
  cout << "Centrality class: " << CentFT0CMin << "-" << CentFT0CMax << endl;
  cout << "InputFile - default: " << Sdef << "\n"
       << endl;

  // efficiency plots
  TCanvas *cEff = new TCanvas("cEff", Form("cEff_Cent%i_%i", CentFT0CMin, CentFT0CMax), 1000, 800);
  StyleCanvas(cEff, 0.15, 0.05, 0.05, 0.15);
  cEff->cd();
  TH1F *hDefault = (TH1F *)fdef->Get(histoName);
  hDefault->SetName("hDefault");
  hDefault->SetLineColor(kBlack);
  hDefault->SetMarkerColor(kBlack);
  hDefault->Draw();

  TLegend *legTrial;
  if (SisSyst == "IR")
    legTrial = new TLegend(0.66, 0.7, 0.96, 0.9);
  else
    legTrial = new TLegend(0.66, 0.2, 0.96, 0.5);
  legTrial->SetBorderSize(0);
  legTrial->SetFillStyle(0);
  legTrial->SetTextSize(0.03);
  legTrial->SetTextFont(42);
  if (SisSyst == "BDT")
    legTrial->AddEntry(hDefault, Form("BDT score > %.3f", DefaultBDTscoreCut), "pl");
  else if (SisSyst == "eta")
    legTrial->AddEntry(hDefault, "|#eta| < 0.8", "pl");
  else if (SisSyst == "IR")
    legTrial->AddEntry(hDefault, "PbPb 2023", "pl");

  TH1F *h[trials];
  TH1F *hRatio[trials];
  Float_t BDTscoreCut = 0;
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;

    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trials * i;
    TString SBDT = Form("_BDT%.3f", BDTscoreCut);
    if (SisSyst == "BDT")
    {
      Svaried = SdefBase;
      Svaried += SBDT;
      if (!useCommonBDTValue)
        Svaried += "_BDTCentDep";
      if (isRun2Binning)
        Svaried += "_Run2Binning";
      Svaried += "_" + RapidityCoverage[isMidRapidity];
    }
    else if (SisSyst == "eta")
      Svaried = Sdef + SEtaSysChoice[i + 1];
    else if (SisSyst == "IR")
    {
      Svaried = "OutputAnalysis/FitEff_" + inputFileNameIR + SIRChoice[i + 1] + "_" + ParticleName[ChosenPart];
      Svaried += IsOneOrTwoGauss[UseTwoGauss];
      Svaried += SIsBkgParab[BkgType];
      Svaried += Form("_Cent%i-%i", CentFT0C[mul], CentFT0C[mul + 1]);
    }

    cout << "InputFile - variation: " << Svaried << endl;

    TFile *fvaried = TFile::Open(Svaried + ".root");
    if (!fvaried)
    {
      cout << "File not found: " << Svaried << endl;
      return;
    }
    TH1F *hVariedCut = (TH1F *)fvaried->Get(histoName);
    hVariedCut->SetName(histoName + Form("_%i", i));
    hRatio[i] = (TH1F *)hVariedCut->Clone(histoName + Form("Ratio_%i", i));
    hRatio[i]->Divide(hDefault);

    // efficiency plots
    cEff->cd();
    hVariedCut->SetLineColor(ColorMult[i]);
    hVariedCut->SetMarkerColor(ColorMult[i]);
    hVariedCut->SetMarkerStyle(MarkerMult[i]);
    if (SisSyst == "BDT")
      legTrial->AddEntry(hVariedCut, Form("BDT score > %.3f", BDTscoreCut), "pl");
    else if (SisSyst == "eta")
      legTrial->AddEntry(hVariedCut, SEtaSysChoice[i + 1], "pl");
    else if (SisSyst == "IR")
      legTrial->AddEntry(hVariedCut, SIRValue[i + 1], "pl");
    hVariedCut->Draw("same");
  }
  legTrial->Draw();
  cEff->SaveAs("Efficiency/EffMultTrial" + Suffix + ".pdf");
  cEff->SaveAs("Efficiency/EffMultTrial" + Suffix + ".png");

  // ratio of varied efficiency to default one
  TCanvas *cEffRatio = new TCanvas("cEffRatio", Form("cEffRatio_Cent%i_%i", CentFT0CMin, CentFT0CMax), 1000, 800);
  StyleCanvas(cEffRatio, 0.15, 0.05, 0.05, 0.15);
  cEffRatio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;
    hRatio[i]->SetLineColor(ColorMult[i]);
    hRatio[i]->SetMarkerColor(ColorMult[i]);
    hRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hRatio[i]->SetTitle("");
    hRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hRatio[i]->GetYaxis()->SetRangeUser(0, 3);
    hRatio[i]->Draw("same");
  }
  TF1 *lineat1 = new TF1("lineat1", "1", 0, 8);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineStyle(7);
  legTrial->Draw();
  lineat1->Draw("same");
  cEffRatio->SaveAs("Efficiency/EfficiencyMultTrial" + Suffix + ".pdf");
  cEffRatio->SaveAs("Efficiency/EfficiencyMultTrial" + Suffix + ".png");
}
