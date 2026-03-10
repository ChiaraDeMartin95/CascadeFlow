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
// #include "CommonVar.h"
#include "CommonVarLambda.h"
#include "StyleFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

Double_t SetEfficiencyError(Int_t k, Int_t n)
{
  return sqrt(((Double_t)k + 1) * ((Double_t)k + 2) / (n + 2) / (n + 3) - pow((Double_t)(k + 1), 2) / pow(n + 2, 2));
}

void ComputeEffLambda(Bool_t isMidRapidity = 0, // 0 for |eta| < 0.8, 1 for |y| < 0.5
                      Int_t ChosenPart = 6,     // 6 for Lambda, 8 for AntiLambda
                      TString inputFileNameEff = SinputFileNameEfficiency)
{

  TString SinputFileReco = "../FileForEfficiency/AnalysisResults_" + inputFileNameEff + ".root";
  cout << "Input file (reco) : " << SinputFileReco << endl;

  TFile *inputFileReco = new TFile(SinputFileReco, "READ");
  if (!inputFileReco)
  {
    cout << "File " << SinputFileReco << " not found" << endl;
    return;
  }

  TDirectoryFile *dirReco = (TDirectoryFile *)inputFileReco->Get("lf-cascade-flow/histosMCReco");
  TH2F *histoRecoLambda = (TH2F *)dirReco->Get("h2DRecoTrue" + ParticleName[ChosenPart]);
  if (!histoRecoLambda)
  {
    cout << "Histogram h2DRecoTrue" << ParticleName[ChosenPart] << " not found" << endl;
    return;
  }
  TH1F *histoRecoPt[numCentLambdaOO + 1];
  TH1F *histoPtEff[numCentLambdaOO + 1];
  TH1F *histoRatioTo0100[numCentLambdaOO + 1];

  TString SinputFileGen = "../FileForEfficiency/AnalysisResults_" + inputFileNameEff + ".root";
  TFile *inputFileGen = new TFile(SinputFileGen, "READ");
  if (!inputFileGen)
  {
    cout << "File " << SinputFileGen << " not found" << endl;
    return;
  }
  cout << "Input file (gen) : " << SinputFileGen << endl;
  TDirectoryFile *dirMCGen = (TDirectoryFile *)inputFileGen->Get("lf-cascade-flow/histosMCGen");
  if (!dirMCGen)
  {
    cout << "Directory histosMCGen not found" << endl;
    return;
  }
  TH2D *histoGen = (TH2D *)dirMCGen->Get("h2DGen" + ParticleName[ChosenPart] + RapidityCoverage[isMidRapidity]);
  if (!histoGen)
  {
    cout << "Histogram h2DGen" << ParticleName[ChosenPart] << " not found" << endl;
    return;
  }
  TH1F *histoGenPt[numCentLambdaOO + 1];

  Float_t NEvents[numCentLambdaOO + 1] = {0};

  TCanvas *cEff = new TCanvas("cEff", "cEff", 900, 700);
  TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 900, 700);
  StyleCanvas(cEff, 0.12, 0.05, 0.02, 0.13);
  StyleCanvas(cRatio, 0.12, 0.05, 0.02, 0.13);

  TLegend *legendMult = new TLegend(0.15, 0.7, 0.5, 0.9);
  legendMult->SetHeader("FT0C Centrality Percentile");
  legendMult->SetNColumns(3);
  legendMult->SetFillStyle(0);
  legendMult->SetTextSize(0.03);
  legendMult->SetBorderSize(0);

  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;

  for (Int_t m = numCentLambdaOO; m >= 0; m--)
  {
    if (m == (numCentLambdaOO))
    {
      CentFT0CMin = 0;
      CentFT0CMax = CentFT0CMaxLambdaOO;
    }
    else
    {
      CentFT0CMin = CentFT0CLambdaOO[m];
      CentFT0CMax = CentFT0CLambdaOO[m + 1];
    }
    cout << "\nCentrality: " << CentFT0CMin << " - " << CentFT0CMax << endl;

    histoGen->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoGenPt[m] = (TH1F *)histoGen->ProjectionY(Form("histoGen" + RapidityCoverage[isMidRapidity] + "Pt_%i-%i", CentFT0CMin, CentFT0CMax));
    histoGenPt[m] = (TH1F *)histoGenPt[m]->Rebin(numPtBinsEff, "", PtBinsEff);
    StyleHistoYield(histoGenPt[m], 0, 1.1 * histoGenPt[numCentLambdaOO]->GetBinContent(histoGenPt[numCentLambdaOO]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1.5, 1.15, 1.6);
    legendMult->AddEntry(histoGenPt[m], Form("%i-%i %%", CentFT0CMin, CentFT0CMax), "pl");
    histoGenPt[m]->GetXaxis()->SetRangeUser(MinPt[ChosenPart], MaxPt[ChosenPart]);
    
    histoRecoLambda->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoRecoPt[m] = (TH1F *)histoRecoLambda->ProjectionY(Form("histoRecoPt_%i-%i", CentFT0CMin, CentFT0CMax));
    histoRecoPt[m] = (TH1F *)histoRecoPt[m]->Rebin(numPtBinsEff, "", PtBinsEff);
    StyleHistoYield(histoRecoPt[m], 0, 1.1 * histoRecoPt[numCentLambdaOO]->GetBinContent(histoRecoPt[numCentLambdaOO]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1.5, 1.15, 1.6);
    histoRecoPt[m]->GetXaxis()->SetRangeUser(MinPt[ChosenPart], MaxPt[ChosenPart]);

    histoPtEff[m] = (TH1F *)histoRecoPt[m]->Clone(Form("histoPtEff_%i-%i", CentFT0CMin, CentFT0CMax));
    histoPtEff[m]->Divide(histoGenPt[m]);
    for (Int_t i = 1; i <= histoPtEff[m]->GetNbinsX(); i++)
    {
      histoPtEff[m]->SetBinError(i, SetEfficiencyError(histoRecoPt[m]->GetBinContent(i), histoGenPt[m]->GetBinContent(i)));
    }
    StyleHistoYield(histoPtEff[m], 0, 0.4, ColorMult[m], MarkerMult[m], TitleXPt, "Efficiency", "", 1, 1.15, 1.2);
    histoPtEff[m]->GetXaxis()->SetRangeUser(MinPt[ChosenPart], MaxPt[ChosenPart]);
    cEff->cd();
    histoPtEff[m]->Draw("same e");
    legendMult->Draw("");

    histoRatioTo0100[m] = (TH1F *)histoPtEff[m]->Clone(Form("histoRatioTo0100_%i-%i", CentFT0CMin, CentFT0CMax));
    histoRatioTo0100[m]->Divide(histoPtEff[numCentLambdaOO]);
    cRatio->cd();
    StyleHistoYield(histoRatioTo0100[m], 0.8, 1.2, ColorMult[m], MarkerMult[m], TitleXPt, "Efficiency / Efficiency 0-100%", "", 1, 1.15, 1.2);
    histoRatioTo0100[m]->GetXaxis()->SetRangeUser(MinPt[ChosenPart], MaxPt[ChosenPart]);
    histoRatioTo0100[m]->Draw("same e");
    legendMult->Draw("");

  } // end loop on centrality

  cEff->SaveAs("../Efficiency/Efficiency_" + inputFileNameEff + "_" + ParticleName[ChosenPart] + RapidityCoverage[isMidRapidity] + ".pdf");
  cEff->SaveAs("../Efficiency/Efficiency_" + inputFileNameEff + "_" + ParticleName[ChosenPart] + RapidityCoverage[isMidRapidity] + ".png");
  cRatio->SaveAs("../Efficiency/EfficiencyRatioTo0100_" + inputFileNameEff + "_" + ParticleName[ChosenPart] + RapidityCoverage[isMidRapidity] + ".pdf");
  cRatio->SaveAs("../Efficiency/EfficiencyRatioTo0100_" + inputFileNameEff + "_" + ParticleName[ChosenPart] + RapidityCoverage[isMidRapidity] + ".png");

  TString SOutputFile = "../Efficiency/Efficiency_" + inputFileNameEff + "_" + ParticleName[ChosenPart] + RapidityCoverage[isMidRapidity] + ".root";
  TFile *file = new TFile(SOutputFile, "RECREATE");
  for (Int_t m = 0; m <= numCentLambdaOO; m++)
  {
    histoPtEff[m]->Write();
  }
  file->Close();

  cout << "I created the file " << file->GetName() << endl;
}
