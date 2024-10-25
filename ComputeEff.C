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

void ComputeEff(Int_t indexMultTrial = 0, Int_t ChosenPart = ChosenParticle, Int_t RebinFactor = 1, TString inputFileNameEff = SinputFileNameEff, Int_t EtaSysChoice = ExtrEtaSysChoice, Bool_t isSysMultTrial = ExtrisSysMultTrial)
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
  SinputFileReco += ".root";
  cout << "Input file: " << SinputFileReco << endl;

  TFile *inputFileReco = new TFile(SinputFileReco, "READ");
  if (!inputFileReco)
  {
    cout << "File " << SinputFileReco << " not found" << endl;
    return;
  }
  TH2F *histoRecoBeforeBDT = (TH2F *)inputFileReco->Get("PtvsCent_BefBDT");
  TH2F *histoRecoAfterBDT = (TH2F *)inputFileReco->Get("PtvsCent_AftBDT");
  TH1F *histoRecoPtAfterBDT[numCent + 1];
  TH1F *histoRecoPtBeforeBDT[numCent + 1];
  TH1F *histoPtEff[numCent + 1]; //reco after BDT / gen
  TH1F *histoPtBDTEff[numCent + 1]; //reco after BDT / reco before BDT

  TString SinputFileGen = "FileForEfficiency/AnalysisResults_" + inputFileNameEff + ".root";
  TFile *inputFileGen = new TFile(SinputFileGen, "READ");
  if (!inputFileGen)
  {
    cout << "File " << SinputFileGen << " not found" << endl;
    return;
  }
  TDirectoryFile *dir = (TDirectoryFile *)inputFileGen->Get("lf-cascade-flow/histosMCGen");
  if (!dir)
  {
    cout << "Directory histosMCGen not found" << endl;
    return;
  }
  TH2D *histoGen = (TH2D *)dir->Get("h2DGen" + ParticleName[part] + "Eta08");
  if (!histoGen)
  {
    cout << "Histogram h2DGen" << ParticleName[part] << " not found" << endl;
    return;
  }
  TH1F *histoGenPt[numCent + 1];

  // project spectra in multiplicity classes
  TCanvas *cGen = new TCanvas("cGen", "cGen", 800, 600);
  StyleCanvas(cGen, 0.1, 0.05, 0.02, 0.1);
  TCanvas *cReco = new TCanvas("cReco", "cReco", 800, 600);
  StyleCanvas(cReco, 0.1, 0.05, 0.02, 0.1);
  TCanvas *cEff = new TCanvas("cEff", "cEff", 800, 600);
  StyleCanvas(cEff, 0.1, 0.05, 0.02, 0.1);
  TCanvas *cEffBDT = new TCanvas("cEffBDT", "cEffBDT", 800, 600);
  StyleCanvas(cEffBDT, 0.1, 0.05, 0.02, 0.1);
  Int_t CentFT0CMax = 0;
  Int_t CentFT0CMin = 0;
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
    histoGen->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoGenPt[m] = (TH1F *)histoGen->ProjectionY(Form("histoGenEta08Pt_%i-%i", CentFT0CMin, CentFT0CMax));
    histoGenPt[m]->Rebin(RebinFactor);
    StyleHistoYield(histoGenPt[m], 0, histoGenPt[numCent]->GetBinContent(histoGenPt[numCent]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1, 1.15, 1.6);
    cGen->cd();
    histoGenPt[m]->Draw("same");

    histoRecoBeforeBDT->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoRecoPtBeforeBDT[m] = (TH1F *)histoRecoBeforeBDT->ProjectionY(Form("histoRecoPtBeforeBDT_%i-%i", CentFT0CMin, CentFT0CMax));
    histoRecoPtBeforeBDT[m]->Rebin(RebinFactor);
    StyleHistoYield(histoRecoPtBeforeBDT[m], 0, histoGenPt[numCent]->GetBinContent(histoGenPt[numCent]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1, 1.15, 1.6);
    cReco->cd();
    histoRecoPtBeforeBDT[m]->Draw("same");

    histoRecoAfterBDT->GetXaxis()->SetRangeUser(CentFT0CMin + 0.001, CentFT0CMax - 0.001);
    histoRecoPtAfterBDT[m] = (TH1F *)histoRecoAfterBDT->ProjectionY(Form("histoRecoPtAfterBDT_%i-%i", CentFT0CMin, CentFT0CMax));
    histoRecoPtAfterBDT[m]->Rebin(RebinFactor);
    StyleHistoYield(histoRecoPtAfterBDT[m], 0, histoGenPt[numCent]->GetBinContent(histoGenPt[numCent]->GetMaximumBin()), ColorMult[m], MarkerMult[m], TitleXPt, "", "", 1, 1.15, 1.6);
    //histoRecoPtAfterBDT[m]->Draw("same");

    histoPtEff[m] = (TH1F *)histoRecoPtAfterBDT[m]->Clone(Form("histoPtEff_%i-%i", CentFT0CMin, CentFT0CMax));
    histoPtEff[m]->Divide(histoGenPt[m]);
    StyleHistoYield(histoPtEff[m], 0, 1, ColorMult[m], MarkerMult[m], TitleXPt, "Efficiency", "", 1, 1.15, 1.6);
    histoPtEff[m]->GetXaxis()->SetRangeUser(0., 10.);
    cEff->cd();
    histoPtEff[m]->Draw("same");

    histoPtBDTEff[m] = (TH1F *)histoRecoPtAfterBDT[m]->Clone(Form("histoPtBDTEff_%i-%i", CentFT0CMin, CentFT0CMax));
    histoPtBDTEff[m]->Divide(histoRecoPtBeforeBDT[m]);
    StyleHistoYield(histoPtBDTEff[m], 0, 1, ColorMult[m], MarkerMult[m], TitleXPt, "BDT efficiency", "", 1, 1.15, 1.6);
    histoPtBDTEff[m]->GetXaxis()->SetRangeUser(0., 10.);
    cEffBDT->cd();
    histoPtBDTEff[m]->Draw("same");
  } // end loop on centrality

  TString SOutputFile = "Efficiency/Efficiency_" + inputFileNameEff + "_" + ParticleName[part] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (!useCommonBDTValue)
    SOutputFile += "_BDTCentDep";
  if (isRun2Binning)
    SOutputFile += "_Run2Binning";
  SOutputFile += ".root";
  TFile *file = new TFile(SOutputFile, "RECREATE");
  for (Int_t cent = 0; cent < numCent + 1; cent++)
  {
    for (Int_t pt = 0; pt < numPtBins + 1; pt++)
    {
    }
  }
  file->Close();
  cout << "I created the file " << file->GetName() << endl;
}
