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

void QAEffWeights(Int_t ChosenPart = ChosenParticle,
                  TString inputFileName = SinputFileName,
                  Int_t EtaSysChoice = ExtrEtaSysChoice)
{

  // Get resolution
  TString fileResoName = "";
  if (v2type == 1)
  {
    if (SinputFileName.Index("CFW") != -1)
      fileResoName = ResoFileName_SPCFW;
    else
      fileResoName = ResoFileName_SPLF;
  }
  else
  {
    if (SinputFileName.Index("CFW") != -1)
      fileResoName = ResoFileName_EPCFW;
    else
      fileResoName = ResoFileName_EPLF;
  }
  fileResoName += ".root";
  TFile *fileResoEP = new TFile(fileResoName, "");
  TH1F *hReso = (TH1F *)fileResoEP->Get("hReso");
  TH1F *hReso080 = (TH1F *)fileResoEP->Get("hReso080");
  Float_t ftcReso[numCent];
  cout << "Reso name: " << fileResoName << endl;

  TString SBDT = "";
  TString SinputFile = "OutputAnalysis/Output_" + inputFileName + "_" + ParticleName[ChosenPart] + SEtaSysChoice[EtaSysChoice] + SBDT;
  if (isApplyWeights)
    SinputFile += "_Weighted";
  if (v2type == 1)
    SinputFile += "_SP";
  if (!useCommonBDTValue)
    SinputFile += "_BDTCentDep";
  if (isRun2Binning)
    SinputFile += "_Run2Binning";
  SinputFile += ".root";
  cout << "Input file: " << SinputFile << endl;
  TFile *inputFile = new TFile(SinputFile);
  TH3D *effWeight3D[numCent];
  TH2D *effWeight2D[numCent];
  TH2F *effWeightVsPt[numCent];
  TH2D *NchVarHisto[numCent];
  TH2D *NchTimesv2Histo[numCent];

  Float_t NchVarMax[numCent];
  Int_t CentFT0CMin, CentFT0CMax;
  TProfile *ProfileNchVar[numCent];
  TProfile *ProfileNchTimesV2[numCent];
  TProfile *ProfileEffWeight[numCent];
  TProfile *ProfileEffWeightVsPt[numCent];

  TH1F *NchVsCentrality = new TH1F("NchVsCentrality", "NchVsCentrality", numCent, fCentFT0C);
  TH1F *NchTimesV2VsCentrality = new TH1F("NchTimesV2VsCentrality", "NchTimesV2VsCentrality", numCent, fCentFT0C);
  TH1F *EffWeightVsCentrality = new TH1F("EffWeightVsCentrality", "EffWeightVsCentrality", numCent, fCentFT0C);
  TH1F *EffWeight[numPtBins];

  for (Int_t pt = 0; pt < numPtBins; pt++)
  {
    EffWeight[pt] = new TH1F(Form("EffWeightPt%i", pt), Form("EffWeightPt%i", pt), numCent, fCentFT0C);
  }
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    if (cent == numCent)
    { // 0-80%
      CentFT0CMin = 0;
      CentFT0CMax = 80;
    }
    else
    {
      CentFT0CMin = CentFT0C[cent];
      CentFT0CMax = CentFT0C[cent + 1];
    }
    effWeight2D[cent] = (TH2D *)inputFile->Get(Form("EffWeight_cent%i-%i", CentFT0CMin, CentFT0CMax));
    effWeight3D[cent] = (TH3D *)inputFile->Get(Form("EffWeight3D_cent%i-%i", CentFT0CMin, CentFT0CMax));
    NchVarHisto[cent] = (TH2D *)inputFile->Get(Form("NchVar_cent%i-%i", CentFT0CMin, CentFT0CMax));
    NchTimesv2Histo[cent] = (TH2D *)inputFile->Get(Form("NchTimesV2_cent%i-%i", CentFT0CMin, CentFT0CMax));
    ProfileEffWeight[cent] = (TProfile *)effWeight2D[cent]->ProfileX(Form("ProfileEffWeight_cent%i-%i", CentFT0CMin, CentFT0CMax));
    ProfileNchVar[cent] = (TProfile *)NchVarHisto[cent]->ProfileX(Form("ProfileNchVar_cent%i-%i", CentFT0CMin, CentFT0CMax));
    ProfileNchTimesV2[cent] = (TProfile *)NchTimesv2Histo[cent]->ProfileX(Form("ProfileNchTimesV2_cent%i-%i", CentFT0CMin, CentFT0CMax));
    NchVarMax[cent] = ProfileNchVar[cent]->GetBinContent(1);
    cout << "Centrality: " << CentFT0CMin << "-" << CentFT0CMax << ", NchVarMax: " << NchVarMax[cent] << endl;
    NchVsCentrality->SetBinContent(cent + 1, NchVarMax[cent]);
    NchVsCentrality->SetBinError(cent + 1, ProfileNchVar[cent]->GetBinError(1));
    NchTimesV2VsCentrality->SetBinContent(cent + 1, ProfileNchTimesV2[cent]->GetBinContent(1));
    NchTimesV2VsCentrality->SetBinError(cent + 1, ProfileNchTimesV2[cent]->GetBinError(1));
    EffWeightVsCentrality->SetBinContent(cent + 1, ProfileEffWeight[cent]->GetBinContent(1));
    EffWeightVsCentrality->SetBinError(cent + 1, ProfileEffWeight[cent]->GetBinError(1));
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      effWeight3D[cent]->GetXaxis()->SetRange(effWeight3D[cent]->GetXaxis()->FindBin(PtBins[pt] + 0.0001), effWeight3D[cent]->GetXaxis()->FindBin(PtBins[pt + 1] - 0.0001));
      TH2F *proj = (TH2F *)effWeight3D[cent]->Project3D("yz");
      proj->SetName(Form("EffWeightPt%i_cent%i-%i", pt, CentFT0CMin, CentFT0CMax));
      TProfile *ProfileEffWeightPt = (TProfile *)proj->ProfileY();
      EffWeight[pt]->SetBinContent(cent + 1, ProfileEffWeightPt->GetBinContent(1));
      EffWeight[pt]->SetBinError(cent + 1, ProfileEffWeightPt->GetBinError(1));
    }
    effWeight3D[cent]->GetXaxis()->SetRangeUser(0, 10);
    effWeight3D[cent]->GetYaxis()->SetRange(1, 1);
    effWeightVsPt[cent] = (TH2F *)effWeight3D[cent]->Project3D("zx");
    ProfileEffWeightVsPt[cent] = (TProfile *)effWeightVsPt[cent]->ProfileX(Form("ProfileEffWeightPt_cent%i-%i", CentFT0CMin, CentFT0CMax));
    effWeight3D[cent]->GetYaxis()->SetRange(-1, -1);
  }

  TCanvas *cNchVsCentrality = new TCanvas("cNchVsCentrality", "cNchVsCentrality", 800, 600);
  StyleCanvas(cNchVsCentrality, 0.12, 0.05, 0.02, 0.13);
  StyleHistoYield(NchVsCentrality, 1, 1.3, ColorMult[0], MarkerMult[0], "Centrality (%)", "N_{in-plane} / N_{ch}", "", 1, 1.15, 1.2);
  NchVsCentrality->Draw("E");

  TCanvas *cNchTimesV2VsCentrality = new TCanvas("cNchTimesV2VsCentrality", "cNchTimesV2VsCentrality", 800, 600);
  StyleCanvas(cNchTimesV2VsCentrality, 0.12, 0.05, 0.02, 0.13);
  StyleHistoYield(NchTimesV2VsCentrality, 1, 1.3, ColorMult[0], MarkerMult[0], "Centrality (%)", "N_{in-plane} #times v_{2}", "", 1, 1.15, 1.2);
  NchTimesV2VsCentrality->Draw("E");

  TCanvas *cEffWeightVsCentrality = new TCanvas("cEffWeightVsCentrality", "cEffWeightVsCentrality", 800, 600);
  StyleCanvas(cEffWeightVsCentrality, 0.12, 0.05, 0.02, 0.13);
  StyleHistoYield(EffWeightVsCentrality, 1, 1.2, ColorMult[0], MarkerMult[0], "Centrality (%)", "Efficiency weight (in-plane)", "", 1, 1.15, 1.2);
  EffWeightVsCentrality->Draw("E");

  TCanvas *cEffWeightVsPt = new TCanvas("cEffWeightVsPt", "cEffWeightVsPt", 800, 600);
  StyleCanvas(cEffWeightVsPt, 0.12, 0.05, 0.02, 0.13);
  TLegend *legEffWeight = new TLegend(0.65, 0.55, 0.9, 0.95);
  legEffWeight->SetBorderSize(0);
  legEffWeight->SetFillStyle(0);
  legEffWeight->SetTextSize(0.03);
  legEffWeight->SetTextFont(42);
  for (Int_t pt = 0; pt < numPtBins; pt++)
  {
    StyleHistoYield(EffWeight[pt], 1, 1.2, ColorMult[pt], MarkerMult[pt], "Centrality (%)", "Efficiency weight (in-plane)", "", 1, 1.15, 1.2);
    legEffWeight->AddEntry(EffWeight[pt], Form("%.1f < p_{T} < %.1f GeV/c", PtBins[pt], PtBins[pt + 1]), "p");
    EffWeight[pt]->Draw("same E");
  }
  legEffWeight->Draw("same");

  TH1F *InducedV2[numPtBins];
  TCanvas *cInducedV2VsPtCent = new TCanvas("cInducedV2VsPtCent", "cInducedV2VsPtCent", 800, 600);
  StyleCanvas(cInducedV2VsPtCent, 0.12, 0.05, 0.02, 0.13);
  for (Int_t pt = 0; pt < numPtBins; pt++)
  {
    InducedV2[pt] = new TH1F(Form("InducedV2Pt%i", pt), Form("InducedV2Pt%i", pt), numCent, fCentFT0C);
    for (Int_t cent = 0; cent < numCent; cent++)
    {
      InducedV2[pt]->SetBinContent(cent + 1, (EffWeight[pt]->GetBinContent(cent + 1) - 1) / 2);
    }
    StyleHistoYield(InducedV2[pt], 0, 0.1, ColorMult[pt], MarkerMult[pt], "Centrality (%)", "Induced v_{2}", "", 1, 1.15, 1.2);
    InducedV2[pt]->Draw("same");
  }
  legEffWeight->Draw("same");

  TH1F *v2Corr[numCent];
  TCanvas *cV2CorrVsCent = new TCanvas("cV2CorrVsCent", "cV2CorrVsCent", 800, 600);
  StyleCanvas(cV2CorrVsCent, 0.12, 0.05, 0.02, 0.13);
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    v2Corr[cent] = new TH1F(Form("v2CorrCent%i", cent), Form("v2CorrCent%i", cent), numPtBins, PtBins);
    for (Int_t pt = 0; pt < numPtBins; pt++)
    {
      v2Corr[cent]->SetBinContent(pt + 1, InducedV2[pt]->GetBinContent(cent + 1));
    }
    StyleHistoYield(v2Corr[cent], 0, 0.1, ColorMult[cent], MarkerMult[cent], "p_{T} (GeV/c)", "correction to v_{2}", "", 1, 1.15, 1.2);
    if (cent == numCent)
      ftcReso[cent] = hReso080->GetBinContent(1);
    else
      ftcReso[cent] = hReso->GetBinContent(hReso->FindBin(CentFT0C[cent] + 0.001));
    cout << "Resolution: " << ftcReso[cent] << endl;
    v2Corr[cent]->Scale(1. / ftcReso[cent]);
    v2Corr[cent]->Draw("same HIST");
  }

  TH1F *V2InducedByEfficiency[numPtBins];
  TCanvas *cV2InducedByEfficiencyVsPtCent = new TCanvas("cV2InducedByEfficiencyVsPtCent", "cV2InducedByEfficiencyVsPtCent", 800, 600);
  StyleCanvas(cV2InducedByEfficiencyVsPtCent, 0.12, 0.05, 0.02, 0.13);
  for (Int_t pt = 0; pt < numPtBins; pt++)
  {
    V2InducedByEfficiency[pt] = new TH1F(Form("V2InducedByEfficiencyPt%i", pt), Form("V2InducedByEfficiencyPt%i", pt), numCent, fCentFT0C);
    for (Int_t cent = 0; cent < numCent; cent++)
    {
      V2InducedByEfficiency[pt]->SetBinContent(cent + 1, (1. / EffWeight[pt]->GetBinContent(cent + 1) - 1) / 2);
    }
    StyleHistoYield(V2InducedByEfficiency[pt], -0.1, 0, ColorMult[pt], MarkerMult[pt], "Centrality (%)", "v_{2} induced by efficiency", "", 1, 1.15, 1.2);
    V2InducedByEfficiency[pt]->Draw("same");
  }
  legEffWeight->Draw("same");

  TCanvas *cEffWeightVsPtCent = new TCanvas("cEffWeightVsPtCent", "cEffWeightVsPtCent", 800, 600);
  StyleCanvas(cEffWeightVsPtCent, 0.12, 0.05, 0.02, 0.13);
  TLegend *legEffWeightCent = new TLegend(0.75, 0.62, 0.9, 0.9);
  legEffWeightCent->SetBorderSize(0);
  legEffWeightCent->SetFillStyle(0);
  legEffWeightCent->SetTextSize(0.03);
  legEffWeightCent->SetTextFont(42);
  for (Int_t cent = 0; cent < numCent; cent++)
  {
    StyleHistoYield((TH1F *)ProfileEffWeightVsPt[cent], 1, 1.2, ColorMult[cent], MarkerMult[cent], "p_{T} (GeV/c)", "Efficiency weight (in-plane)", "", 1, 1.15, 1.2);
    legEffWeightCent->AddEntry(ProfileEffWeightVsPt[cent], Form("%i-%i %%", CentFT0C[cent], CentFT0C[cent + 1]), "p");
    ProfileEffWeightVsPt[cent]->GetXaxis()->SetRangeUser(0.8, 8);
    ProfileEffWeightVsPt[cent]->Draw("same E");
  }
  legEffWeightCent->Draw("same");

  TFile *fileout = new TFile("V2Corr.root", "RECREATE");
  for (Int_t cent = 0; cent < numCent; cent++)
    v2Corr[cent]->Write();
  fileout->Close();

  cout << "Output file: V2Corr.root" << endl;
}
