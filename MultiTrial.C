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

double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr)
{
  // Error in a Ratio
  if (B != 0)
  {
    Double_t errorfromtop = Aerr * Aerr / (B * B);
    Double_t errorfrombottom = ((A * A) / (B * B * B * B)) * Berr * Berr;
    return TMath::Sqrt(TMath::Abs(errorfromtop - errorfrombottom));
  }
  return 1.;
}

void DivideAndComputeRogerBarlow(TH1F *hVar, TH1F *hDef)
{
  // Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lhVarNBins = hVar->GetNbinsX();
  Double_t lhDefNBins = hDef->GetNbinsX();

  if (lhVarNBins != lhDefNBins)
  {
    cout << "Problem! Number of bins doesn't match! " << endl;
    return;
  }

  Double_t lSigmaDelta[100];
  for (Int_t i = 1; i < hVar->GetNbinsX() + 1; i++)
  {
    // Computation of roger barlow sigma_{delta}
    lSigmaDelta[i] = TMath::Sqrt(TMath::Abs(TMath::Power(hVar->GetBinError(i), 2) - TMath::Power(hDef->GetBinError(i), 2)));
    // Computation of relationship to hDef for plotting in ratio plot
    if (hDef->GetBinContent(i) > 1e-12)
      lSigmaDelta[i] /= hDef->GetBinContent(i);
    else
      lSigmaDelta[i] = 0;
  }
  // Regular Division
  hVar->Divide(hDef);
  // Replace Errors
  for (Int_t i = 1; i < hVar->GetNbinsX() + 1; i++)
  {
    hVar->SetBinError(i, lSigmaDelta[i]);
  }
}

double PassRogerBarlowCriterion(int nsigmas, Double_t dev, Double_t RBsigma)
{

  if (TMath::Abs(dev) > (nsigmas * RBsigma))
  {
    return dev;
  }
  else
  {
    return 0.;
  }
}

TH1F *makeSystPlots(int num = 1, TString Sdef = "", TString Svaried = "")
{

  TFile *fdef = TFile::Open(Sdef);
  TFile *fvaried = TFile::Open(Svaried);

  TH1F *hVariedCut = (TH1F *)fvaried->Get("histoV2");
  hVariedCut->SetName(Form("hVarCutSet%i", num));
  TH1F *hDefault = (TH1F *)fdef->Get("histoV2");
  hDefault->SetName(Form("hDefault"));

  TH1F *hDev = (TH1F *)hDefault->Clone("hDev");
  hDev->Reset();

  DivideAndComputeRogerBarlow(hVariedCut, hDefault);
  // assigns to hVariedCut: hR = hvar/hdef as bin content, sigmaBarlow (sB) as error

  for (int i = 1; i <= hDev->GetNbinsX(); i++)
  {
    double dev = hVariedCut->GetBinContent(i) - 1; // hR - 1
    double err = hVariedCut->GetBinError(i);       // sB

    hDev->SetBinContent(i, PassRogerBarlowCriterion(nsigmaBarlow, dev, err)); // rel. syst. error = hR-1 if |hR-1| > 1*sB
  }

  hDev->GetYaxis()->SetRangeUser(0., 0.5);

  // histo with rel. syst. error (if variation passes the Barlow cut)
  return hDev;
}

void MultiTrial(
    Int_t mul = 0,
    Bool_t isXi = ChosenParticleXi,
    Int_t EtaSysChoice = ExtrEtaSysChoice,
    Int_t BkgType = ExtrBkgType,
    TString inputFileName = SinputFileNameSyst,
    Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  TString Sdef = "OutputAnalysis/FitV2_" + inputFileName + "_" + ParticleName[!isXi];
  Sdef += IsOneOrTwoGauss[UseTwoGauss];
  Sdef += SIsBkgParab[BkgType];
  Sdef += Form("_Cent%i-%i", CentFT0C[mul], CentFT0C[mul + 1]);

  TString Svaried = "";

  TFile *fdef = TFile::Open(Sdef + ".root");

  // v2 plots
  TCanvas *cv2 = new TCanvas("cv2", "cv2", 1000, 800);
  StyleCanvas(cv2, 0.15, 0.05, 0.05, 0.15);
  cv2->cd();
  TH1F *hDefault = (TH1F *)fdef->Get("histoV2");
  hDefault->SetName("hDefault");
  hDefault->SetLineColor(kBlack);
  hDefault->Draw();
  TH1F *hDefaultYield = (TH1F *)fdef->Get("histoYield");
  hDefaultYield->SetName("hDefaultYield");
  TH1F *hDefaultPurity = (TH1F *)fdef->Get("histoPurity");
  hDefaultPurity->SetName("hDefaultPurity");
  TLegend *legTrial = new TLegend(0.66, 0.3, 0.96, 0.6);
  legTrial->SetBorderSize(0);
  legTrial->SetFillStyle(0);
  legTrial->SetTextSize(0.03);
  legTrial->SetTextFont(42);
  legTrial->AddEntry(hDefault, Form("BDT score > %.3f", DefaultBDTscoreCut), "pl");

  TH1F *h[trials];
  TH1F *hRawYield[trials];
  TH1F *hRawYieldRatio[trials];
  TH1F *hPurity[trials];
  TH1F *hPurityRatio[trials];
  Float_t BDTscoreCut = 0;
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trials * i;
    TString SBDT = Form("_BDT%.3f", BDTscoreCut);
    Svaried = Sdef + SBDT;
    cout << "InputFile - variation: " << Svaried << endl;

    TFile *fvaried = TFile::Open(Svaried + ".root");
    TH1F *hVariedCut = (TH1F *)fvaried->Get("histoV2");
    hVariedCut->SetName(Form("histoV2_%i", i));

    hRawYield[i] = (TH1F *)fvaried->Get("histoYield");
    hRawYield[i]->SetName(Form("histoYield_%i", i));
    hRawYieldRatio[i] = (TH1F *)hRawYield[i]->Clone(Form("histoYieldRatio_%i", i));
    hRawYieldRatio[i]->Divide(hDefaultYield);

    hPurity[i] = (TH1F *)fvaried->Get("histoPurity");
    hPurity[i]->SetName(Form("histoPurity_%i", i));
    hPurityRatio[i] = (TH1F *)hPurity[i]->Clone(Form("histoPurity_%i", i));
    hPurityRatio[i]->Divide(hDefaultPurity);

    // v2 plots
    cv2->cd();
    hVariedCut->SetLineColor(ColorMult[i]);
    hVariedCut->SetMarkerColor(ColorMult[i]);
    hVariedCut->SetMarkerStyle(MarkerMult[i]);
    legTrial->AddEntry(hVariedCut, Form("BDT score > %.3f", BDTscoreCut), "pl");
    hVariedCut->Draw("same");

    h[i] = makeSystPlots(i + 1, Sdef + ".root", Svaried + ".root");
  }
  legTrial->Draw();
  cv2->SaveAs("Systematics/V2MultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".pdf");
  cv2->SaveAs("Systematics/V2MultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".png");

  // raw yields plots
  TCanvas *cYield = new TCanvas("cYield", "cYield", 1000, 800);
  StyleCanvas(cYield, 0.15, 0.05, 0.05, 0.15);
  cYield->cd();
  hDefaultYield->SetTitle("");
  hDefaultYield->SetLineColor(kBlack);
  hDefaultYield->GetYaxis()->SetRangeUser(0, 1.5 * hDefaultYield->GetMaximum());
  hDefaultYield->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;
    hRawYield[i]->SetLineColor(ColorMult[i]);
    hRawYield[i]->SetMarkerColor(ColorMult[i]);
    hRawYield[i]->SetMarkerStyle(MarkerMult[i]);
    hRawYield[i]->Draw("same");
  }
  legTrial->Draw();
  cYield->SaveAs("Systematics/YieldMultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".pdf");
  cYield->SaveAs("Systematics/YieldMultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".png");

  // ratio of raw yields to default one
  TCanvas *cYieldRatio = new TCanvas("cYieldRatio", "cYieldRatio", 1000, 800);
  StyleCanvas(cYieldRatio, 0.15, 0.05, 0.05, 0.15);
  cYieldRatio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;
    hRawYieldRatio[i]->SetLineColor(ColorMult[i]);
    hRawYieldRatio[i]->SetMarkerColor(ColorMult[i]);
    hRawYieldRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hRawYieldRatio[i]->SetTitle("");
    hRawYieldRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hRawYieldRatio[i]->GetYaxis()->SetRangeUser(0, 2);
    hRawYieldRatio[i]->Draw("same");
  }
  TF1 *lineat1 = new TF1("lineat1", "1", 0, 5);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineStyle(7);
  legTrial->Draw();
  lineat1->Draw("same");
  cYieldRatio->SaveAs("Systematics/YieldRatioMultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".pdf");
  cYieldRatio->SaveAs("Systematics/YieldRatioMultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".png");

  // purity plots
  TCanvas *cPurity = new TCanvas("cPurity", "cPurity", 1000, 800);
  StyleCanvas(cPurity, 0.15, 0.05, 0.05, 0.15);
  cPurity->cd();
  hDefaultPurity->SetTitle("");
  hDefaultPurity->SetLineColor(kBlack);
  hDefaultPurity->GetYaxis()->SetRangeUser(0, 1);
  hDefaultPurity->Draw();
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;
    hPurity[i]->SetLineColor(ColorMult[i]);
    hPurity[i]->SetMarkerColor(ColorMult[i]);
    hPurity[i]->SetMarkerStyle(MarkerMult[i]);
    hPurity[i]->Draw("same");
  }
  legTrial->Draw();
  cPurity->SaveAs("Systematics/PurityMultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".pdf");
  cPurity->SaveAs("Systematics/PurityMultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".png");

  // ratio of purities to default one
  TCanvas *cPurityRatio = new TCanvas("cPurityRatio", "cPurityRatio", 1000, 800);
  StyleCanvas(cPurityRatio, 0.15, 0.05, 0.05, 0.15);
  cPurityRatio->cd();
  for (int i = 0; i < trials; i++)
  {
    if (i == 8)
      continue;
    hPurityRatio[i]->SetLineColor(ColorMult[i]);
    hPurityRatio[i]->SetMarkerColor(ColorMult[i]);
    hPurityRatio[i]->SetMarkerStyle(MarkerMult[i]);
    hPurityRatio[i]->SetTitle("");
    hPurityRatio[i]->GetYaxis()->SetTitle("Ratio to default");
    hPurityRatio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
    hPurityRatio[i]->Draw("same");
  }
  lineat1->Draw("same");
  legTrial->Draw();
  cPurityRatio->SaveAs("Systematics/PurityRatioMultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".pdf");
  cPurityRatio->SaveAs("Systematics/PurityRatioMultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".png");

  // gaussian distributions + fits
  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1200);
  c2->Divide(3, 4);

  const int bins = h[0]->GetNbinsX();
  TH1F *hPtDev[bins];
  TF1 *fgaus[bins];

  for (int pt = 0; pt < bins; pt++) // loop over pT bins
  {
    hPtDev[pt] = new TH1F(Form("hPtDev%i", pt), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(pt + 1), h[0]->GetBinLowEdge(pt + 2)), 15, -0.5, +0.5);
    hPtDev[pt]->SetStats(0);
    fgaus[pt] = new TF1(Form("fgaus%i", pt), "gaus", -0.5, +0.5);
  }

  for (int i = 1; i < trials; i++)
  {
    if (i == 8)
      continue;
    for (int pt = 0; pt < bins; pt++)
    {
      hPtDev[pt]->Fill(h[i]->GetBinContent(pt + 1));
    }
  }

  for (int pt = 0; pt < bins; pt++)
  {
    fgaus[pt]->SetParameter(0, hPtDev[pt]->GetMaximum());
    hPtDev[pt]->Fit(fgaus[pt], "R");
  }

  TLatex *ltx = new TLatex();
  ltx->SetTextSize(0.05);
  ltx->SetTextColor(kRed);

  for (int i = 2; i < bins; i++)
  {
    c2->cd(i - 2 + 1);
    c2->cd(i - 2 + 1)->SetBottomMargin(0.15);
    hPtDev[i]->GetYaxis()->SetRangeUser(0., hPtDev[i]->GetMaximum() * 1.2);
    hPtDev[i]->GetXaxis()->SetTitleSize(0.06);
    hPtDev[i]->GetXaxis()->SetTitleOffset(1.);
    hPtDev[i]->Draw("EP");
    hPtDev[i]->SetMarkerStyle(kFullCircle);
    hPtDev[i]->Draw("EP SAME");
    fgaus[i]->Draw("same");

    ltx->DrawLatexNDC(0.6, 0.8, Form("#mu = %.3f", fgaus[i]->GetParameter(1)));
    ltx->DrawLatexNDC(0.6, 0.7, Form("#sigma = %.3f", fgaus[i]->GetParameter(2)));
    ltx->DrawLatexNDC(0.6, 0.6, Form("#chi^{2}/ndf = %.1f", fgaus[i]->GetChisquare() / fgaus[i]->GetNDF()));
  }
  c2->SaveAs("Systematics/MultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".pdf");
  c2->SaveAs("Systematics/MultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".png");

  // relative syst. uncertainty
  TCanvas *cSyst = new TCanvas("cSyst", "cSyst", 1000, 800);
  TH1F *hSystMultiTrial = (TH1F *)h[0]->Clone("hSystMultiTrial");
  hSystMultiTrial->Reset();
  for (int pt = 0; pt < bins; pt++)
  {
    hSystMultiTrial->SetBinContent(pt + 1, fgaus[pt]->GetParameter(2));
    hSystMultiTrial->SetBinError(pt + 1, 0);
  }
  hSystMultiTrial->Draw();

  // save histos in output files
  TString OutputFile = "Systematics/SystMultiTrial_" + inputFileName + "_" + ParticleName[!isXi] + ".root";
  TFile *Write = new TFile(OutputFile, "RECREATE");
  for (int pt = 0; pt < bins; pt++)
  {
    hPtDev[pt]->Write();
    fgaus[pt]->Write();
  }
  hSystMultiTrial->Write();
}
