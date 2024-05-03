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
    TString inputFileName = SinputFileName,
    Bool_t UseTwoGauss = ExtrUseTwoGauss)
{

  TString Sdef = "OutputAnalysis/FitV2_" + inputFileName + "_" + ParticleName[!isXi];
  Sdef += IsOneOrTwoGauss[UseTwoGauss];
  Sdef += SIsBkgParab[BkgType];
  Sdef += Form("_Cent%i-%i", CentFT0C[mul], CentFT0C[mul + 1]);

  TString Svaried = "";

  // plot histos
  TFile *fdef = TFile::Open(Sdef + ".root");
  TH1F *hDefault = (TH1F *)fdef->Get("histoV2");
  hDefault->SetName(Form("hDefault"));

  TCanvas *cv2 = new TCanvas("cv2", "cv2", 1000, 1200);
  cv2->cd();
  hDefault->Draw();

  TH1F *h[trials];
  Float_t BDTscoreCut = 0;
  for (int i = 0; i < trials; i++)
  {
    if (i > 1)
      continue;
    BDTscoreCut = LowerlimitBDTscoreCut + (UpperlimitBDTscoreCut - LowerlimitBDTscoreCut) * 1. / trials * i;
    TString SBDT = Form("_BDT%.3f", BDTscoreCut);
    Svaried = Sdef + SBDT;
    cout << "InputFile - variation: " << Svaried << endl;

    TFile *fvaried = TFile::Open(Svaried + ".root");
    TH1F *hVariedCut = (TH1F *)fvaried->Get("histoV2");
    hVariedCut->SetName(Form("histoV2_%i", i));
    //cv2->cd();
    //hVariedCut->Draw("same");

    h[i] = makeSystPlots(i + 1, Sdef + ".root", Svaried + ".root");
  }

  cv2->SaveAs("Systematics/V2MultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".root");

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
    if (i > 1)
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

  TH1F *hSystMultiTrial = (TH1F *)h[0]->Clone("hSystMultiTrial");
  hSystMultiTrial->Reset();
  for (int pt = 0; pt < bins; pt++)
  {
    hSystMultiTrial->SetBinContent(pt + 1, fgaus[pt]->GetParameter(2));
    hSystMultiTrial->SetBinError(pt + 1, 0);
  }

  TString OutputFile = "Systematics/SystMultiTrial_" + inputFileName + "_" + ParticleName[!isXi] + ".root";
  TFile *Write = new TFile(OutputFile, "RECREATE");
  for (int pt = 0; pt < bins; pt++)
  {
    hPtDev[pt]->Write();
    fgaus[pt]->Write();
  }
  hSystMultiTrial->Write();

  TLatex *ltx = new TLatex();
  ltx->SetTextSize(0.05);
  ltx->SetTextColor(kRed);

  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1200);
  c2->Divide(3, 4);
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
  c2->SaveAs("Systematics/MultTrial" + inputFileName + "_" + ParticleName[!isXi] + ".root");

}
