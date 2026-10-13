// pT spectrum of Lambda daughters of Xi, from Xi fit functions (TF1).
//
// Every TF1 found in the input file (optionally filtered by name) is taken as
// the fully corrected Xi spectrum d2N/(dpT dy) in |y| < 0.5 for one centrality
// class. Xi are generated with pT sampled from the TF1, flat in rapidity in
// |y| < yMaxXi (boost invariance), uniform in phi, and decayed with
// TGenPhaseSpace into Lambda + pi. The Lambda spectrum is filled in |y| < yMaxLambda
// with the PtBins binning from CommonVarLambda.h.
//
// Output, for each TF1 <name>:
//   hLambdaFromXi_<name> : d2N/(dpT dy) of Lambda from Xi, same absolute
//                          normalisation as the TF1 (e.g. per event), incl. BR
//   hXiCheck_<name>      : generated Xi in |y| < yMaxLambda, same binning;
//                          must reproduce the TF1 (closure check)
//   fXi_<name>           : copy of the input TF1
//
// Usage:
//   root -l -b -q 'LambdaFromXi.C+("XiFits.root")'
//   root -l -b -q 'LambdaFromXi.C+("XiFits.root","LambdaFromXi.root",5e6,1.0,0.5,"fLevy")'

#include "TFile.h"
#include "TKey.h"
#include "TF1.h"
#include "TH1D.h"
#include "TString.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include <iostream>
#include <vector>

#include "CommonVarPub.h"
#include "CommonVarLambda.h" // numPtBins, PtBins

TString inputFileNameS = "../XiSpectraSara/fitfunctions.root";
const float yMaxXi = 0.5;
void LambdaFromXi(double yMaxLambda = 10, TString inFileName = inputFileNameS, const char *outFileName = "../LambdaFromXi.root",
                  Long64_t nGen = 5000000,
                  const char *nameFilter = "", double ptMinXi = 0., double ptMaxXi = 15.,
                  bool isXi0 = false)
{
    // PDG masses (GeV/c^2) and branching ratio Xi -> Lambda pi
    const double mXi = isXi0 ? 1.31486 : 1.32171;
    const double mLambda = 1.115683;
    const double mPi = isXi0 ? 0.1349768 : 0.13957039;
    const double BR = isXi0 ? 0.99524 : 0.99887;

    TFile *fin = TFile::Open(inFileName);
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open " << inFileName << "\n";
        return;
    }

    std::vector<TF1 *> funcs;
    TIter nextKey(fin->GetListOfKeys());
    while (TKey *key = (TKey *)nextKey())
    {
        TString name = key->GetName();
        if (nameFilter && nameFilter[0] && !name.Contains(nameFilter))
            continue;
        TF1 *f = dynamic_cast<TF1 *>(key->ReadObj());
        if (f)
            funcs.push_back(f);
    }
    if (funcs.empty())
    {
        std::cerr << "No TF1 found in " << inFileName << "\n";
        return;
    }

    gRandom = new TRandom3(0);
    TGenPhaseSpace gen;
    double masses[2] = {mLambda, mPi};
    TLorentzVector pXi;

    TFile fout(outFileName, "RECREATE");
    std::vector<TH1D *> hLambdas;

    for (TF1 *fXi : funcs)
    {
        TString name = fXi->GetName();
        fXi->SetNpx(5000); // fine grid for GetRandom/Integral
        cout << "f1 eval at 0.5: " << fXi->Eval(0.5) << "\n";

        // Xi per event generated in |y| < yMaxXi
        double dNdy = fXi->Integral(ptMinXi, ptMaxXi);
        double w = dNdy * 2 * yMaxXi / nGen; // weight of one generated Xi (in this way hXiCheck will have the proper normalization)

        TH1D *hLambda = new TH1D("hLambdaFromXi_" + name,
                                 "#Lambda from #Xi, " + name + ";#it{p}_{T} (GeV/#it{c});d^{2}#it{N}/(d#it{p}_{T}d#it{y})",
                                 numPtBins, PtBins);
        TH1D *hXiCheck = (TH1D *)hLambda->Clone("hXiCheck_" + name);
        hXiCheck->SetTitle("Generated #Xi, " + name);
        hLambda->Sumw2();
        hXiCheck->Sumw2();
        Long64_t counter = 0;
        for (Long64_t i = 0; i < nGen; ++i)
        {
            double pt = fXi->GetRandom(ptMinXi, ptMaxXi);
            double y = gRandom->Uniform(-yMaxXi, yMaxXi);
            double phi = gRandom->Uniform(0, TMath::TwoPi());
            double mt = TMath::Sqrt(pt * pt + mXi * mXi);
            pXi.SetPxPyPzE(pt * TMath::Cos(phi), pt * TMath::Sin(phi),
                           mt * TMath::SinH(y), mt * TMath::CosH(y));

            hXiCheck->Fill(pt, w);

            if (!gen.SetDecay(pXi, 2, masses))
            {
                counter++;
                continue;
            }
            double wDecay = gen.Generate(); // constant for 2-body decays (=1)
            TLorentzVector *pL = gen.GetDecay(0);
            if (TMath::Abs(pL->Rapidity()) < yMaxLambda)
                hLambda->Fill(pL->Pt(), w * wDecay * BR);
        }
        cout << "Failed decays: " << counter << "\n";
        
        hLambda->Scale(1., "width");
        // -> d2N/(dpT dy)
        hXiCheck->Scale(1. / (2 * yMaxXi), "width");

        // Closure: generated Xi vs TF1 integral in the same range
        double xiGen = hXiCheck->Integral("width");
        double xiFun = fXi->Integral(PtBins[0], PtBins[numPtBins]);
        std::cout << name << ": dN/dy(Xi) = " << dNdy
                  << ", Xi closure (gen/TF1 in " << PtBins[0] << "-" << PtBins[numPtBins]
                  << " GeV/c) = " << xiGen / xiFun
                  << ", Lambda-from-Xi dN/dy in range = " << hLambda->Integral("width") << "\n";

        fout.cd();
        hLambda->Write();
        hXiCheck->Write();
        fXi->Write("fXi_" + name);
        hLambdas.push_back(hLambda);
        cout << "f1 eval at 0.5 (end): " << fXi->Eval(0.5) << "\n";
    }

    // Overview plot
    TCanvas *c = new TCanvas("cLambdaFromXi", "", 800, 600);
    c->SetLogy();
    TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.88);
    for (size_t i = 0; i < hLambdas.size(); ++i)
    {
        hLambdas[i]->SetLineColor(i + 1 == 5 ? 28 : i + 1); // skip yellow
        hLambdas[i]->SetMarkerColor(hLambdas[i]->GetLineColor());
        hLambdas[i]->SetMarkerStyle(20);
        hLambdas[i]->Draw(i == 0 ? "E" : "E same");
        leg->AddEntry(hLambdas[i], funcs[i]->GetName(), "lp");
    }
    leg->Draw();
    c->Write();
    c->SaveAs(TString(outFileName).ReplaceAll(".root", ".pdf"));

    fout.Close();
    fin->Close();
    std::cout << "Saved " << outFileName << "\n";
}
