#include <Riostream.h>
#include "TProfile2D.h"
#include <string>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TCutG.h>
#include "TFitResult.h"
#include "TLegend.h"

void BulkViscosityParam()
{

    const Int_t npoints = 58;
    Double_t TValues[npoints];
    Double_t BulkValues[npoints];

    ifstream file("../BulkViscosityParamIII.txt");

    if (!file)
    {
        cout << "Error opening file." << endl;
        return;
    }

    int count = 0;

    while (file >> TValues[count] >> BulkValues[count])
    {
        count++;
    }

    file.close();
    for (Int_t i = 0; i < npoints; i++)
    {
        //cout << "T: " << TValues[i] << ", Bulk: " << BulkValues[i] << endl;
    }

    TGraph *gBulkParamIII = new TGraph(npoints);
    for (Int_t i = 0; i < gBulkParamIII->GetN(); i++)
    {
        gBulkParamIII->SetPoint(i, TValues[i], BulkValues[i]);
    }

    const Int_t npointsB = 74;
    Double_t TValuesBayes[npointsB];
    Double_t BulkValuesBayes[npointsB];

    ifstream fileB("../BulkViscosityBayes2019.txt");

    if (!fileB)
    {
        cout << "Error opening file." << endl;
        return;
    }

    int countB = 0;

    while (fileB >> TValuesBayes[countB] >> BulkValuesBayes[countB])
    {
        countB++;
    }

    fileB.close();
    for (Int_t i = 0; i < npointsB; i++)
    {
        //cout << "T: " << TValuesBayes[i] << ", Bulk: " << BulkValuesBayes[i] << endl;
    }

    TGraph *gBulkParamBayes = new TGraph(npointsB);
    for (Int_t i = 0; i < gBulkParamBayes->GetN(); i++)
    {
        gBulkParamBayes->SetPoint(i, TValuesBayes[i], BulkValuesBayes[i]);
    }

    const Int_t npointsUpperCL = 30;
    Double_t TValuesBayesUpperCL[npointsUpperCL];
    Double_t BulkValuesBayesUpperCL[npointsUpperCL];

    ifstream fileBUpperCL("../BulkViscosityBayes2019UpperCL.txt");

    if (!fileBUpperCL)
    {
        cout << "Error opening file." << endl;
        return;
    }

    int countUpperCL = 0;

    while (fileBUpperCL >> TValuesBayesUpperCL[countUpperCL] >> BulkValuesBayesUpperCL[countUpperCL])
    {
        countUpperCL++;
    }

    fileBUpperCL.close();
    for (Int_t i = 0; i < npointsUpperCL; i++)
    {
        //cout << "T: " << TValuesBayesUpperCL[i] << ", Bulk: " << BulkValuesBayesUpperCL[i] << endl;
    }

    TGraph *gBulkParamBayesUpperCL = new TGraph(npointsUpperCL);
    for (Int_t i = 0; i < gBulkParamBayesUpperCL->GetN(); i++)
    {
        gBulkParamBayesUpperCL->SetPoint(i, TValuesBayesUpperCL[i], BulkValuesBayesUpperCL[i]);
    }

    const Int_t npointsLowerCL = 30;
    Double_t TValuesBayesLowerCL[npointsLowerCL];
    Double_t BulkValuesBayesLowerCL[npointsLowerCL];

    ifstream fileBLowerCL("../BulkViscosityBayes2019LowerCL.txt");

    if (!fileBLowerCL)
    {
        cout << "Error opening file." << endl;
        return;
    }

    int countLowerCL = 0;

    while (fileBLowerCL >> TValuesBayesLowerCL[countLowerCL] >> BulkValuesBayesLowerCL[countLowerCL])
    {
        countLowerCL++;
    }

    fileBLowerCL.close();
    for (Int_t i = 0; i < npointsLowerCL; i++)
    {
        //cout << "T: " << TValuesBayesLowerCL[i] << ", Bulk: " << BulkValuesBayesLowerCL[i] << endl;
    }

    TGraph *gBulkParamBayesLowerCL = new TGraph(npointsLowerCL);
    for (Int_t i = 0; i < gBulkParamBayesLowerCL->GetN(); i++)
    {
        gBulkParamBayesLowerCL->SetPoint(i, TValuesBayesLowerCL[i], BulkValuesBayesLowerCL[i]);
    }

    TCanvas *canvasBulkParam = new TCanvas("canvasBulkParam", "canvasBulkParam", 800, 600);
    gBulkParamIII->SetLineColor(kBlue);
    gBulkParamIII->SetLineWidth(2);
    gBulkParamIII->SetTitle("");
    gBulkParamIII->GetXaxis()->SetTitle("Temperature (GeV)");
    gBulkParamIII->GetYaxis()->SetTitle("Bulk Viscosity / Entropy Density");
    gBulkParamIII->Draw("AL");
    gBulkParamBayes->SetLineColor(kRed);
    gBulkParamBayes->SetLineWidth(2);
    gBulkParamBayes->Draw("L SAME");
    gBulkParamBayesUpperCL->SetLineColor(kRed + 2);
    gBulkParamBayesUpperCL->SetLineWidth(2);
    gBulkParamBayesUpperCL->Draw("L SAME");
    gBulkParamBayesLowerCL->SetLineColor(kRed + 2);
    gBulkParamBayesLowerCL->SetLineWidth(2);
    gBulkParamBayesLowerCL->Draw("L SAME");
}