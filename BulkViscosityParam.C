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

    //from paper https://arxiv.org/abs/2404.14295
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
    TGraph *gBulkParamIII = new TGraph(npoints);
    for (Int_t i = 0; i < gBulkParamIII->GetN(); i++)
    {
        gBulkParamIII->SetPoint(i, TValues[i], BulkValues[i]);
    }


    //from paper https://www.nature.com/articles/s41567-019-0611-8
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
    TGraph *gBulkParamBayesLowerCL = new TGraph(npointsLowerCL);
    for (Int_t i = 0; i < gBulkParamBayesLowerCL->GetN(); i++)
    {
        gBulkParamBayesLowerCL->SetPoint(i, TValuesBayesLowerCL[i], BulkValuesBayesLowerCL[i]);
    }

    //from paper https://arxiv.org/pdf/2605.29383
    const Int_t npointsBec = 25;
    Double_t TValuesBayesBec[npointsBec];
    Double_t BulkValuesBayesBec[npointsBec];

    ifstream fileBec("../BulkViscosityBecattiniNoPz.txt");

    if (!fileBec)
    {
        cout << "Error opening file." << endl;
        return;
    }

    int countBec = 0;

    while (fileBec >> TValuesBayesBec[countBec] >> BulkValuesBayesBec[countBec])
    {
        countBec++;
    }

    fileBec.close();

    TGraph *gBulkParamBayesBec = new TGraph(npointsBec);
    for (Int_t i = 0; i < gBulkParamBayesBec->GetN(); i++)
    {
        gBulkParamBayesBec->SetPoint(i, TValuesBayesBec[i], BulkValuesBayesBec[i]);
    }

    const Int_t npointsBecUpperCL = 13;
    Double_t TValuesBayesBecUpperCL[npointsBecUpperCL];
    Double_t BulkValuesBayesBecUpperCL[npointsBecUpperCL];

    ifstream fileBecUpperCL("../BulkViscosityBecattiniNoPzUpperCL.txt");

    if (!fileBecUpperCL)
    {
        cout << "Error opening file." << endl;
        return;
    }

    int countBecUpperCL = 0;

    while (fileBecUpperCL >> TValuesBayesBecUpperCL[countBecUpperCL] >> BulkValuesBayesBecUpperCL[countBecUpperCL])
    {
        countBecUpperCL++;
    }

    fileBecUpperCL.close();

    TGraph *gBulkParamBayesBecUpperCL = new TGraph(npointsBecUpperCL);
    for (Int_t i = 0; i < gBulkParamBayesBecUpperCL->GetN(); i++)
    {
        gBulkParamBayesBecUpperCL->SetPoint(i, TValuesBayesBecUpperCL[i], BulkValuesBayesBecUpperCL[i]);
    }

    const Int_t npointsBecLowerCL = 22;
    Double_t TValuesBayesBecLowerCL[npointsBecLowerCL];
    Double_t BulkValuesBayesBecLowerCL[npointsBecLowerCL];

    ifstream fileBecLowerCL("../BulkViscosityBecattiniNoPzLowerCL.txt");
    if (!fileBecLowerCL)
    {
        cout << "Error opening file." << endl;
        return;
    }
    int countBecLowerCL = 0;
    while (fileBecLowerCL >> TValuesBayesBecLowerCL[countBecLowerCL] >> BulkValuesBayesBecLowerCL[countBecLowerCL])
    {
        countBecLowerCL++;
    }
    fileBecLowerCL.close();
    TGraph *gBulkParamBayesBecLowerCL = new TGraph(npointsBecLowerCL);
    for (Int_t i = 0; i < gBulkParamBayesBecLowerCL->GetN(); i++)
    {
        gBulkParamBayesBecLowerCL->SetPoint(i, TValuesBayesBecLowerCL[i], BulkValuesBayesBecLowerCL[i]);
    }

    //from paper https://arxiv.org/pdf/2010.15130
    const Int_t npointsGovert = 25;
    Double_t TValuesBayesGovert[npointsGovert];
    Double_t BulkValuesBayesGovert[npointsGovert];

    ifstream fileGovert("../BulkViscosityGovert.txt");

    if (!fileGovert)
    {
        cout << "Error opening file." << endl;
        return;
    }

    int countGovert = 0;

    while (fileGovert >> TValuesBayesGovert[countGovert] >> BulkValuesBayesGovert[countGovert])
    {
        countGovert++;
    }

    fileGovert.close();
    TGraph *gBulkParamBayesGovert = new TGraph(npointsGovert);
    for (Int_t i = 0; i < gBulkParamBayesGovert->GetN(); i++)
    {
        gBulkParamBayesGovert->SetPoint(i, TValuesBayesGovert[i], BulkValuesBayesGovert[i]);
    }

    const Int_t npointsGovertUpperCL = 14;
    Double_t TValuesBayesGovertUpperCL[npointsGovertUpperCL];
    Double_t BulkValuesBayesGovertUpperCL[npointsGovertUpperCL];

    ifstream fileGovertUpperCL("../BulkViscosityGovertUpperCL.txt");
    if (!fileGovertUpperCL)
    {
        cout << "Error opening file." << endl;
        return;
    }
    int countGovertUpperCL = 0;
    while (fileGovertUpperCL >> TValuesBayesGovertUpperCL[countGovertUpperCL] >> BulkValuesBayesGovertUpperCL[countGovertUpperCL])
    {
        countGovertUpperCL++;
    }
    fileGovertUpperCL.close();
    TGraph *gBulkParamBayesGovertUpperCL = new TGraph(npointsGovertUpperCL);
    for (Int_t i = 0; i < gBulkParamBayesGovertUpperCL->GetN(); i++)
    {
        gBulkParamBayesGovertUpperCL->SetPoint(i, TValuesBayesGovertUpperCL[i], BulkValuesBayesGovertUpperCL[i]);
    }


    TCanvas *canvasBulkParam = new TCanvas("canvasBulkParam", "canvasBulkParam", 800, 600);
    gBulkParamIII->SetLineColor(kBlue);
    gBulkParamIII->SetLineWidth(2);
    gBulkParamIII->SetTitle("");
    gBulkParamIII->GetXaxis()->SetTitle("Temperature (GeV)");
    gBulkParamIII->GetYaxis()->SetTitle("Bulk Viscosity / Entropy Density");
    gBulkParamIII->Draw("AL SAME");
    gBulkParamBayes->SetLineColor(kRed);
    gBulkParamBayes->SetLineWidth(2);
    gBulkParamBayes->Draw("L SAME");
    gBulkParamBayesUpperCL->SetLineColor(kRed + 2);
    gBulkParamBayesUpperCL->SetLineWidth(2);
    gBulkParamBayesUpperCL->Draw("L SAME");
    gBulkParamBayesLowerCL->SetLineColor(kRed + 2);
    gBulkParamBayesLowerCL->SetLineWidth(2);
    gBulkParamBayesLowerCL->Draw("L SAME");
    gBulkParamBayesBec->SetLineColor(kViolet + 1);
    gBulkParamBayesBec->SetLineWidth(2);
    gBulkParamBayesBec->Draw("L SAME");
    gBulkParamBayesBecUpperCL->SetLineColor(kViolet);
    gBulkParamBayesBecUpperCL->SetLineWidth(2);
    gBulkParamBayesBecUpperCL->Draw("L SAME");
    gBulkParamBayesBecLowerCL->SetLineColor(kViolet);
    gBulkParamBayesBecLowerCL->SetLineWidth(2);
    gBulkParamBayesBecLowerCL->Draw("L SAME");
    gBulkParamBayesGovert->SetLineColor(kGreen + 2);
    gBulkParamBayesGovert->SetLineWidth(2);
    gBulkParamBayesGovert->Draw("L SAME");
    gBulkParamBayesGovertUpperCL->SetLineColor(kGreen + 1);
    gBulkParamBayesGovertUpperCL->SetLineWidth(2);
    gBulkParamBayesGovertUpperCL->Draw("L SAME");

    TLegend *legendBulkParam = new TLegend(0.52, 0.7, 0.87, 0.9);
    legendBulkParam->SetHeader("Bulk Viscosity");
    legendBulkParam->SetNColumns(1);
    legendBulkParam->SetFillStyle(0);
    legendBulkParam->SetTextSize(0.03);
    legendBulkParam->SetBorderSize(0);
    legendBulkParam->AddEntry(gBulkParamIII, "Param III (Palermo)", "l");
    legendBulkParam->AddEntry(gBulkParamBayes, "Bayes 2019 (90\% CL)", "l");
    legendBulkParam->AddEntry(gBulkParamBayesBec, "Becattini 2605.29383 (90\% CL)", "l");
    legendBulkParam->AddEntry(gBulkParamBayesGovert, "Govert 2010.15130 (90\% CL)", "l");
    legendBulkParam->Draw();
}