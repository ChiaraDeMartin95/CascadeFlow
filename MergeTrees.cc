#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>
#include <string>

// loop over AO2D directories and merge the trees into a single one

void MergeTrees(Int_t ChosenPart = 0, Int_t isRedFactor = 1, const std::string inputFileName = "AnalysisResults_trees_LHC25_OO_Train487953", bool isMC = 0, bool isEff = 0)
{
    std::string folderName = "TreeForAnalysis";
    if (isEff)
        folderName = "FileForEfficiency";
    TString inputFileNameNew = folderName + "/" + inputFileName + ".root";
    cout << "Input file name: " << inputFileNameNew << endl;
    TFile inputFile(inputFileNameNew);
    TDirectory *dir = gDirectory;
    TList *keys = dir->GetListOfKeys();
    TList *outlist = new TList;
    std::string treeName = (isMC) ? "O2casctraining" : "O2cascanalysis";
    if (isEff)
        treeName = "O2cascanalysis";
    if (ChosenPart == 6)
        treeName = "O2lambdaanalysis";
    for (int i = 0; i < keys->GetEntries(); i++)
    {
        TKey *key = (TKey *)keys->At(i);
        TString keyName = key->GetName();
        key->Print();
        if (i % isRedFactor != 0)
            continue;
        cout << "Index i: " << i << endl;
        if (key->IsFolder() && keyName.BeginsWith("DF"))
        {
            TDirectory *dir = (TDirectory *)inputFile.Get(keyName);
            TTree *tree = (TTree *)dir->Get(treeName.c_str());

            if (!tree)
            {
                std::cerr << "Error: could not get TTree from directory " << keyName << std::endl;
                continue;
            }
            else
            {
                // tree->SetDirectory(0);
                outlist->Add(tree);
                std::cout << "Tree Added!" << std::endl;
            }
        }
    }

    bool isOneDF = 0;
    TString outputFileNameNew = folderName + "/" + inputFileName + "_New.root";
    if (isRedFactor != 1)
        outputFileNameNew = folderName + "/" + inputFileName + "_New_RedFactor" + std::to_string(isRedFactor) + ".root";
    cout << "Output file name: " << outputFileNameNew << endl;
    TFile outputFile(outputFileNameNew, "RECREATE");
    if (outlist->GetSize() == 0)
    {
        std::cerr << "Error: nothing was found" << std::endl;
        return;
    }
    if (outlist->GetSize() == 1)
    {
        isOneDF = 1;
        std::cerr << "1 DF found was found" << std::endl;
    }

    TTree *newtree = TTree::MergeTrees(outlist);
    // if (isOneDF)
    //    newtree = (TTree *)outlist->At(0);
    newtree->SetName(treeName.c_str());
    newtree->Write();
    outputFile.Close();
}