#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>
#include <string>

// loop over AO2D directories and merge the trees into a single one

void MergeTrees(Int_t isRedFactor = 1, const std::string folderName = "TreeForTrainingBkg", const std::string inputFileName = "AnalysisResultsTree_Bkg_LHC23_PbPb_pass3_Train207099", bool isMC = 1)
{
    TString inputFileNameNew = folderName + "/" + inputFileName + ".root";
    cout <<"Input file name: " << inputFileNameNew << endl;
    TFile inputFile(inputFileNameNew);
    TDirectory *dir = gDirectory;
    TList *keys = dir->GetListOfKeys();
    TList *outlist = new TList;
    std::string treeName = (isMC) ? "O2casctraining" : "O2cascanalysis";
    for (int i = 0; i < keys->GetEntries(); i++)
    {
        TKey *key = (TKey *)keys->At(i);
        TString keyName = key->GetName();
        key->Print();
        if (i%isRedFactor != 0) continue;
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

    TString outputFileNameNew = folderName + "/" + inputFileName + "_New.root";
    if (isRedFactor!=1) outputFileNameNew = folderName + "/" + inputFileName+ "_New_RedFactor" + std::to_string(isRedFactor) + ".root";
    cout <<"Output file name: " << outputFileNameNew << endl;
    TFile outputFile(outputFileNameNew, "RECREATE");
    if(outlist->GetSize() == 0){
        std::cerr << "Error: nothing was found" << std::endl;
        return;
    }
    TTree *newtree = TTree::MergeTrees(outlist);
    newtree->SetName(treeName.c_str());
    newtree->Write();
}