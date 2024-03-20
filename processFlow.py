import argparse
import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml import plot_utils
from hipe4ml_converter.h4ml_converter import H4MLConverter

parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--isOmega', action='store_true', help="True for Omegas", default=False)
parser.add_argument('--isIntegratedPt', action='store_true', help="True for pt integrated training", default=False)
args = parser.parse_args()
print("is Omega: ", args.isOmega)
isXi = not args.isOmega
print("isIntegratedPt: ", args.isIntegratedPt)

FileBkg="TreeForTrainingBkg/AnalysisResultsTree_Bkg_Train180896.root"#170556.root"
FileSig="TreeForTrainingSignal/AnalysisResultsTree_Signal_Train181283.root"#176512.root"
bkgCandidates= TreeHandler(FileBkg,'O2casctraining', folder_name='DF_*')
sigCandidates= TreeHandler(FileSig,'O2casctraining', folder_name='DF_*')

preselection_string_common = 'abs(fBachBaryonDCAxyToPV) < 10 and fCascRadius < 33 and fV0Radius < 35 and abs(fDCABachToPV) < 20 and fV0CosPA > 0.9'

preselection_string_sig = 'abs(fMcPdgCode) == 3312'
preselection_string_mass = 'fMassXi > 1.28 and fMassXi < 1.36'
preselection_sidebands = 'fMassXi > 1.305 and fMassXi < 1.31 or fMassXi > 1.332 and fMassXi < 1.338'
if not isXi:
    preselection_string_sig = 'abs(fMcPdgCode) == 3334'
    preselection_string_mass = 'fMassOmega > 1.63 and fMassOmega < 1.73'
    preselection_sidebands = 'fMassOmega > 1.656 and fMassOmega < 1.662 or fMassOmega > 1.683 and fMassOmega < 1.689'
    
#competing mass rejection
sigCandidates.apply_preselections(preselection_string_sig)
sigCandidates.apply_preselections(preselection_string_common)
bkgCandidates.apply_preselections(preselection_string_common)
bkgCandidates.apply_preselections(preselection_string_mass)
sigCandidates.apply_preselections(preselection_string_mass)
bkgCandidates.apply_preselections(preselection_sidebands)

vars_to_draw = sigCandidates.get_var_names()
features_for_train = vars_to_draw.copy()
#features_for_train = ['fCpa', 'fCpaXY', 'fDecayLength', 'fDecayLengthXY', 'fDeltaMassPhi', 'fImpactParameterXY', 'fAbsCos3PiK',
#                     'fMaxNormalisedDeltaIP'] 
features_for_train.remove('fMassXi')
features_for_train.remove('fMassOmega')
features_for_train.remove('fPt')
features_for_train.remove('fMcPdgCode')
features_for_train.remove('fSign')
features_for_train.remove('fEta')
features_for_train.remove('fMassLambdaDau')
features_for_train.remove('fMultFT0M')

print('Number of signal candidates: ', sigCandidates.get_n_cand())
print('Number of background candidates: ', bkgCandidates.get_n_cand())
print('Ratio of signal to background: ', sigCandidates.get_n_cand()/bkgCandidates.get_n_cand())
if (10*sigCandidates.get_n_cand()/bkgCandidates.get_n_cand() < 1): 
    bkgCandidatesRed = bkgCandidates.get_subset('', 10*sigCandidates.get_n_cand()/bkgCandidates.get_n_cand())
else:
    bkgCandidatesRed = bkgCandidates.get_subset('', 1)
print('Number of signal candidates (after): ', sigCandidates.get_n_cand())
print('Number of background candidates (after): ', bkgCandidatesRed.get_n_cand())

leg_labels = ['background', 'signal']
plot_utils.plot_distr([bkgCandidatesRed, sigCandidates], vars_to_draw, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig("TrainingPlots/Distributions.png")

vars_to_draw_mass = ['fMassXi']
if not isXi: 
    vars_to_draw_mass = ['fMassOmega']
plot_utils.plot_distr([bkgCandidatesRed, sigCandidates], vars_to_draw_mass, bins=100, labels=leg_labels, log=False, density=True, figsize=(12, 7), alpha=0.3, grid=False)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig("TrainingPlots/MassXi.png")

vars_to_draw_bis = [
                'fCascRadius', 'fV0Radius', 'fCascCosPA', 
                'fV0CosPA', 'fDCAPosToPV', 'fDCANegToPV', 
                'fDCABachToPV', 'fDCACascDaughters', 'fDCAV0Daughters', 
                'fDCAV0ToPV','fBachBaryonCosPA', 'fBachBaryonDCAxyToPV']
#plot_utils.plot_distr([bkgCandidatesRed, sigCandidates], vars_to_draw_bis, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
plot_utils.plot_distr([bkgCandidatesRed, sigCandidates], features_for_train, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig("TrainingPlots/DistributionsInputOnly.png")

plot_utils.plot_corr([bkgCandidatesRed, sigCandidates], vars_to_draw, leg_labels)
plt.savefig("TrainingPlots/Correlations.png")
#plt.show()

npt = 3
minpt = 0.6
if not isXi: 
    minpt = 0.8
#ptbin = [minpt, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0]
ptbin = [minpt, 1.0, 2.0, 3.0, 10.0]
x = slice(1, npt+2)
ptbinMax = ptbin[x] 
nsig = [0, 0, 0, 0, 0, 0, 0]
nbkg = [0, 0, 0, 0, 0, 0, 0]

if (args.isIntegratedPt):
    ptbin[0] = minpt
    ptbinMax[0] = 10
print('ptbinMax: ', ptbinMax)

if isXi:
    Cascade_string = 'Xi'
else:
    Cascade_string = 'Omega'

for ptbin, ptbinMax, nsig, nbkg in zip(ptbin, ptbinMax, nsig, nbkg):
    if ptbin == 10: break
    if (args.isIntegratedPt and ptbin != minpt): break
    print('pt bin: ', ptbin)
    print('pt+1 bin: ', ptbinMax)
    preselection_string_pt = 'fPt > ' + str(ptbin) + ' and fPt < ' + str(ptbinMax)
    sigCandidatesNew = sigCandidates.apply_preselections(preselection_string_pt, False)
    bkgCandidatesPtSel = bkgCandidates.apply_preselections(preselection_string_pt, False)

    print('Number of signal candidates: ', sigCandidatesNew.get_n_cand())
    print('Number of background candidates: ', bkgCandidatesPtSel.get_n_cand())
    print('Ratio of signal to background: ', sigCandidatesNew.get_n_cand()/bkgCandidatesPtSel.get_n_cand())
    if (10*sigCandidatesNew.get_n_cand()/bkgCandidatesPtSel.get_n_cand() < 1): 
        bkgCandidatesNew = bkgCandidatesPtSel.get_subset('', 10*sigCandidatesNew.get_n_cand()/bkgCandidatesPtSel.get_n_cand())
    else: 
        bkgCandidatesNew = bkgCandidatesPtSel.get_subset('', 1)
    print('Number of signal candidates (after): ', sigCandidatesNew.get_n_cand())
    print('Number of background candidates (after): ', bkgCandidatesNew.get_n_cand())

    nsig = sigCandidatesNew.get_n_cand()
    nbkg = bkgCandidatesNew.get_n_cand()
    print('nsig: ', nsig)
    print('nbkg: ', nbkg)
    #here I associate the signal and background candidates to the label 1 and 0, respectively
    train_test_data = train_test_generator([sigCandidatesNew, bkgCandidatesNew], [1,0], test_size=0.5, random_state=42)

    #Draw input features in pt intervals
    plot_utils.plot_distr([bkgCandidatesNew, sigCandidatesNew], vars_to_draw_bis, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig("TrainingPlots/DistributionsInputOnly_" + Cascade_string + str(ptbin)+ "_" + str(ptbinMax)+".png")

    #Draw all features in pt intervals
    #label 0 is for background and 1 for signal
    leg_labels = ['background', 'signal']
    plot_utils.plot_distr([bkgCandidatesNew, sigCandidatesNew], vars_to_draw, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig("TrainingPlots/Distributions_" + Cascade_string + str(ptbin)+ "_" + str(ptbinMax)+".png")
    
    plot_utils.plot_distr([bkgCandidatesNew, sigCandidatesNew], vars_to_draw_mass, bins=100, labels=leg_labels, log=False, density=True, figsize=(12, 7), alpha=0.3, grid=False)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig("TrainingPlots/Mass" + Cascade_string + str(ptbin)+ "_" + str(ptbinMax)+".png")

    model_clf = xgb.XGBClassifier()
    model_hdl = ModelHandler(model_clf, features_for_train)
    
    #optimization of the parameters
    hyper_pars_ranges = {'n_estimators': (200, 1000), 'max_depth': (2, 4), 'learning_rate': (0.01, 0.1)}
    #model_hdl.optimize_params_optuna(train_test_data, hyper_pars_ranges, cross_val_scoring='roc_auc', timeout=120,
     #                                n_jobs=-1, n_trials=100, direction='maximize')
    
    #training of the model
    model_hdl.train_test_model(train_test_data)
    
    y_pred_train = model_hdl.predict(train_test_data[0], False) #0
    y_pred_test = model_hdl.predict(train_test_data[2], False) #2
    
    #distributions of the BDT score
    plt.rcParams["figure.figsize"] = (10, 7)
    ml_out_fig = plot_utils.plot_output_train_test(model_hdl, train_test_data, 100, 
                                                   False, leg_labels, True, density=True)
    
    plt.savefig("TrainingPlots/BDTScore" + Cascade_string + str(ptbin)+ "_" + str(ptbinMax)+".png")
    
    #ROC curve
    roc_train_test_fig = plot_utils.plot_roc_train_test(train_test_data[3], y_pred_test,
                                                        train_test_data[1], y_pred_train, None, leg_labels)
    
    plt.savefig("TrainingPlots/ROC" + Cascade_string + str(ptbin)+ "_" + str(ptbinMax)+".png")
    
    #feature importance
    plot_utils.plot_feature_imp(train_test_data[2], train_test_data[3], model_hdl) 
    plt.savefig("TrainingPlots/FeatureImportance"+ Cascade_string + str(ptbin)+ "_" + str(ptbinMax)+".png")
    
    #dataCandidates.apply_model_handler(model_hdl, False)
    #selected_data_hndl = dataCandidates.get_subset('model_output>0.7')
    ##print('Number of selected candidates: ', selected_data_hndl.get_n_entries())
    #labels_list = ["after selection","before selection"]
    #colors_list = ['orangered', 'cornflowerblue']
    #plot_utils.plot_distr([selected_data_hndl, bkgCandidates], column='fMassXi', bins=200, labels=labels_list, colors=colors_list, density=True,fill=True, histtype='step', alpha=0.5)
    #ax = plt.gca()
    #ax.set_xlabel(r'm (GeV/$c^2$)')
    #ax.margins(x=0)
    #ax.xaxis.set_label_coords(0.9, -0.075)
    ##plt.show()
    #plt.savefig("InvMassBkgXi" + str(pt) +".png")
    
    #store the trained model
    model_hdl.dump_model_handler("ModelHandler/ModelHandler_" + Cascade_string + str(ptbin)+ "_" + str(ptbinMax)+".pickle")
    model_hdl.dump_original_model("ModelHandler/XGBoostModel_" + Cascade_string + str(ptbin)+ "_" + str(ptbinMax)+".pickle")
    
    #convert model in ONNX
    model_converter = H4MLConverter(model_hdl) # create the converter object
    model_onnx = model_converter.convert_model_onnx(1, len(features_for_train))
    model_converter.dump_model_onnx("Onnx/model_onnx"+ Cascade_string + "_" + str(ptbin)+ "_" + str(ptbinMax)+".onnx") # dump the model in ONNX format

    plt.close()

plt.bar(ptbin,nbkg)
plt.savefig("TrainingPlots/Nbkg"+ Cascade_string +".png")
plt.bar(ptbin,nsig)
plt.savefig("TrainingPlots/Nsig"+ Cascade_string +".png")