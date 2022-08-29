#!/usr/bin/env python

from __future__ import print_function, division

import ROOT
from my_unfolder import unpickle_region

for algo in ["ak4","ak8"]:
 for var in ["jet_puppiMultiplicity", "jet_pTD","jet_LHA", "jet_width", "jet_thrust",
             "jet_puppiMultiplicity_charged", "jet_pTD_charged", "jet_LHA_charged", "jet_width_charged", "jet_thrust_charged"]:
  for region in ["ZPlusJets", "Dijet_central", "Dijet_forward"]:
    outfile="unfolding_info_"+algo+"_"+region+"_"+var+".root"
    print("CREATE", outfile)

    f="/nfs/dust/cms/user/hinzmann/qganalysis/CMSSW_10_2_17/src/trees/workdir_102X_v3data_v2mc_"+algo+"puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6_PUDownX3/unfolding_output/unfolding_regularizeNone_DATA_subFakes_densityModeBinWidth_constraintNone_dijet_experimentalSystAsAltResponse_scaleSyst_pdfSyst_zpj_ExperimentalSystFromFile_scalSystFromFile_pdfSystFromFile_mergeLastPtBin_zeroBadResponseBins_onlyZPJmulti_signalRegionOnly/"+region+"/"+var+"/unfolding_result.pkl"
    region = unpickle_region(f)

    o=ROOT.TFile(outfile,"RECREATE")

    #print(region.keys())
    #print(region["unfolder"])
    #print(dir(region["unfolder"]))
    print(region["unfolder"].get_rhoij_total())
    print(region["unfolder"].binning_handler)
    region["unfolder"].get_ematrix_input().Write()
    region["unfolder"].get_ematrix_stat_response().Write()
    region["unfolder"].get_ematrix_stat().Write()
    region["unfolder"].get_ematrix_tunfold_total().Write()
    region["unfolder"].get_rhoij_total().Write()
    region["unfolder"].get_probability_matrix().Write()
    region["unfolder"].get_unfolded_with_ematrix_stat().Write()
    region["unfolder"].get_unfolded_with_ematrix_response().Write()
    region["unfolder"].get_folded_unfolded().Write()
    region["unfolder"].get_folded_mc_truth().Write()
    region["unfolder"].generator_binning.Write()
    region["unfolder"].detector_binning.Write()
