from __future__ import print_function
import hepdata_lib
from hepdata_lib import Submission
from hepdata_lib.helpers import *

def round_value_to_decimals(cont, decimals=3):
    """
    round all values in a dictionary to some decimals in one go
    default round to 3 digits after period
    possible use case: correlations where typical values are within -1,1
    : param cont : dictionary as returned e.g. by RootFileReader::read_hist_1d()
    : type  cont : dictionary
    : param decimals: how many decimals for the rounding
    : type  decimals: integer
    """

    decimals = int(decimals)

    for i, val in enumerate(cont):
        if isinstance(val, tuple):
            cont[i] = (round(val[0], decimals), round(val[1], decimals))
        else:
            cont[i] = round(val, decimals)

def round_value_and_uncertainty(cont_val, cont_unc, sig_digits_unc=2):
    """
    round values and uncertainty according to the precision of the uncertainty,
    and also round uncertainty to a given number of significant digits
    Typical usage:
         reader = RootFileReader("rootfile.root")
         data = reader.read_hist_1d("histogramName")
         round_value_and_uncertainty(data,"y","dy",2)
    will round data["y"] to match the precision of data["dy"] for each element, after
    rounding each element of data["dy"] to 2 significant digits
    e.g. 26.5345 +/- 1.3456 --> 26.5 +/- 1.3
    : param cont : dictionary as returned e.g. by RootFileReader::read_hist_1d()
    : type  cont : dictionary
    : param sig_digits_unc: how many significant digits used to round the uncertainty
    : type  sig_digits_unc: integer
    """

    sig_digits_unc = int(sig_digits_unc)

    for i, (val, unc) in enumerate(zip(cont_val, cont_unc)):
        if isinstance(unc, tuple):
            # case for TGraphAsymmErrors with unc = (elow,ehigh), the central value is rounded
            # using the significant digits of the largest of the two uncertainties,
            # the smaller uncertainty would be rounded accordingly (at least 1 digit)
            # usually lower and higher uncertainties will be of the same order of magnitude
            # or at most different by 1 order (like +0.132  -0.083), in which case,
            # if choosing 2 significant digits, the rounding should result in +0.13  -0.08
            max_absunc = 0.0
            index_min_unc = 0
            # set default precision for both sides of uncertainty
            sig_digits_unc_ntuple = [sig_digits_unc, sig_digits_unc]
            if abs(unc[0]) < abs(unc[1]):
                max_absunc = abs(unc[1])
                index_min_unc = 0
                relative_precision = get_value_precision_wrt_reference(unc[0], unc[1])
            else:
                max_absunc = abs(unc[0])
                index_min_unc = 1
                relative_precision = get_value_precision_wrt_reference(unc[1], unc[0])
            # update precision on smaller uncertainty (at least 1 significant digit)
            sig_digits_unc_ntuple[index_min_unc] = int(max(1, sig_digits_unc + relative_precision))
            cont_unc[i] = (relative_round(unc[0], sig_digits_unc_ntuple[0]),
                                relative_round(unc[1], sig_digits_unc_ntuple[1]))
            val_precision = get_value_precision_wrt_reference(val, max_absunc)
            sig_digits_val = int(sig_digits_unc + val_precision)
            cont_val[i] = relative_round(val, sig_digits_val)
        else:
            # standard case for TH1 or TGraphErrors, uncertainty is a single value
            cont_unc[i] = relative_round(unc, sig_digits_unc)
            val_precision = get_value_precision_wrt_reference(val, unc)
            sig_digits_val = int(sig_digits_unc + val_precision)
            cont_key[i] = relative_round(val, sig_digits_val)

submission = Submission()
submission.read_abstract("abstract.txt")
submission.add_link("Webpage with all figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-20-010/")
submission.add_link("arXiv", "http://arxiv.org/abs/arXiv:2109.03340")
submission.add_record_id(1920187, "inspire")

paper_plots={}
paper_plots["multiplicity","ungroomed","120 < PT < 150 GeV","Z+jet","AK4",None]="Figure 6 (upper left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_puppiMultiplicity/ZPlusJets/ak4/unfolded_ZPlusJets_jet_puppiMultiplicity_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["multiplicity","ungroomed","120 < PT < 150 GeV","central dijet","AK4",None]="Figure 6 (upper right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_puppiMultiplicity/Dijet_central/ak4/unfolded_Dijet_central_jet_puppiMultiplicity_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["pTD2","ungroomed","120 < PT < 150 GeV","Z+jet","AK4",None]="Figure 6 (lower left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_pTD/ZPlusJets/ak4/unfolded_ZPlusJets_jet_pTD_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["pTD2","ungroomed","120 < PT < 150 GeV","central dijet","AK4",None]="Figure 6 (lower right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_pTD/Dijet_central/ak4/unfolded_Dijet_central_jet_pTD_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["thrust","ungroomed","120 < PT < 150 GeV","Z+jet","AK4",None]="Figure 7 (upper left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_thrust/ZPlusJets/ak4/unfolded_ZPlusJets_jet_thrust_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["thrust","ungroomed","120 < PT < 150 GeV","central dijet","AK4",None]="Figure 7 (upper right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_thrust/Dijet_central/ak4/unfolded_Dijet_central_jet_thrust_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["width","ungroomed","120 < PT < 150 GeV","Z+jet","AK4",None]="Figure 7 (lower left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_width/ZPlusJets/ak4/unfolded_ZPlusJets_jet_width_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["width","ungroomed","120 < PT < 150 GeV","central dijet","AK4",None]="Figure 7 (lower right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_width/Dijet_central/ak4/unfolded_Dijet_central_jet_width_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["LHA","ungroomed","120 < PT < 150 GeV","Z+jet","AK4",None]="Figure 8 (left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA/ZPlusJets/ak4/unfolded_ZPlusJets_jet_LHA_alt_truth_bin_3_divBinWidth_paper_ak4.pdf"
paper_plots["LHA","ungroomed","120 < PT < 150 GeV","central dijet","AK4",None]="Figure 8 (right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA/Dijet_central/ak4/unfolded_Dijet_central_jet_LHA_alt_truth_bin_3_divBinWidth_paper_ak4.pdf"
paper_plots["LHA","ungroomed","408 < PT < 1500 GeV","Z+jet","AK4",None]="Figure 9 (upper left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA/ZPlusJets/ak4/unfolded_ZPlusJets_jet_LHA_alt_truth_bin_8_divBinWidth_paper.pdf"
paper_plots["LHA","ungroomed","1000 < PT < 4000 GeV","central dijet","AK4",None]="Figure 9 (upper right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA/Dijet_central/ak4/unfolded_Dijet_central_jet_LHA_alt_truth_bin_12_divBinWidth_paper.pdf"
paper_plots["LHA","ungroomed","120 < PT < 150 GeV","Z+jet","AK8",None]="Figure 9 (lower left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA/ZPlusJets/ak8/unfolded_ZPlusJets_jet_LHA_alt_truth_bin_3_divBinWidth_paper_ak8.pdf"
paper_plots["LHA","ungroomed","120 < PT < 150 GeV","central dijet","AK8",None]="Figure 9 (lower right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA/Dijet_central/ak8/unfolded_Dijet_central_jet_LHA_alt_truth_bin_3_divBinWidth_paper_ak8.pdf"
paper_plots["LHA (charged-only)","ungroomed","120 < PT < 150 GeV","Z+jet","AK4",None]="Figure 10 (upper left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA_charged/ZPlusJets/ak4/unfolded_ZPlusJets_jet_LHA_charged_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["LHA (charged-only)","ungroomed","120 < PT < 150 GeV","central dijet","AK4",None]="Figure 10 (upper right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA_charged/Dijet_central/ak4/unfolded_Dijet_central_jet_LHA_charged_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["LHA","groomed","120 < PT < 150 GeV","Z+jet","AK4",None]="Figure 10 (lower left)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA/ZPlusJets/ak4/unfolded_ZPlusJets_groomed_jet_LHA_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["LHA","groomed","120 < PT < 150 GeV","central dijet","AK4",None]="Figure 10 (lower right)","SMP-20-010/Figures/unfolding_plots_data_unreg_differential/jet_LHA/Dijet_central/ak4/unfolded_Dijet_central_groomed_jet_LHA_alt_truth_bin_3_divBinWidth_paper.pdf"
paper_plots["LHA","ungroomed","vsPt","Z+jet","AK4",None]="Figure 11 (left)","SMP-20-010/Figures/summary_plots/plot_zpj_means_vs_pt_all/zpj_mean_vs_pt_jet_LHA_ak4puppi_onlyDataNomMC_paper.pdf"
paper_plots["LHA","ungroomed","vsPt","central dijet","AK4",None]="Figure 11 (right)","SMP-20-010/Figures/summary_plots/plot_dijet_cen_means_vs_pt_all/dijet_cen_mean_vs_pt_jet_LHA_ak4puppi_onlyDataNomMC_paper.pdf"
#print(paper_plots)

plots=[]
plot=[None,None,None,None,None,None,None,None]
for f in [open("CMS_2018_PAS_SMP_18_QGX_ZPJ.yoda"),open("CMS_2018_PAS_SMP_18_QGX_DIJET.yoda")]:
  i=0
  for l in f.readlines():
    if "Multiplicity" in l: plot[0]="multiplicity"; plot[1]="ungroomed"
    if "p_{T}^{D}" in l: plot[0]="pTD2"; plot[1]="ungroomed"
    if "Thrust" in l: plot[0]="thrust"; plot[1]="ungroomed"
    if "Width" in l: plot[0]="width"; plot[1]="ungroomed"
    if "LHA" in l: plot[0]="LHA"; plot[1]="ungroomed"
    if "Groomed" in l: plot[1]="groomed"
    if "charged" in l: plot[0]+=" (charged-only)"
    if " p_{T}" in l: plot[2]=l.split("<")[0].split("$")[-1].strip()+" < PT < "+l.split("<")[-1].split("$")[0].strip()+" GeV"
    if "Z+jet" in l: plot[3]="Z+jet"
    if "Central dijet" in l: plot[3]="central dijet"
    if "Forward dijet" in l: plot[3]="forward dijet"
    if "AK4" in l: plot[4]="AK4"
    if "AK8" in l: plot[4]="AK8"
    if "# xval" in l: plot[5]=i+1
    if "END" in l:
      plot[6]=i-plot[5]
      plots+=[plot]
      plot=[None,None,None,None,None,None,None,None]
      #print(plots[-1])
    i+=1

#plots=plots[:2]
plot=[None,None,"vsPt",None,None,None,None,None]
for f in [open("plotsvspt.txt")]:
  i=0
  reading=False
  for l in f.readlines():
    if ("NAME" in l or "..." in l or "h_min" in l) and reading:
      reading=False
      plot[6]=i-plot[5]
      plots+=[plot]
      plot=[None,None,"vsPt",None,None,None,None,None]
      #print(plots[-1])
    if "NAME" in l and "puppiMultiplicity" in l: plot[0]="multiplicity"; plot[1]="ungroomed"
    if "NAME" in l and "pTD" in l: plot[0]="pTD2"; plot[1]="ungroomed"
    if "NAME" in l and "thrust" in l: plot[0]="thrust"; plot[1]="ungroomed"
    if "NAME" in l and "width" in l: plot[0]="width"; plot[1]="ungroomed"
    if "NAME" in l and "LHA" in l: plot[0]="LHA"; plot[1]="ungroomed"
    if "NAME" in l and "groomed" in l: plot[1]="groomed"
    if "NAME" in l and "charged" in l: plot[0]+=" (charged-only)"
    if "NAME" in l and "ZPlusJets" in l: plot[3]="Z+jet"
    if "NAME" in l and "Dijet_central" in l: plot[3]="central dijet"
    if "NAME" in l and "Dijet_forward" in l: plot[3]="forward dijet"
    if "NAME" in l and "AK4" in l: plot[4]="AK4"
    if "NAME" in l and "AK8" in l: plot[4]="AK8"
    if "NAME" in l:
      reading=True
      plot[5]=i+1
    i+=1
#print(plots)

#plots=plots[:4]
plot=[None,None,None,None,None,None,None,"Correlation"]
for f in [open("correlations.txt")]:
  i=0
  reading=False
  for l in f.readlines():
    if ("...doing uncert fraction" in l or "Correlation" in l) and reading:
      reading=False
      plot[6]=i-plot[5]
      plots+=[plot]
      plot=plot[:]
      #print(plots[-1])
    if "Algo/Region/Var" in l and "puppiMultiplicity" in l: plot[0]="multiplicity"; plot[1]="ungroomed"
    if "Algo/Region/Var" in l and "pTD" in l: plot[0]="pTD2"; plot[1]="ungroomed"
    if "Algo/Region/Var" in l and "thrust" in l: plot[0]="thrust"; plot[1]="ungroomed"
    if "Algo/Region/Var" in l and "width" in l: plot[0]="width"; plot[1]="ungroomed"
    if "Algo/Region/Var" in l and "LHA" in l: plot[0]="LHA"; plot[1]="ungroomed"
    if "Algo/Region/Var" in l and "groomed" in l: plot[1]="groomed"
    if "Algo/Region/Var" in l and "charged" in l: plot[0]+=" (charged-only)"
    if "Algo/Region/Var" in l and "ZPlusJets" in l: plot[3]="Z+jet"
    if "Algo/Region/Var" in l and "Dijet_central" in l: plot[3]="central dijet"
    if "Algo/Region/Var" in l and "Dijet_forward" in l: plot[3]="forward dijet"
    if "Algo/Region/Var" in l and "AK4" in l: plot[4]="AK4"
    if "Algo/Region/Var" in l and "AK8" in l: plot[4]="AK8"
    if "Correlation" in l:
      reading=True
      plot[2]=l.split(",")[1].strip().replace(".0","")
      plot[5]=i+1
    i+=1
print(plots)

ptbins=[50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000]
table_index={}
table_index2={}
tables_paper=[]
tables_additional=[]
already_done=[]
for plot in plots:
    from hepdata_lib import Table
    pl=(plot[0],plot[1],plot[2],plot[3],plot[4],plot[7])
    if pl in already_done: continue
    else: already_done+=[pl]
    #print(pl)
    if plot[2]=="vsPt":
      additional_name=plot[1][0]+plot[4][-1]+("c" if "charged" in plot[0] else "")+plot[0][0]+"["+plot[3][0]+"] (vs.PT)"
    else:
      additional_name=plot[1][0]+plot[4][-1]+("c" if "charged" in plot[0] else "")+plot[0][0]+"["+plot[3][0]+plot[2].split(" ")[0]+"]"
    if plot[7]=="Correlation":
      additional_name+=" (corr)"
    name,pdf=(paper_plots[pl][0]+" "+additional_name,paper_plots[pl][1]) if (pl in paper_plots.keys() and plot[7]!="Correlation") else ("Extra Figure "+additional_name,"")
    #if pdf=="": continue
    print(name)
    table = Table(name)
    if plot[2]!="vsPt" and plot[7]!="Correlation":
      table_index[name]=1*(0*(plot[3]=="Z+jet")+1*(plot[3]=="central dijet")+2*(plot[3]=="forward dijet"))+\
                        10*int(plot[1][0]=="g")+\
                        100*int(plot[4][-1]=="8")+\
                        1000*int(5*("charged" in plot[0])+0*("multiplicity" in plot[0])+1*("pTD2" in plot[0])+2*("LHA" in plot[0])+3*("width" in plot[0])+4*("thrust" in plot[0]))+\
                        10000*ptbins.index(int(plot[2].split(" ")[0]))
      table_index2[name]=pl
      #print(table_index[name])
    if plot[7]=="Correlation":
      table.description = "Correlation matrix of the particle-level distributions of "+plot[1]+" "+plot[4]+" "+plot[0]+" in "+plot[2]+" in the "+plot[3]+" region."
    elif plot[2]=="vsPt":
      table.description = "Mean of "+plot[1]+" "+plot[0]+" for "+plot[4]+" jets as a function of PT in the "+plot[3]+" region."
    else:
      table.description = "Particle-level distributions of "+plot[1]+" "+plot[4]+" "+plot[0]+" in "+plot[2]+" in the "+plot[3]+" region."
    table.description += " For easier navigation in HEPData, figures are named with a code with each letter/number representing in the following order: ungroomed/groomed, AK4/AK8, charged, multiplicity/pTD2/thrust/width/LHA, [central dijet/forward dijet/Z+jet, PT]."
    table.location = "Data from "+name
    table.keywords["observables"] = ["1/SIG DSIG/DLAMBDA"]
    import numpy as np
    if plot[2]=="vsPt":
      data = np.loadtxt("plotsvspt.txt", skiprows=plot[5], max_rows=plot[6])
    elif plot[7]=="Correlation":
      databins = np.loadtxt("correlations.txt", skiprows=plot[5], max_rows=int((plot[6]-2)/2))
      data = np.loadtxt("correlations.txt", skiprows=plot[5]+int((plot[6]-2)/2), max_rows=int((plot[6]-2)/2))
    elif plot[3]=="Z+jet":
      data = np.loadtxt("CMS_2018_PAS_SMP_18_QGX_ZPJ.yoda", skiprows=plot[5], max_rows=plot[6])
    else:
      data = np.loadtxt("CMS_2018_PAS_SMP_18_QGX_DIJET.yoda", skiprows=plot[5], max_rows=plot[6])
    if plot[3]=="Z+jet":
      table.keywords["reactions"] = ["P P --> Z0 JET"]
    else:
      table.keywords["reactions"] = ["P P --> JET JET"]
    if plot[7]!="Correlation":
      data=np.array([d for d in data if d[3]!=0])
    #print(data)
    from hepdata_lib import Variable
    if plot[7]=="Correlation":
      x = Variable("LAMBDA", is_independent=True, is_binned=True, units="")
      x2 = Variable("LAMBDA", is_independent=True, is_binned=True, units="")
    elif plot[2]=="vsPt":
      x = Variable("PT", is_independent=True, is_binned=True, units="GeV")
    else:
      x = Variable(plot[1].replace("ungroomed","").replace("groomed","groomed ")+plot[0], is_independent=True, is_binned=True, units="")
    if plot[7]!="Correlation":
      x.values = np.column_stack((data[:,(0)]-data[:,(1)],data[:,(0)]+data[:,(2)]))
    else:
      bins = [(min(0,databins[0][0]),databins[0][1])]
      bins.extend(databins[1:])
      x.values = np.array([])
      x2.values = np.array([])
      for d in bins:
       for d2 in bins:
        x.values.extend([(float(d[0]),float(d[1]))])
        x2.values.extend([(float(d2[0]),float(d2[1]))])
    if plot[7]=="Correlation":
      y = Variable("Correlation", is_independent=False, is_binned=False, units="")
    elif plot[2]=="vsPt":
      y = Variable("<LAMBDA>", is_independent=False, is_binned=False, units="")
    else:
      y = Variable("1/SIG DSIG/DLAMBDA", is_independent=False, is_binned=False, units="")
    if plot[7]!="Correlation":
      y.values = data[:,3]
    else:
      y.values = np.array([])
      i=0
      for d in data:
       for dd in d:
        y.values.extend([float(dd)])
        i+=1
      round_value_to_decimals(y.values)
      #print(y.values)
    #print(x.values)
    #print(y.values)
    y.add_qualifier("observable", plot[0])
    y.add_qualifier("jet algorithm", plot[1]+" "+plot[4])
    if plot[2]!="vsPt":
      y.add_qualifier("PT", plot[2])
    y.add_qualifier("region", plot[3])
    y.add_qualifier("SQRT(S)", 13, "TeV")
    table.add_variable(x)
    if plot[7]=="Correlation":
      table.add_variable(x2)
    table.add_variable(y)
    if plot[7]!="Correlation":
      from hepdata_lib import Uncertainty
      unc = Uncertainty("STAT+SYS", is_symmetric=False)
      unc.values = np.column_stack((-data[:,(4)],+data[:,(5)]))
      round_value_and_uncertainty(y.values,unc.values)
      y.add_uncertainty(unc)
    if pdf!="":
      table.add_image(pdf)
    if "xtra" in name:
      tables_additional+=[table]
    else:
      tables_paper+=[table]

def get_name(table): return ("a" if "PT" in table.name else "b")+table.name.replace("u4","aa").replace("g4","ab").replace("u8","ba").replace("g8","bb").replace("0","").replace("1","").replace("2","").replace("3","").replace("4","").replace("5","").replace("6","").replace("7","").replace("8","").replace("9","").replace("upper","apper")
def get_num(table): return int(table.name.split()[1])
tables_paper.sort(key=get_name)
tables_paper.sort(key=get_num)
tables_additional.sort(key=get_name)
for table in tables_paper+tables_additional:
    submission.add_table(table)

### Print map from HepData index to plot names
#i=1
#print("std::map<int,int> hepdata_index {", end='')
#for table in tables_paper+tables_additional:
#  if table.name in table_index.keys():
#    print("{"+str(table_index[table.name])+","+str(i)+"},", end='')
#  i+=1
# print("};")

### Print Plot file    
#i=1
#for table in tables_paper+tables_additional:
#  if table.name in table_index2.keys():
#    print("# BEGIN PLOT /CMS_2018_PAS_SMP_18_QGX_ZPJ/d"+("0"+str(i) if i<10 else str(i))+"-x01-y01")
#    print("Title="+table_index2[table.name][4]+" jets, "+table_index2[table.name][3]+" region, "+table_index2[table.name][2].replace("PT","$p_{T}^{\\text{jet}}$")).replace("< $p_{T}^{\text{jet}}$ <","$<p_{T}^{\text{jet}}<$")
#    print("XLabel="+("groomed " if table_index2[table.name][1]=="groomed" else "")+table_index2[table.name][0])
#    print("YLabel=$\\frac{1}{\mathrm d N / \mathrm d p_{T}}\\frac{\mathrm d^{2}N}{\mathrm d p_{T}~\mathrm d \lambda}$")
#    print("LeftMargin=1.7")
#    print("XMin=0.0")
#    print("XMax="+("150." if "multiplicity" in table_index2[table.name][0] else "1.0"))
#    print("LogY=0")
#    print("XTwosidedTicks=1")
#    print("YTwosidedTicks=1")
#    print("NormalizeToIntegral=1")
#    print("RatioPlotSameStyle=1")
#    print("LegendXPos=0.5")
#    print("# END PLOT")
#    print("")
#  i+=1
    
for table in submission.tables:
  table.keywords["cmenergies"] = [13000]
outdir = "output"
submission.create_files(outdir, remove_old=True)
