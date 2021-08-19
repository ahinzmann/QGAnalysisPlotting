import ROOT


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()

ROOT.TGaxis.SetExponentOffset(-0.075, 0.01, "y")

My_Style = ROOT.TStyle("myStyle", "Blah")

# For the canvas:
My_Style.SetCanvasBorderMode(0)
My_Style.SetCanvasColor(ROOT.kWhite)
My_Style.SetCanvasDefH(600) #Height of canvas
My_Style.SetCanvasDefW(800) #Width of canvas
My_Style.SetCanvasDefX(0)   #POsition on screen
My_Style.SetCanvasDefY(0)

# For the Pad:
My_Style.SetPadBorderMode(0)
# My_Style.SetPadBorderSize(Width_t size = 1)
My_Style.SetPadColor(ROOT.kWhite)
My_Style.SetPadGridX(False)
My_Style.SetPadGridY(False)
My_Style.SetGridColor(0)
My_Style.SetGridStyle(3)
My_Style.SetGridWidth(1)

# For the frame:
My_Style.SetFrameBorderMode(0)
My_Style.SetFrameBorderSize(1)
My_Style.SetFrameFillColor(0)
My_Style.SetFrameFillStyle(0)
My_Style.SetFrameLineColor(1)
My_Style.SetFrameLineStyle(1)
My_Style.SetFrameLineWidth(1)

# Margins:
# My_Style.SetPadTopMargin(0.05)
My_Style.SetPadBottomMargin(0.13)
My_Style.SetPadLeftMargin(0.15)
My_Style.SetPadRightMargin(0.04)

# For the Global title:
# My_Style.SetOptTitle(0)
My_Style.SetTitleFont(42)
My_Style.SetTitleColor(1)
My_Style.SetTitleTextColor(1)
My_Style.SetTitleFillColor(10)
My_Style.SetTitleFontSize(0.045)
# My_Style.SetTitleH(0) # Set the height of the title box
# My_Style.SetTitleW(0) # Set the width of the title box
My_Style.SetTitleX(0.5) # Set the position of the title box
My_Style.SetTitleAlign(23) # Set the position of the title box
# My_Style.SetTitleY(0.985) # Set the position of the title box
# My_Style.SetTitleStyle(Style_t style = 1001)
My_Style.SetTitleBorderSize(0)


# For the axis titles:
My_Style.SetTitleColor(1, "XYZ")
My_Style.SetTitleFont(42, "XYZ")
My_Style.SetTitleSize(0.038, "XYZ")
# My_Style.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
# My_Style.SetTitleYSize(Float_t size = 0.02)
My_Style.SetTitleXOffset(1.2)
My_Style.SetTitleYOffset(1.5)

# For the statistics box:
My_Style.SetOptFile(0)  # what does this even do?
My_Style.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
My_Style.SetOptFit(1111)  # doesn't work with THStack unfortunately

# For the legend
My_Style.SetLegendBorderSize(0)

My_Style.SetPaintTextFormat(".3g")

My_Style.SetHistTopMargin(0.)

My_Style.SetNumberContours(512)

final_reading_factor=1.2
My_Style.SetLegendTextSize(0.03*final_reading_factor)
