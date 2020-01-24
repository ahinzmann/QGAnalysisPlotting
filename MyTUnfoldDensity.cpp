// Hacky code to access protected methods
// This converts protected methods/vars into public ones, so that they can be
// used in PyROOT (which can't see protected things)
// 
// IMPORTANT: don't start with a block comment...wont work

class MyTUnfoldDensity : public TUnfoldDensity {
public:
    MyTUnfoldDensity(const TH2 *hist_A,
                     EHistMap histmap,
                     ERegMode regmode = kRegModeCurvature,
                     EConstraint constraint=kEConstraintArea,
                     EDensityMode densityMode=kDensityModeBinWidthAndUser,
                     const TUnfoldBinning *outputBins=0,
                     const TUnfoldBinning *inputBins=0,
                     const char *regularisationDistribution=0,
                     const char *regularisationAxisSteering="*[UOB]"
    ) {
        
        // TDOMParser parser;
        // Int_t error = parser.ParseFile(detector_binning_xml.c_str());
        // if (error) cout<<"error="<<error<<" from TDOMParser" << endl;
        // TXMLDocument const *XMLdocument = parser.GetXMLDocument();
        // detectorBinning = TUnfoldBinningXML::ImportXML(XMLdocument, "detector");
        // if (!detectorBinning) {cout<<"error: can not read detector binning (document empty?)" << endl;}

        // const TUnfoldBinning * node = detectorBinning->FindNode("detectordistribution");
        // variable_name = node->GetDistributionAxisLabel(0);
        // cout << "Setup for " << variable_name << endl;

        // error = parser.ParseFile(generator_binning_xml.c_str());
        // if (error) cout<<"error="<<error<<" from TDOMParser" << endl;
        // XMLdocument = parser.GetXMLDocument();
        // generatorBinning = TUnfoldBinningXML::ImportXML(XMLdocument, "generator");
        // if (!generatorBinning) {cout<<"error: can not read generator binning (document empty?)" << endl;}

        // string this_variable_name = "LHA";
        // if (!detectorBinning) {
        //     detectorBinning = new TUnfoldBinning("detector");

        //     std::vector<double> variable_bin_edges_reco = {0, 0.07, 0.14, 0.18, 0.22, 0.255, 0.29, 0.32, 0.35, 0.385, 0.42, 0.455, 0.49, 0.525, 0.56, 0.6, 0.64, 0.695, 0.75, 0.875, 1};
        //     int nbins_variable_reco = variable_bin_edges_reco.size()-1;

        //     std::vector<double> pt_bin_edges_underflow_reco = {30, 34, 38, 44, 50};
        //     int nbins_pt_underflow_reco = pt_bin_edges_underflow_reco.size()-1;

        //     std::vector<double> pt_bin_edges_reco = {50, 57.5, 65, 76.5, 88, 104, 120, 135, 150, 168, 186, 220, 254, 290, 326, 367, 408, 444.5, 481, 547.5, 614, 707, 800, 900, 1000, 1250, 1500, 1750, 2000, 4250, 6500};
        //     //std::vector<double> pt_bin_edges_reco = {50, 57.5, 65, 76.5, 88, 104, 120, 135, 150, 168, 186, 220, 254, 290, 326, 367, 408, 444.5, 481, 547.5, 614, 707, 800, 3557, 6500};
        //     int nbins_pt_reco = pt_bin_edges_reco.size()-1;

        //     TUnfoldBinning * detector_distribution_underflow = detectorBinning->AddBinning("detectordistribution_underflow");
        //     detector_distribution_underflow->AddAxis(this_variable_name.c_str(), nbins_variable_reco, &variable_bin_edges_reco[0], false, false);
        //     detector_distribution_underflow->AddAxis("pt", nbins_pt_underflow_reco, &pt_bin_edges_underflow_reco[0], false, false);

        //     TUnfoldBinning * detector_distribution = detectorBinning->AddBinning("detectordistribution");
        //     detector_distribution->AddAxis(this_variable_name.c_str(), nbins_variable_reco, &variable_bin_edges_reco[0], false, false);
        //     detector_distribution->AddAxis("pt", nbins_pt_reco, &pt_bin_edges_reco[0], false, false);
        // }

        // if (!generatorBinning) {

        //     generatorBinning = new TUnfoldBinning("generator");
        //     std::vector<double> variable_bin_edges_gen = {0.0, 0.14, 0.22, 0.29, 0.35, 0.42, 0.49, 0.56, 0.64, 0.75, 1.0};
        //     int nbins_variable_gen = variable_bin_edges_gen.size()-1;

        //     std::vector<double> pt_bin_edges_underflow_gen = {30, 38, 50};
        //     int nbins_pt_underflow_gen = pt_bin_edges_underflow_gen.size()-1;

        //     std::vector<double> pt_bin_edges_gen = {50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 1500, 2000, 6500};
        //     //std::vector<double> pt_bin_edges_gen = {50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 6500};
        //     int nbins_pt_gen = pt_bin_edges_gen.size()-1;

        //     TUnfoldBinning * generator_distribution_underflow = generatorBinning->AddBinning("generatordistribution_underflow");
        //     generator_distribution_underflow->AddAxis(this_variable_name.c_str(), nbins_variable_gen, &variable_bin_edges_gen[0], false, false);
        //     generator_distribution_underflow->AddAxis("pt", nbins_pt_underflow_gen, &pt_bin_edges_underflow_gen[0], false, false);

        //     TUnfoldBinning * generator_distribution = generatorBinning->AddBinning("generatordistribution");
        //     generator_distribution->AddAxis(this_variable_name.c_str(), nbins_variable_gen, &variable_bin_edges_gen[0], false, false);
        //     generator_distribution->AddAxis("pt", nbins_pt_gen, &pt_bin_edges_gen[0], false, false);
        // }
        
        TUnfoldDensity(hist_A,
                       histmap,
                       regmode,
                       constraint,
                       densityMode,
                       outputBins,
                       inputBins,
                       regularisationDistribution,
                       regularisationAxisSteering);
       std::cout << "Setup (My)TUnfoldDensity" << std::endl;
    }

    // TUnfoldBinning *generatorBinning, *detectorBinning;
    // TString variable_name;

    // Here are methods that were protected in TUnfoldDensity, that are now public for us to use

    /// vector of the unfolding result
    inline const TMatrixD *GetX(void) const { return TUnfoldDensity::GetX(); }
    /// covariance matrix of the result
    inline const TMatrixDSparse *GetVxx(void) const { return TUnfoldDensity::GetVxx(); }
    /// inverse of covariance matrix of the result
    inline const TMatrixDSparse *GetVxxInv(void) const { return TUnfoldDensity::GetVxxInv(); }
    /// vector of folded-back result
    inline const TMatrixDSparse *GetAx(void) const { return TUnfoldDensity::GetAx(); }
    /// matrix of derivatives dx/dy
    inline const TMatrixDSparse *GetDXDY(void) const { return TUnfoldDensity::GetDXDY(); }
    /// matrix contributions of the derivative dx/dA
    inline const TMatrixDSparse *GetDXDAM(int i) const { return TUnfoldDensity::GetDXDAM(i); }
    /// vector contributions of the derivative dx/dA
    inline const TMatrixDSparse *GetDXDAZ(int i) const { return TUnfoldDensity::GetDXDAZ(i); }
    /// matrix E<sup>-1</sup>, using internal bin counting
    inline const TMatrixDSparse *GetEinv(void) const { return TUnfoldDensity::GetEinv(); }
    /// matrix E, using internal bin counting
    inline const TMatrixDSparse *GetE(void) const { return TUnfoldDensity::GetE(); }
    /// inverse of covariance matrix of the data y
    inline const TMatrixDSparse *GetVyyInv(void) const { return TUnfoldDensity::GetVyyInv(); }
    /// response matrix A
    inline const TMatrixDSparse *GetA(void) const { return fA; }
    /// mapping array to hist
    inline const TArrayI & GetXToHist(void) const { return fXToHist; }
    /// mapping hist to array
    inline const TArrayI & GetHistToX(void) const { return fHistToX; }

    Bool_t AddRegularisationCondition(Int_t i0, Double_t f0, Int_t i1=-1, Double_t f1=0., Int_t i2=-1, Double_t f2=0.) { return TUnfold::AddRegularisationCondition(i0,f0,i1,f1,i2,f2); } // add regularisation condition for a triplet of output bins
    Bool_t AddRegularisationCondition(Int_t nEle, const Int_t *indices, const Double_t *rowData) { return TUnfold::AddRegularisationCondition(nEle, indices, rowData); } // add a regularisation condition

};