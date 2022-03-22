/*
  ==============================================================================

    DynamicStiffString.h
    Created: 7 Feb 2022 5:09:58pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include <fstream>
#include "Global.h"
//==============================================================================
/*
*/
class DynamicStiffString  : public juce::Component
{
public:
    DynamicStiffString (NamedValueSet& parameters, double k);
    ~DynamicStiffString() override;

    void paint (juce::Graphics&) override;
    void resized() override;
    
    // function to draw the state of the string
    Path visualiseState (Graphics& g, double visualScaling, double& length);

    void calculateScheme();
    void updateStates();
    
    void refreshParameter (int changedParameterIdx, double changedParameterValue);
    
    void refreshCoefficients (bool init = false);
    
    //return u at the current sample at a location given by the length ratio

    double getOutput (double Lratio)
    {
        return  v[1][static_cast<int> (round(N * Lratio))];
    }
    
    void excite (int loc = -1);
    void mouseDown (const MouseEvent& e) override;
    
    bool shouldExcite() { return excitationFlag; };
    
    void addRemovePoint();
    void refreshCustomIp();
private:
    
    // Model parameters
    double L, rho, r, A, T, E, I, cSq, kappaSq, sigma0, sigma1, lambdaSq, muSq, h, k;
    double origR, origL, origE, origT, origRho;
    std::vector<double*> parameterPtrs; // to easily locate parameters
    std::vector<double> parametersToGoTo;
    std::vector<bool> parameterChanged;

    // Number of intervals (N+1 is number of points including boundaries)
    int N, Nmax, Nprev = 0;
    
    // Number of intervals of subsystems
    int Mv, Mw;
    
    // Fractional number of intervals used for dynaic grid
    double Nfrac, NfracPrev;
    double alf, Iterm, A0, A1, A2, A3, AA;
    
    // (N+1) x 3 'matrices' containing the state of the left and right system at all time-steps
    std::vector<std::vector<double>> vStates;
    std::vector<std::vector<double>> wStates;

    // vectors of pointers that point to state vectors
    std::vector<double*> v;
    std::vector<double*> w;

    /* Scheme variables
        - Adiv for u^{n+1} (that all terms get divided by)
        - B for u^n
        - C for u^{n-1}
        - S for precalculated sigma terms
    */
    double Adiv, B0, Bss, B1, B2, C0, C1, S0, S1;
    
    // flag to tell MainComponent whether to excite the scheme or not
    bool excitationFlag = false;
    
    // initialise location of excitation
    double excitationLoc = 0.5;
    
    bool clamped = true;
    int numFromRightBound = 2;
    
    std::vector<double> customIp;
    
    std::ofstream uSave, MvSave, MwSave, alfSave;
    
#ifdef RECORD
    int counter = 0;
#endif
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (DynamicStiffString)
};
