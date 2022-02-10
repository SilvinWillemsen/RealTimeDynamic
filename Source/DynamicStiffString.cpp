/*
  ==============================================================================

    DynamicStiffString.cpp
    Created: 7 Feb 2022 5:09:58pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include <JuceHeader.h>
#include "DynamicStiffString.h"

//==============================================================================
DynamicStiffString::DynamicStiffString (NamedValueSet& parameters, double k) : k (k)
{
    
    // Initialise member variables using the parameter set
    L = *parameters.getVarPointer ("L");
    rho = *parameters.getVarPointer ("rho");
    A = *parameters.getVarPointer ("A");
    T = *parameters.getVarPointer ("T");
    E = *parameters.getVarPointer ("E");
    I = *parameters.getVarPointer ("I");
    sigma0 = *parameters.getVarPointer ("sigma0");
    sigma1 = *parameters.getVarPointer ("sigma1");
    
    parameterPtrs.reserve(8);
    parameterPtrs.push_back (&L);
    parameterPtrs.push_back (&rho);
    parameterPtrs.push_back (&A);
    parameterPtrs.push_back (&T);
    parameterPtrs.push_back (&E);
    parameterPtrs.push_back (&I);
    parameterPtrs.push_back (&sigma0);
    parameterPtrs.push_back (&sigma1);
    
    parametersToGoTo.resize (parameterPtrs.size(), 0);
    for (int i = 0; i < parameterPtrs.size(); ++i)
        parametersToGoTo[i] = *parameterPtrs[i];
    
    double cSqMin = 0.5 * T / (2.0 * rho * 2.0 * A);
    
    double hMin = sqrt(cSqMin) * k;
    int Nmax = floor (2.0 * L / hMin);
    
    // only add to left system (v)
    int MvMax = Nmax - numFromRightBound;
    
    Mw = numFromRightBound;
    // Initialise vectors (excluding outer boundaries
    vStates = std::vector<std::vector<double>> (3,
                                        std::vector<double>(MvMax+1, 0));
    
    wStates = std::vector<std::vector<double>> (3,
                                        std::vector<double>(numFromRightBound+1, 0));

    /*  Make u pointers point to the first index of the state vectors.
        To use u (and obtain a vector from the state vectors) use indices like u[n][l] where,
             - n = 0 is u^{n+1},
             - n = 1 is u^n, and
             - n = 2 is u^{n-1}.
        Also see calculateScheme()
     */
    
    // Initialise pointer vector
    v.resize (3, nullptr);
    w.resize (3, nullptr);

    // Make set memory addresses to first index of the state vectors.
    for (int i = 0; i < 3; ++i)
    {
        v[i] = &vStates[i][0];
        w[i] = &wStates[i][0];
    }
    customIp.resize(4, 0);
    
    refreshCoefficients();
    
    Nprev = N;
    NfracPrev = Nfrac;
    
    uSave.open("uSaveDynamic.csv");
    MvSave.open("MvSave.csv");
    alfSave.open("alfSave.csv");
    excite();
    
}

DynamicStiffString::~DynamicStiffString()
{
}

void DynamicStiffString::paint (juce::Graphics& g)
{
    // clear the background
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
    
    // choose your favourite colour
    g.setColour(Colours::cyan);
    
    // draw the state
    g.strokePath(visualiseState (g, 100), PathStrokeType(2.0f));

}

Path DynamicStiffString::visualiseState (Graphics& g, double visualScaling)
{
    // String-boundaries are in the vertical middle of the component
    double stringBoundaries = getHeight() / 2.0;
    
    // initialise path
    Path stringPath;
    
    // start path
    stringPath.startNewSubPath (0, -v[1][0] * visualScaling + stringBoundaries);
    
    double spacing = static_cast<double>(getWidth()) / Nfrac;
    double x = spacing;
    
    for (int l = 1; l <= Mv; ++l)
    {
        // Needs to be -u, because a positive u would visually go down
        float newY = -v[1][l] * visualScaling + stringBoundaries;
        
        // if we get NAN values, make sure that we don't get an exception
        if (isnan(newY))
            newY = 0;
        
        stringPath.lineTo (x, newY);
        x += spacing;
    }
    x -= spacing;
    x += alf;
    for (int l = 0; l < Mw; ++l)
    {
        float newY = -w[1][l] * visualScaling + stringBoundaries;
        
        // if we get NAN values, make sure that we don't get an exception
        if (isnan(newY))
            newY = 0;
        
        stringPath.lineTo (x, newY);
        x += spacing;

    }

    return stringPath;
}

void DynamicStiffString::resized()
{

}

void DynamicStiffString::calculateScheme()
{
    
    // simply supported left boundary
    v[0][1] = Bss * v[1][1] + B1 * (v[1][2] + v[1][0]) + B2 * v[1][3]
            + C0 * v[2][1] + C1 * (v[2][2] + v[2][0]);
    
    // main left scheme
    for (int l = 2; l < Mv-1; ++l)
        v[0][l] = B0 * v[1][l] + B1 * (v[1][l + 1] + v[1][l - 1]) + B2 * (v[1][l + 2] + v[1][l - 2])
                + C0 * v[2][l] + C1 * (v[2][l + 1] + v[2][l - 1]);
    
    // next-to-boundary point
    v[0][Mv-1] = (2.0 * v[1][Mv-1] - v[2][Mv-1]
        + lambdaSq * (v[1][Mv] - 2.0 * v[1][Mv - 1] + v[1][Mv - 2])
        - muSq * (v[1][Mv - 3] - 4 * v[1][Mv - 2] + 6 * v[1][Mv-1] + (Iterm - 4) * v[1][Mv] + w[1][0])
        + S0 * v[2][Mv-1]
        + S1 * (v[1][Mv] - 2.0 * v[1][Mv - 1] + v[1][Mv - 2])
        - S1 * (v[2][Mv] - 2.0 * v[2][Mv - 1] + v[2][Mv - 2])) / (1.0 + S0);
    
    // boundary points
    v[0][Mv] = (2.0 * v[1][Mv]
        + lambdaSq * (w[1][0] + (Iterm - 2.0) * v[1][Mv] + v[1][Mv - 1])
        - muSq * (v[1][Mv - 2] - 4 * v[1][Mv - 1] + ((Iterm - 2) * (Iterm - 2) + 2) * v[1][Mv] + (2 * Iterm - 4) * w[1][0])
        + (-1.0 + S0) * v[2][Mv]
        + S1 * (w[1][0] + (Iterm - 2.0) * v[1][Mv] + v[1][Mv - 1])
        - S1 * (w[2][0] + (Iterm - 2.0) * v[2][Mv] + v[2][Mv - 1])) / (1.0 + S0);

    // right system (single point now)
    w[0][0] = (2.0 * w[1][0]
        + lambdaSq * (-Iterm * v[1][Mv-1] + v[1][Mv] + (Iterm - 2.0) * w[1][0] + w[1][1]) // w[1][1] is 0
        - muSq * (-Iterm * v[1][Mv-2] + (-Iterm * Iterm + 4 * Iterm + 1) * v[1][Mv-1] + (Iterm - 4) * v[1][Mv]
                  + ((Iterm-2) * (Iterm-2) + 1) * w[1][0] - 4 * w[1][1])  // w[1][1] is 0
        + (-1.0 + S0) * w[2][0]
        + S1 * (-Iterm * v[1][Mv-1] + v[1][Mv] + (Iterm - 2.0) * w[1][0] + w[1][1])
        - S1 * (-Iterm * v[2][Mv-1] + v[2][Mv] + (Iterm - 2.0) * w[2][0] + w[2][1])) / (1.0 + S0);

    for (int i = 0; i <= vStates[0].size(); ++i)
        uSave << v[1][i] << ",";

    uSave << w[1][0] << "," << w[1][1] << "\n;";
    
    
    
    
}

void DynamicStiffString::updateStates()
{
    // Do a pointer-switch. MUCH quicker than copying two entire state vectors every time-step.
    double* vTmp = v[2];
    v[2] = v[1];
    v[1] = v[0];
    v[0] = vTmp;
    
    double* wTmp = w[2];
    w[2] = w[1];
    w[1] = w[0];
    w[0] = wTmp;

    
    NfracPrev = Nfrac;
    Nprev = N;
    
    ++counter;
    if (counter > 500)
    {
        uSave.close();
        MvSave.close();
        alfSave.close();
    }

}

void DynamicStiffString::excite()
{
    //// Arbitrary excitation function (raised cosine) ////
    
    // width (in grid points) of the excitation
    double width = 10;
    
    // make sure we're not going out of bounds at the left boundary
    int start = std::max (floor((N+1) * excitationLoc) - floor(width * 0.5), 1.0);

    for (int l = 0; l < width; ++l)
    {
        // make sure we're not going out of bounds at the right boundary (this does 'cut off' the raised cosine)
        if (l+start > (clamped ? N - 2 : N - 1))
            break;
        
        v[1][l+start] += 0.5 * (1 - cos(2.0 * double_Pi * l / (width-1.0)));
        v[2][l+start] += 0.5 * (1 - cos(2.0 * double_Pi * l / (width-1.0)));
    }
    // Disable the excitation flag to only excite once
    excitationFlag = false;
}

void DynamicStiffString::mouseDown (const MouseEvent& e)
{
    // Get the excitation location as a ratio between the x-location of the mouse-click and the width of the app
    excitationLoc = e.x / static_cast<double> (getWidth());
    
    // Activate the excitation flag to be used by the MainComponent to excite the string
    excitationFlag = true;
}

void DynamicStiffString::refreshParameter (int changedParameterIdx, double changedParameterValue)
{
    parametersToGoTo[changedParameterIdx] = changedParameterValue;
}

void DynamicStiffString::refreshCoefficients()
{
    
    for (int i = 0; i < parameterPtrs.size(); ++i)
        *parameterPtrs[i] = 0.9999 * (*parameterPtrs[i]) + 0.0001 * parametersToGoTo[i];
    
    // Calculate wave speed (squared)
    cSq = T / (rho * A);
    
    // Calculate stiffness coefficient (squared)
    kappaSq = E * I / (rho * A);

    double stabilityTerm = cSq * k * k + 4.0 * sigma1 * k; // just easier to write down below
    
    h =  sqrt (0.5 * (stabilityTerm + sqrt ((stabilityTerm * stabilityTerm) + 16.0 * kappaSq * k * k)));
    Nfrac = L / h;
    
    // check if the change does not surpass a limit
    N = floor (Nfrac);
    alf = Nfrac - N;
    if (Nprev == 0)
        Nprev = N;
    
    // Check whether a grid point needs to be added or removed
    if (Nprev != N)
        addRemovePoint();
    
    Mv = N - numFromRightBound;
    MvSave << Mv << ";\n";
    alfSave << alf << ";\n";

    Iterm = (alf - 1.0) / (alf + 1.0);
    
    lambdaSq = cSq * k * k / (h * h);
    muSq = kappaSq * k * k / (h * h * h * h);
//    std::cout << lambdaSq + 4.0 * muSq + 4.0 * sigma1 * k / (h*h) << std::endl;
    // Coefficients used for damping
    S0 = sigma0 * k;
    S1 = (2.0 * sigma1 * k) / (h * h);
    
    // Scheme coefficients
    B0 = 2.0 - 2.0 * lambdaSq - 6.0 * muSq - 2.0 * S1; // u_l^n
    Bss = 2.0 - 2.0 * lambdaSq - 5.0 * muSq - 2.0 * S1;
    B1 = lambdaSq + 4.0 * muSq + S1;                   // u_{l+-1}^n
    B2 = -muSq;                                        // u_{l+-2}^n
    C0 = -1.0 + S0 + 2.0 * S1;                         // u_l^{n-1}
    C1 = -S1;                                          // u_{l+-1}^{n-1}
    
    Adiv = 1.0 / (1.0 + S0);                           // u_l^{n+1}
    
    // Divide by u_l^{n+1} term
    B0 *= Adiv;
    Bss *= Adiv;
    B1 *= Adiv;
    B2 *= Adiv;
    C0 *= Adiv;
    C1 *= Adiv;
}

void DynamicStiffString::addRemovePoint()
{
    jassert (abs (N-Nprev) <= 1);
    refreshCustomIp();
    if (N > NfracPrev)
    {
        // possibly unnecessary to update up[0]
        v[0][Mv + 1] = customIp[0] * v[0][Mv-1]
            + customIp[1] * v[0][Mv]
            + customIp[2] * w[0][0];
        
        v[1][Mv + 1] = customIp[0] * v[1][Mv-1]
            + customIp[1] * v[1][Mv]
            + customIp[2] * w[1][0];

        v[2][Mv + 1] = customIp[0] * v[2][Mv-1]
           + customIp[1] * v[2][Mv]
            + customIp[2] * w[2][0];
    
    } else {
        v[0][Mv] = 0;
        v[1][Mv] = 0;
        v[2][Mv] = 0;
    }
}

void DynamicStiffString::refreshCustomIp()
{
    customIp[0] = -alf * (alf + 1.0) / ((alf + 2.0) * (alf + 3.0));
    customIp[1] = 2.0 * alf / (alf + 2.0);
    customIp[2] = 2.0 / (alf + 2.0);
    customIp[3] = -2.0 * alf / ((alf + 3.0) * (alf + 2.0));
}
