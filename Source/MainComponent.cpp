#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent()
{

    // Some platforms require permissions to open input channels so request that here
    if (juce::RuntimePermissions::isRequired (juce::RuntimePermissions::recordAudio)
        && ! juce::RuntimePermissions::isGranted (juce::RuntimePermissions::recordAudio))
    {
        juce::RuntimePermissions::request (juce::RuntimePermissions::recordAudio,
                                           [&] (bool granted) { setAudioChannels (granted ? 2 : 0, 2); });
    }
    else
    {
        // Specify the number of input and output channels that we want to open
        setAudioChannels (0, 2);
    }
}

MainComponent::~MainComponent()
{
    // Stop the graphics update
    stopTimer();

    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    //// Set the paramters ///
    NamedValueSet parameters;
    
    // parameters you'll use to initialise more than one other parameter should be defined here
    double r = 0.0005;
    
    parameters.set ("L", 2);
    parameters.set ("rho", 7850);
    parameters.set ("r", 0.0005);
    parameters.set ("T", 299.75);
    parameters.set ("E", 2e11);
    parameters.set ("sigma0", 2);
    parameters.set ("sigma1", 0.005);
    
    //// Initialise an instance of the DynamicStiffString class ////
    dynamicStiffString = std::make_unique<DynamicStiffString> (parameters, 1.0 / sampleRate);
    
    addAndMakeVisible (dynamicStiffString.get()); // add the string to the application
    
    controlPanel = std::make_unique<ControlPanel> (this);
    
    addAndMakeVisible (controlPanel.get());
    
    controlPanel->refreshSliders (parameters);
    // Call resized again as our components need a sample rate before they can get initialised.
    setSize (800, 600);

    startTimerHz (15); // start the timer (15 Hz is a nice tradeoff between CPU usage and update speed)
    
//     dynamicStiffString->refreshParameter (0, 1);
}

void MainComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    bufferToFill.clearActiveBufferRegion();

    int numChannels = bufferToFill.buffer->getNumChannels();
    
    // Get pointers to output locations
    float* const channelData1 = bufferToFill.buffer->getWritePointer (0, bufferToFill.startSample);
    float* const channelData2 = numChannels > 1 ? bufferToFill.buffer->getWritePointer (1, bufferToFill.startSample) : nullptr;

    float output = 0.0;

    std::vector<float* const*> curChannel {&channelData1, &channelData2};
    
    // only do control stuff out of the buffer (at least work with flags so that control doesn't interfere with the scheme calculation)
    if (dynamicStiffString->shouldExcite())
        dynamicStiffString->excite();
        
    audioMutex.lock();
    for (int i = 0; i < bufferToFill.numSamples; ++i)
    {
#ifndef RECORD
        if (parameterChangedFlag)
        {
            dynamicStiffString->refreshParameter(controlPanel->getChangedParameterIdx(), controlPanel->getChangedParameterValue());
            parameterChangedFlag = false;
        }
#endif
        dynamicStiffString->refreshCoefficients(); // for every loop for now
        dynamicStiffString->calculateScheme();
        dynamicStiffString->updateStates();
        
        output = dynamicStiffString->getOutput (0.13); // get output at 0.13L of the string
        for (int channel = 0; channel < numChannels; ++channel)
            curChannel[channel][0][i] = Global::limit (output, -1.0, 1.0);
    }
    
    audioMutex.unlock();

}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint (juce::Graphics& g)
{
}

void MainComponent::resized()
{
    
    Rectangle<int> totalArea = getLocalBounds();
    
    if (controlPanel != nullptr)
        controlPanel->setBounds (totalArea.removeFromBottom (100));
    
    // put the string in the application
    if (dynamicStiffString != nullptr)
        dynamicStiffString->setBounds (totalArea);
    
        
}

void MainComponent::timerCallback()
{
    repaint(); // update the graphics X times a second
}

void MainComponent::changeListenerCallback (ChangeBroadcaster* changeBroadcaster)
{
    if (changeBroadcaster == controlPanel.get())
    {
        parameterChangedFlag = true;
        std::cout << controlPanel->getChangedParameterName() << ": " << controlPanel->getChangedParameterValue() << std::endl;

    }
}
