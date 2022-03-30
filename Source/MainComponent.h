#pragma once

#include <JuceHeader.h>
#include "Global.h"
#include "ControlPanel.h"
#include "DynamicStiffString.h"

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent  : public AudioAppComponent,
                       public ChangeListener,
                       public Timer, // for graphics refresh
                       public KeyListener
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent() override;

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint (juce::Graphics& g) override;
    void resized() override;

    void timerCallback() override;
    
    void changeListenerCallback (ChangeBroadcaster* changeBroadcaster) override;
    
    bool keyPressed (const KeyPress& key, Component* originatingComponent) override;
private:
    //==============================================================================
    // Your private member variables go here...
    std::unique_ptr<ControlPanel> controlPanel;
    
    std::unique_ptr<DynamicStiffString> dynamicStiffString;

    std::mutex audioMutex;
    bool parameterChangedFlag = false;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
