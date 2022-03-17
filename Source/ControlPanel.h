/*
  ==============================================================================

    ControlPanel.h
    Created: 7 Feb 2022 5:19:05pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
*/
class ControlPanel  : public juce::Component, public Slider::Listener, public ChangeBroadcaster
{
public:
    ControlPanel (ChangeListener* changeListener);
    ~ControlPanel() override;

    void paint (juce::Graphics&) override;
    void resized() override;
    
    void refreshSliders (NamedValueSet& parameters);
    
    void sliderValueChanged (Slider* slider) override;

    String getChangedParameterName() { return changedParameterName; };
    double getChangedParameterValue() { return changedParameterValue; };
    int getChangedParameterIdx() { return changedParameterIdx; };
private:
    OwnedArray<Slider> sliders;
    OwnedArray<Label> labels;

    double changedParameterValue = -1.0;
    String changedParameterName = "";
    int changedParameterIdx = -1;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ControlPanel)
};
