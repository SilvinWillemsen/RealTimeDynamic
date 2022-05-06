/*
  ==============================================================================

    ControlPanel.cpp
    Created: 7 Feb 2022 5:19:05pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include <JuceHeader.h>
#include "ControlPanel.h"

//==============================================================================
ControlPanel::ControlPanel (ChangeListener* changeListener)
{
    // In your constructor, you should add any child components, and
    // initialise any special settings that your component needs.
    addChangeListener (changeListener);
}

ControlPanel::~ControlPanel()
{
}

void ControlPanel::paint (juce::Graphics& g)
{
    /* This demo code just fills the component's background and
       draws some placeholder text to get you started.

       You should replace everything in this method with your own
       drawing code..
    */

//    g.fillAll (Colours::green);   // clear the background

}

void ControlPanel::resized()
{
    // This method is where you should set the bounds of any child
    // components that your component contains..
    Rectangle<int> controlArea = getLocalBounds().reduced (Global::margin);
    int sliderWidth = getWidth() / sliders.size();
    for (int i = 0; i < sliders.size(); ++i)
    {
        Rectangle<int> paramArea = controlArea.removeFromLeft (sliderWidth).reduced (Global::margin, 0);
        labels[i]->setBounds (paramArea.removeFromTop (paramArea.getHeight()*0.25));
        sliders[i]->setBounds (paramArea);
    }
}

void ControlPanel::refreshSliders (NamedValueSet& parameters)
{
    sliders.clearQuick (true);
    labels.clearQuick (true);
    for (int i = 0; i < parameters.size(); ++i)
    {
        sliders.add (new Slider (Slider::LinearHorizontal, Slider::TextBoxBelow));
        Slider* newSlider = sliders[sliders.size()-1];
        
        labels.add (new Label (parameters.getName(i).toString(), parameters.getName(i).toString()));
        Label* newLabel = labels[labels.size()-1];
        newLabel->setColour(Label::textColourId, Colours::white);
        Font font (18.0f);
        newLabel->setJustificationType (Justification::centred);
        newLabel->setFont (font);
        addAndMakeVisible (newLabel);
        
        const String sliderName = parameters.getName(i).getCharPointer();
        newSlider->setName (sliderName);
        
        double val = *parameters.getVarPointerAt(i);
        if (val != 0)
        {
            if (parameters.getName(i).toString() == "sigma0")
                newSlider->setRange (0, val * 2.0);
            else if (parameters.getName(i).toString() == "E")
                newSlider->setRange (1e9, 4e13);
            else if (parameters.getName(i).toString() == "T")
                newSlider->setRange (0, val * 2.0);
            else if (parameters.getName(i).toString() == "L")
                newSlider->setRange (0.1, val * 2.0);
            else if (parameters.getName(i).toString() == "sigma1") // have a minimum value for sigma1 to prevent artefacts
                newSlider->setRange (Global::sig1min, val * 10);
            else if (parameters.getName(i).toString() == "r") // have a minimum value for sigma1 to prevent artefacts
                newSlider->setRange (val * 0.5, val * 2.0);
            else
                newSlider->setRange (val * 0.5, val * 2.0);
            
            newSlider->setValue (val);
            newSlider->setSkewFactorFromMidPoint (val);
        }
        newSlider->addListener (this);
        addAndMakeVisible (newSlider);
    }
    resized();
}

void ControlPanel::sliderValueChanged (Slider* slider)
{
    changedParameterName = slider->getName();
    changedParameterValue = slider->getValue();
    changedParameterIdx = sliders.indexOf (slider);
    sendChangeMessage();
}
