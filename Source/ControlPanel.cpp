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
    Rectangle<int> controlArea = getLocalBounds();
    int sliderHeight = getHeight() / sliders.size();
    for (auto s : sliders)
        s->setBounds(controlArea.removeFromTop(sliderHeight));
}

void ControlPanel::refreshSliders (NamedValueSet& parameters)
{
    sliders.clearQuick (true);
    for (int i = 0; i < parameters.size(); ++i)
    {
        sliders.add (new Slider (Slider::LinearHorizontal, Slider::TextBoxAbove));
        Slider* newSlider = sliders[sliders.size()-1];
        
        const String sliderName = parameters.getName(i).getCharPointer();
        newSlider->setName (sliderName);
        
        double val = *parameters.getVarPointerAt(i);
        if (val != 0)
        {
            if (parameters.getName(i).toString() == "E" ||
                parameters.getName(i).toString() == "I" ||
                parameters.getName(i).toString() == "sigma0" ||
                parameters.getName(i).toString() == "sigma1")
                newSlider->setRange (0, val * 2.0);
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
