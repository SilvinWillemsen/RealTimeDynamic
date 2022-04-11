/*
  ==============================================================================

    Global.h
    Created: 10 Feb 2022 9:39:22am
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once
//#define RECORD
namespace Global {
#ifdef RECORD
    static int samplesToRecord = 500;

#endif
    static int margin = 10;
    static double NmaxChange = 1.0 / 20.0;
    static double sig1min = 0.0002;
    static float boundaryEllRad = 10.0;
    // limiter
    static double limit (double val, double min, double max)
    {
        if (val < min)
        {
            val = min;
            return val;
        }
        else if (val > max)
        {
            val = max;
            return val;
        }
        return val;
    }
}
