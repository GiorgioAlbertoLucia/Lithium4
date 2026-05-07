#pragma once

#include <TH1.h>
#include <TGraph.h>

float readFromTH1(TH1& histogram, const float x)
{
    return histogram.GetBinContent(histogram.FindBin(x));
}

float readFromTGraph(TGraph * graph, const float x0)
{
    double x, y;
    int best = 0;
    double dxmin = 1e300;

    for (int i = 0; i < graph->GetN(); ++i) {
        graph->GetPoint(i, x, y);
        double dx = std::abs(x - x0);
        if (dx < dxmin) {
            dxmin = dx;
            best = i;
        }
    }

    graph->GetPoint(best, x, y);
    return y;
}

