#ifndef PlotDelta_hpp
#define PlotDelta_hpp

#include <iostream>

/*
Accept days points and pipe them into gnuplot.
The goal is to plot log delta vs log log n,
as in figure 3 (top) of Briggs 2006.
*/
struct PlotDeltaStruct
{
    FILE* plotpipe;

    PlotDeltaStruct()
    {
        plotpipe = popen("gnuplot", "w");
        if(not plotpipe)
        {
            std::cerr << "Failed to open gnuplot pipe." << std::endl;
            throw std::runtime_error("Failed to open gnuplot pipe.");
        }
        fprintf(plotpipe, "set output 'DeltaPlot.png'\n");
        fprintf(plotpipe, "set xlabel 'log(log(N))'\n");
        fprintf(plotpipe, "set ylabel 'log \\delta(N)\n");
        fprintf(plotpipe, "plot '-'\n");
    }

    void AddPoint(double loglogn, double logdelta)
    {
        fprintf(plotpipe, "%.5g %.5g\n", loglogn, logdelta);
    }

    ~PlotDeltaStruct()
    {
        fprintf(plotpipe, "e\n");
        fflush(plotpipe);
        pclose(plotpipe);
    }
};

#endif
