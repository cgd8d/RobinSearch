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
        fprintf(plotpipe, "set xlabel 'log(log(N))'\n");
        fprintf(plotpipe, "set ylabel '\\delta(N)\n");
        fprintf(plotpipe, "plot '-'\n");
    }

    ~PlotDeltaStruct()
    {
        fprintf(plotpipe, "e\n");
        fflush(plotpipe);
        pclose(plotpipe);
    }
};

#endif
