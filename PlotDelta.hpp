#ifndef PlotDelta_hpp
#define PlotDelta_hpp

#include <iostream>
#include <vector>
#include <utility>

/*
Accept days points and pipe them into gnuplot.
The goal is to plot log delta vs log log n,
as in figure 3 (top) of Briggs 2006.
*/
struct PlotDeltaStruct
{

    double loglogn_step;
    double loglogn_next;
    std::vector<std::pair<double, double>> data;

    PlotDeltaStruct()
    
    {

        loglogn_step = 0.001;
        loglogn_next = 5;
    }

    void AddPoint(double loglogn, double logdelta)
    {
        if(loglogn >= loglogn_next)
        {
            data.push_back(loglogn, logdelta);
            loglogn_next = loglogn + loglogn_step;
        }
    }

    ~PlotDeltaStruct()
    {
        FILE* plotpipe = popen("gnuplot", "w");
        if(not plotpipe)
        {
            std::cerr << "Failed to open gnuplot pipe." << std::endl;
            throw std::runtime_error("Failed to open gnuplot pipe.");
        }
        fprintf(plotpipe, "set term pngcairo\n");
        fprintf(plotpipe, "set output 'DeltaPlot.png'\n");
        fprintf(plotpipe, "set xlabel 'log(log(N))'\n");
        fprintf(plotpipe, "set ylabel 'log \\delta(N)\n");
        fprintf(plotpipe, "set xrange [15:*]\n");
        fprintf(plotpipe, "set xtics 2\n");
        fprintf(plotpipe, "set mxtics 2\n");
        fprintf(plotpipe, "set ytics 1\n");
        fprintf(plotpipe, "set mytics 2\n");
        fprintf(plotpipe, "set grid xtics mxtics ytics mytics\n");
        fprintf(plotpipe, "plot '-' with points pt 0\n");

        FILE* plotpipe2 = popen("gnuplot", "w");
        if(not plotpipe2)
        {
            std::cerr << "Failed to open gnuplot pipe." << std::endl;
            throw std::runtime_error("Failed to open gnuplot pipe.");
        }
        fprintf(plotpipe2, "set term pngcairo\n");
        fprintf(plotpipe2, "set output 'DeltaPlotScaled.png'\n");
        fprintf(plotpipe2, "set xlabel 'log(log(N))'\n");
        fprintf(plotpipe2, "set ylabel 'log \\delta(N)+log(log(N))/2-0.323336\n");
        fprintf(plotpipe2, "set xrange [11:*]\n");
        fprintf(plotpipe2, "set yrange [*:0.08]\n");
        fprintf(plotpipe2, "set xtics 5\n");
        fprintf(plotpipe2, "set mxtics 5\n");
        fprintf(plotpipe2, "set ytics .02\n");
        fprintf(plotpipe2, "set mytics 2\n");
        fprintf(plotpipe2, "set grid xtics mxtics ytics mytics\n");
        fprintf(plotpipe2, "plot '-' with points pt 0\n");

        for(size_t i = 0;
            i < data.size();
            i++)
        {
            double loglogn = data[i].first;
            double logdelta = data[i].second;
            fprintf(plotpipe, "%.5g %.5g\n", loglogn, logdelta);
            fprintf(plotpipe2, "%.5g %.5g\n", loglogn, logdelta+loglogn/2-0.323336);
        }
        
        
        
        fprintf(plotpipe, "e\n");
        fflush(plotpipe);
        pclose(plotpipe);
        fprintf(plotpipe2, "e\n");
        fflush(plotpipe2);
        pclose(plotpipe2);
    }
};

#endif
