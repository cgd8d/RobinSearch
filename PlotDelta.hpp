#ifndef PlotDelta_hpp
#define PlotDelta_hpp

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

/*
Accept data points and pipe them into gnuplot.
The goal is to plot log delta vs log log n,
as in figure 3 (top) of Briggs 2006.
*/
struct PlotDeltaStruct
{

    double loglogn_step;
    double loglogn_next;
    double last_yval; // log delta + loglogn/2
    double max_yval_change;

    // Store values to float precision.
    // This is only for plotting, and
    // the size of the serialized file
    // can become annoyingly big.
    std::vector<std::pair<float, float>> data;

    PlotDeltaStruct()

    {
        // Configure to save 20
        // points per unit of loglogn.
        loglogn_step = 0.05;
        loglogn_next = 5;

        // Additionally, we want to save
        // a point whenever
        // log delta + loglogn/2 has changed
        // by enough.
        // Configure to save 20
        // points per y-axis grid mark (0.01).
        last_yval = -1000; // dummy value.
        max_yval_change = 0.0005;
    }

    void AddPoint(double loglogn, double logdelta)
    {
        double this_yval = logdelta+loglogn/2;
        if(loglogn >= loglogn_next
            or std::fabs(this_yval-last_yval) >= max_yval_change)
        {
            data.emplace_back(loglogn, logdelta);
            loglogn_next = loglogn + loglogn_step;
            last_yval = this_yval;
        }
    }

    ~PlotDeltaStruct()
    {
        FILE* plotpipe = popen("gnuplot", "w");
        if(not plotpipe)
        {
            std::cerr << "Failed to open gnuplot pipe." << std::endl;
            //can't throw exception from destructor, no matter
            //how much we want to.
            //throw std::runtime_error("Failed to open gnuplot pipe.");
            return;
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
            //can't throw exception from destructor, no matter
            //how much we want to.
            //throw std::runtime_error("Failed to open gnuplot pipe.");
            pclose(plotpipe);
            return;
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
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & loglogn_step;
        ar & loglogn_next;
        ar & data;
    }
};

#endif
