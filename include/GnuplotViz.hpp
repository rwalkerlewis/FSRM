#ifndef GNUPLOT_VIZ_HPP
#define GNUPLOT_VIZ_HPP

#include <vector>
#include <string>
#include <map>

namespace ResSim {

/**
 * @brief High-quality plotting using gnuplot-iostream
 * 
 * Generates publication-quality plots with proper:
 * - Colormaps (jet, viridis, plasma, etc.)
 * - Axis labels and titles
 * - Colorbars with tick labels
 * - Well symbols and annotations
 * - Grid lines and contours
 */
class GnuplotViz {
public:
    GnuplotViz(const std::string& output_dir = "output");
    ~GnuplotViz();
    
    // 2D field plots
    void plot2DField(const std::vector<std::vector<double>>& data,
                    double Lx, double Ly,
                    const std::string& title,
                    const std::string& xlabel,
                    const std::string& ylabel,
                    const std::string& colormap,
                    const std::string& filename,
                    double vmin = 0.0, double vmax = 1.0);
    
    // Add wells to last plot
    struct WellMarker {
        double x, y;
        std::string label;
        bool is_injector;
    };
    void addWells(const std::vector<WellMarker>& wells);
    
    // Line plots
    void plotLines(const std::vector<double>& x,
                  const std::map<std::string, std::vector<double>>& y_data,
                  const std::string& title,
                  const std::string& xlabel,
                  const std::string& ylabel,
                  const std::string& filename);
    
    // Multi-panel plots (2x2 grid)
    void plotMultiPanel(const std::vector<std::vector<std::vector<double>>>& data_panels,
                       const std::vector<std::string>& titles,
                       double Lx, double Ly,
                       const std::string& colormap,
                       const std::string& filename,
                       const std::vector<std::pair<double,double>>& ranges);
    
    void setOutputDir(const std::string& dir) { output_dir_ = dir; }
    
private:
    std::string output_dir_;
    std::vector<WellMarker> current_wells_;
    
    // Helper to write gnuplot data files
    void writeDataMatrix(const std::vector<std::vector<double>>& data,
                        double Lx, double Ly,
                        const std::string& filename);
    
    std::string getPalette(const std::string& colormap);
};

} // namespace ResSim

#endif // GNUPLOT_VIZ_HPP
