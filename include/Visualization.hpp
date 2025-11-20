#ifndef VISUALIZATION_HPP
#define VISUALIZATION_HPP

#include "ReservoirSim.hpp"
#include <string>
#include <vector>
#include <map>
#include <array>

namespace ResSim {

// ============================================================================
// PPM Image Generation Classes
// ============================================================================

struct Color {
    unsigned char r, g, b;
    
    Color(unsigned char red = 0, unsigned char green = 0, unsigned char blue = 0)
        : r(red), g(green), b(blue) {}
    
    static Color lerp(const Color& c1, const Color& c2, double t);
};

class ColorMap {
public:
    virtual Color getColor(double value) const = 0;
    virtual ~ColorMap() = default;
};

class JetColorMap : public ColorMap {
public:
    Color getColor(double value) const override;
};

class ViridisColorMap : public ColorMap {
public:
    Color getColor(double value) const override;
};

class PressureColorMap : public ColorMap {
public:
    Color getColor(double value) const override;
};

class SubsidenceColorMap : public ColorMap {
public:
    Color getColor(double value) const override;
};

// Scalable bitmap font for text rendering (base 6x8, scalable)
class BitmapFont {
public:
    static void drawChar(std::vector<std::vector<Color>>& image, 
                        int x, int y, char c, const Color& color, int scale = 3);
    static void drawString(std::vector<std::vector<Color>>& image,
                          int x, int y, const std::string& text, 
                          const Color& color, int scale = 3);
    static int getCharWidth(int scale = 3) { return 6 * scale; }
    static int getCharHeight(int scale = 3) { return 8 * scale; }
};

// ============================================================================
// Plot Generator for 2D Raster Images
// ============================================================================

class PlotGenerator2D {
public:
    PlotGenerator2D(int width, int height, int margin_left = 120, int margin_right = 220,
                   int margin_top = 100, int margin_bottom = 80);
    
    // Draw 2D field data
    void drawField(const std::vector<std::vector<double>>& data,
                   const ColorMap& cmap, double vmin, double vmax);
    
    // Draw colorbar with labels
    void drawColorbar(const ColorMap& cmap, double vmin, double vmax,
                     const std::string& label = "");
    
    // Draw well symbols
    void drawWell(double x_pos, double y_pos, double Lx, double Ly,
                  bool is_injector, const std::string& label = "");
    
    // Draw grid lines
    void drawGrid(int num_x_lines, int num_y_lines);
    
    // Draw rectangles and borders
    void drawRect(int x0, int y0, int x1, int y1, const Color& color, int thickness = 1);
    
    // Add title and axis labels
    void setTitle(const std::string& title);
    void setXLabel(const std::string& label);
    void setYLabel(const std::string& label);
    
    // Write to file
    void writePPM(const std::string& filename) const;
    bool convertToPNG(const std::string& ppm_file, const std::string& png_file) const;
    void writeImage(const std::string& base_filename) const;  // Tries PNG, falls back to PPM
    
    // Get dimensions
    int getFieldWidth() const { return field_width_; }
    int getFieldHeight() const { return field_height_; }
    
private:
    int plot_width_, plot_height_;
    int margin_l_, margin_r_, margin_t_, margin_b_;
    int field_width_, field_height_;
    std::vector<std::vector<Color>> image_;
    std::string title_;
    std::string xlabel_;
    std::string ylabel_;
    
    void drawText(int x, int y, const std::string& text, const Color& color);
};

// ============================================================================
// Multi-panel Plot Generator
// ============================================================================

class MultiPanelPlot {
public:
    struct Panel {
        std::vector<std::vector<double>> data;
        const ColorMap* colormap;
        double vmin, vmax;
        std::string title;
    };
    
    MultiPanelPlot(int num_rows, int num_cols, int panel_width, int panel_height);
    
    void setPanel(int row, int col, const Panel& panel);
    void addWellToPanel(int row, int col, double x, double y, double Lx, double Ly,
                       bool is_injector, const std::string& label = "");
    
    void setMainTitle(const std::string& title);
    void writeImage(const std::string& filename);
    
private:
    int num_rows_, num_cols_;
    int panel_width_, panel_height_;
    std::vector<std::vector<Panel>> panels_;
    std::string main_title_;
};

// ============================================================================
// Original Visualization Class (Extended)
// ============================================================================

class Visualization {
public:
    Visualization();
    ~Visualization();
    
    // Output formats
    enum class OutputFormat {
        VTK,
        VTU,
        HDF5,
        ECLIPSE,
        ENSIGHT,
        XDMF
    };
    
    // Set output format
    void setFormat(OutputFormat format);
    void setOutputDirectory(const std::string& dir);
    
    // Write solutions
    PetscErrorCode writeSolution(Vec solution, DM dm, double time, int step,
                                const std::string& filename);
    
    // Write specific fields
    PetscErrorCode writePressureField(Vec solution, DM dm, int step);
    PetscErrorCode writeSaturationField(Vec solution, DM dm, int step);
    PetscErrorCode writeDisplacementField(Vec solution, DM dm, int step);
    PetscErrorCode writeTemperatureField(Vec solution, DM dm, int step);
    
    // Write derived quantities
    PetscErrorCode writeVelocityField(Vec solution, DM dm, int step);
    PetscErrorCode writeStressField(Vec solution, DM dm, int step);
    PetscErrorCode writeStrainField(Vec solution, DM dm, int step);
    
    // 2D plots and charts
    void plotPressureHistory(const std::vector<double>& times,
                            const std::vector<double>& pressures,
                            const std::string& filename);
    
    void plotProductionHistory(const std::vector<double>& times,
                              const std::map<std::string, std::vector<double>>& rates,
                              const std::string& filename);
    
    void plotConvergenceHistory(const std::vector<int>& iterations,
                               const std::vector<double>& residuals,
                               const std::string& filename);
    
    // Contour plots
    void plotContour2D(const std::vector<double>& x,
                      const std::vector<double>& y,
                      const std::vector<double>& values,
                      const std::string& title,
                      const std::string& filename);
    
    // 3D visualization
    void plot3DField(Vec solution, DM dm, const std::string& field_name,
                    const std::string& filename);
    
    // Generate summary tables
    void writeSummaryTable(const std::map<std::string, double>& metrics,
                          const std::string& filename);
    
    // Performance plots
    void plotScalability(const std::vector<int>& num_procs,
                        const std::vector<double>& solve_times,
                        const std::string& filename);
    
    void plotConvergenceRate(const std::vector<double>& mesh_sizes,
                            const std::vector<double>& errors,
                            double expected_rate,
                            const std::string& filename);
    
    // Animation generation
    void generateAnimation(const std::vector<std::string>& vtk_files,
                          const std::string& output_file);
    
    // PPM/PNG 2D visualization
    void plot2DField(const std::vector<std::vector<double>>& data,
                    double Lx, double Ly,
                    const std::string& title,
                    const std::string& colormap_name,
                    const std::string& filename,
                    double vmin = 0.0, double vmax = 1.0);
    
private:
    OutputFormat format;
    std::string output_dir;
    
    // VTK writing helpers
    PetscErrorCode writeVTK(Vec solution, DM dm, const std::string& filename);
    PetscErrorCode writeVTU(Vec solution, DM dm, const std::string& filename);
    
    // HDF5 writing helpers
    PetscErrorCode writeHDF5(Vec solution, DM dm, const std::string& filename);
    
    // Plotting helpers (using gnuplot or matplotlib via system calls)
    void generateGnuplotScript(const std::string& script_content,
                              const std::string& output_file);
};

} // namespace ResSim

#endif // VISUALIZATION_HPP
