#include "Visualization.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <petscviewerhdf5.h>

namespace FSRM {

// ============================================================================
// Color and ColorMap Implementations
// ============================================================================

Color Color::lerp(const Color& c1, const Color& c2, double t) {
    t = std::max(0.0, std::min(1.0, t));
    return Color(
        static_cast<unsigned char>(c1.r + t * (c2.r - c1.r)),
        static_cast<unsigned char>(c1.g + t * (c2.g - c1.g)),
        static_cast<unsigned char>(c1.b + t * (c2.b - c1.b))
    );
}

Color JetColorMap::getColor(double value) const {
    value = std::max(0.0, std::min(1.0, value));
    
    if (value < 0.25) {
        double t = value / 0.25;
        return Color::lerp(Color(0, 0, 128), Color(0, 0, 255), t);
    } else if (value < 0.5) {
        double t = (value - 0.25) / 0.25;
        return Color::lerp(Color(0, 0, 255), Color(0, 255, 255), t);
    } else if (value < 0.75) {
        double t = (value - 0.5) / 0.25;
        return Color::lerp(Color(0, 255, 255), Color(255, 255, 0), t);
    } else {
        double t = (value - 0.75) / 0.25;
        return Color::lerp(Color(255, 255, 0), Color(255, 0, 0), t);
    }
}

Color ViridisColorMap::getColor(double value) const {
    value = std::max(0.0, std::min(1.0, value));
    
    const std::vector<Color> palette = {
        Color(68, 1, 84),
        Color(59, 82, 139),
        Color(33, 145, 140),
        Color(94, 201, 98),
        Color(253, 231, 37)
    };
    
    int n = palette.size() - 1;
    double scaled = value * n;
    int idx = static_cast<int>(scaled);
    idx = std::max(0, std::min(n - 1, idx));
    
    double t = scaled - idx;
    return Color::lerp(palette[idx], palette[idx + 1], t);
}

Color PressureColorMap::getColor(double value) const {
    value = std::max(0.0, std::min(1.0, value));
    
    if (value < 0.25) {
        double t = value * 4.0;
        return Color::lerp(Color(0, 0, 255), Color(0, 255, 255), t);
    } else if (value < 0.5) {
        double t = (value - 0.25) * 4.0;
        return Color::lerp(Color(0, 255, 255), Color(0, 255, 0), t);
    } else if (value < 0.75) {
        double t = (value - 0.5) * 4.0;
        return Color::lerp(Color(0, 255, 0), Color(255, 255, 0), t);
    } else {
        double t = (value - 0.75) * 4.0;
        return Color::lerp(Color(255, 255, 0), Color(255, 0, 0), t);
    }
}

Color SubsidenceColorMap::getColor(double value) const {
    value = std::max(0.0, std::min(1.0, value));
    
    if (value < 0.25) {
        double t = value * 4.0;
        return Color::lerp(Color(255, 255, 255), Color(255, 255, 100), t);
    } else if (value < 0.5) {
        double t = (value - 0.25) * 4.0;
        return Color::lerp(Color(255, 255, 100), Color(255, 165, 0), t);
    } else if (value < 0.75) {
        double t = (value - 0.5) * 4.0;
        return Color::lerp(Color(255, 165, 0), Color(255, 50, 50), t);
    } else {
        double t = (value - 0.75) * 4.0;
        return Color::lerp(Color(255, 50, 50), Color(139, 0, 0), t);
    }
}

// ============================================================================
// Bitmap Font for Labels with Scaling Support
// ============================================================================

void BitmapFont::drawChar(std::vector<std::vector<Color>>& image,
                         int x, int y, char c, const Color& color, int scale) {
    // Simple 6x8 bitmap font patterns
    static const std::map<char, std::vector<uint8_t>> patterns = {
        {' ', {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00}},
        {'0', {0x3C,0x66,0x6E,0x76,0x66,0x66,0x3C,0x00}},
        {'1', {0x18,0x38,0x18,0x18,0x18,0x18,0x7E,0x00}},
        {'2', {0x3C,0x66,0x06,0x0C,0x18,0x30,0x7E,0x00}},
        {'3', {0x3C,0x66,0x06,0x1C,0x06,0x66,0x3C,0x00}},
        {'4', {0x0C,0x1C,0x2C,0x4C,0x7E,0x0C,0x0C,0x00}},
        {'5', {0x7E,0x60,0x7C,0x06,0x06,0x66,0x3C,0x00}},
        {'6', {0x3C,0x60,0x60,0x7C,0x66,0x66,0x3C,0x00}},
        {'7', {0x7E,0x06,0x0C,0x18,0x30,0x30,0x30,0x00}},
        {'8', {0x3C,0x66,0x66,0x3C,0x66,0x66,0x3C,0x00}},
        {'9', {0x3C,0x66,0x66,0x3E,0x06,0x0C,0x38,0x00}},
        {'.', {0x00,0x00,0x00,0x00,0x00,0x18,0x18,0x00}},
        {'-', {0x00,0x00,0x00,0x7E,0x00,0x00,0x00,0x00}},
        {'(', {0x0C,0x18,0x30,0x30,0x30,0x18,0x0C,0x00}},
        {')', {0x30,0x18,0x0C,0x0C,0x0C,0x18,0x30,0x00}},
        {':', {0x00,0x18,0x18,0x00,0x18,0x18,0x00,0x00}},
        {'/', {0x02,0x06,0x0C,0x18,0x30,0x60,0x40,0x00}},
        {'e', {0x00,0x00,0x3C,0x66,0x7E,0x60,0x3C,0x00}},
        {'a', {0x00,0x00,0x3C,0x06,0x3E,0x66,0x3E,0x00}},
        {'P', {0x7C,0x66,0x66,0x7C,0x60,0x60,0x60,0x00}},
        {'r', {0x00,0x00,0x5C,0x66,0x60,0x60,0x60,0x00}},
        {'s', {0x00,0x00,0x3E,0x60,0x3C,0x06,0x7C,0x00}},
        {'u', {0x00,0x00,0x66,0x66,0x66,0x66,0x3E,0x00}},
        {'M', {0x63,0x77,0x7F,0x6B,0x63,0x63,0x63,0x00}},
        {'D', {0x78,0x6C,0x66,0x66,0x66,0x6C,0x78,0x00}},
        {'W', {0x63,0x63,0x63,0x6B,0x7F,0x77,0x63,0x00}},
        {'H', {0x66,0x66,0x66,0x7E,0x66,0x66,0x66,0x00}},
        {'o', {0x00,0x00,0x3C,0x66,0x66,0x66,0x3C,0x00}},
        {'i', {0x18,0x00,0x38,0x18,0x18,0x18,0x3C,0x00}},
        {'z', {0x00,0x00,0x7E,0x0C,0x18,0x30,0x7E,0x00}},
        {'t', {0x10,0x30,0x7C,0x30,0x30,0x34,0x18,0x00}},
        {'l', {0x38,0x18,0x18,0x18,0x18,0x18,0x3C,0x00}},
        {'n', {0x00,0x00,0x5C,0x66,0x66,0x66,0x66,0x00}},
        {'d', {0x06,0x06,0x3E,0x66,0x66,0x66,0x3E,0x00}},
        {'T', {0x7E,0x18,0x18,0x18,0x18,0x18,0x18,0x00}},
        {'m', {0x00,0x00,0x6C,0x7E,0x6B,0x63,0x63,0x00}},
        {'S', {0x3C,0x66,0x60,0x3C,0x06,0x66,0x3C,0x00}},
        {'b', {0x60,0x60,0x7C,0x66,0x66,0x66,0x7C,0x00}},
        {'c', {0x00,0x00,0x3C,0x66,0x60,0x66,0x3C,0x00}},
        {'F', {0x7E,0x60,0x60,0x7C,0x60,0x60,0x60,0x00}},
        {'I', {0x3C,0x18,0x18,0x18,0x18,0x18,0x3C,0x00}},
        {'V', {0x66,0x66,0x66,0x66,0x66,0x3C,0x18,0x00}},
        {'w', {0x00,0x00,0x63,0x63,0x6B,0x7F,0x36,0x00}},
        {'N', {0x66,0x76,0x7E,0x6E,0x66,0x66,0x66,0x00}},
    };
    
    auto it = patterns.find(c);
    if (it == patterns.end()) return;
    
    const auto& pattern = it->second;
    for (int row = 0; row < 8; ++row) {
        uint8_t bits = pattern[row];
        for (int col = 0; col < 6; ++col) {
            if (bits & (1 << (5 - col))) {
                // Draw scaled pixel
                for (int sy = 0; sy < scale; ++sy) {
                    for (int sx = 0; sx < scale; ++sx) {
                        int px = x + col * scale + sx;
                        int py = y + row * scale + sy;
                        if (py >= 0 && py < (int)image.size() && 
                            px >= 0 && px < (int)image[0].size()) {
                            image[py][px] = color;
                        }
                    }
                }
            }
        }
    }
}

void BitmapFont::drawString(std::vector<std::vector<Color>>& image,
                           int x, int y, const std::string& text,
                           const Color& color, int scale) {
    int cur_x = x;
    for (char c : text) {
        drawChar(image, cur_x, y, c, color, scale);
        cur_x += getCharWidth(scale);
    }
}

// ============================================================================\n// PlotGenerator2D Implementation
// ============================================================================

PlotGenerator2D::PlotGenerator2D(int width, int height, int margin_left,
                                 int margin_right, int margin_top, int margin_bottom)
    : plot_width_(width), plot_height_(height),
      margin_l_(margin_left), margin_r_(margin_right),
      margin_t_(margin_top), margin_b_(margin_bottom) {
    
    field_width_ = plot_width_ - margin_l_ - margin_r_;
    field_height_ = plot_height_ - margin_t_ - margin_b_;
    
    image_.resize(plot_height_, std::vector<Color>(plot_width_, Color(255, 255, 255)));
}

Visualization::Visualization() 
    : format(OutputFormat::HDF5), output_dir("output") {}

Visualization::~Visualization() {}

void Visualization::setFormat(OutputFormat fmt) {
    format = fmt;
}

void Visualization::setOutputDirectory(const std::string& dir) {
    output_dir = dir;
}

PetscErrorCode Visualization::writeSolution(Vec solution, DM dm, double time, int step,
                                           const std::string& filename) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    switch (format) {
        case OutputFormat::VTK:
        case OutputFormat::VTU:
            ierr = writeVTK(solution, dm, filename); CHKERRQ(ierr);
            break;
        case OutputFormat::HDF5:
            ierr = writeHDF5(solution, dm, filename); CHKERRQ(ierr);
            break;
        default:
            ierr = writeVTK(solution, dm, filename); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writeVTK(Vec solution, DM dm, const std::string& filename) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    PetscViewer viewer;
    std::string full_path = output_dir + "/" + filename;
    
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, full_path.c_str(), 
                             FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = VecView(solution, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writeVTU(Vec solution, DM dm, const std::string& filename) {
    return writeVTK(solution, dm, filename);
}

PetscErrorCode Visualization::writeHDF5(Vec solution, DM dm, const std::string& filename) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    PetscViewer viewer;
    std::string full_path = output_dir + "/" + filename + ".h5";
    
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, full_path.c_str(),
                              FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = VecView(solution, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writePressureField(Vec solution, DM dm, int step) {
    PetscFunctionBeginUser;
    
    std::ostringstream filename;
    filename << "pressure_" << std::setfill('0') << std::setw(6) << step << ".vtu";
    
    PetscErrorCode ierr = writeSolution(solution, dm, 0.0, step, filename.str());
    CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writeSaturationField(Vec solution, DM dm, int step) {
    PetscFunctionBeginUser;
    
    std::ostringstream filename;
    filename << "saturation_" << std::setfill('0') << std::setw(6) << step << ".vtu";
    
    PetscErrorCode ierr = writeSolution(solution, dm, 0.0, step, filename.str());
    CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writeDisplacementField(Vec solution, DM dm, int step) {
    PetscFunctionBeginUser;
    
    std::ostringstream filename;
    filename << "displacement_" << std::setfill('0') << std::setw(6) << step << ".vtu";
    
    PetscErrorCode ierr = writeSolution(solution, dm, 0.0, step, filename.str());
    CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writeTemperatureField(Vec solution, DM dm, int step) {
    PetscFunctionBeginUser;
    
    std::ostringstream filename;
    filename << "temperature_" << std::setfill('0') << std::setw(6) << step << ".vtu";
    
    PetscErrorCode ierr = writeSolution(solution, dm, 0.0, step, filename.str());
    CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writeVelocityField(Vec solution, DM dm, int step) {
    PetscFunctionBeginUser;
    // Compute velocity from pressure gradient and Darcy's law
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writeStressField(Vec solution, DM dm, int step) {
    PetscFunctionBeginUser;
    // Compute stress from displacement field
    PetscFunctionReturn(0);
}

PetscErrorCode Visualization::writeStrainField(Vec solution, DM dm, int step) {
    PetscFunctionBeginUser;
    // Compute strain from displacement field
    PetscFunctionReturn(0);
}

void Visualization::plotPressureHistory(const std::vector<double>& times,
                                       const std::vector<double>& pressures,
                                       const std::string& filename) {
    std::string full_path = output_dir + "/" + filename;
    
    // Write data
    std::ofstream data(full_path + ".dat");
    for (size_t i = 0; i < times.size(); ++i) {
        data << times[i] << " " << pressures[i] << "\n";
    }
    data.close();
    
    // Generate gnuplot script
    std::ofstream script(full_path + ".gp");
    script << "set terminal png size 800,600\n";
    script << "set output '" << full_path << ".png'\n";
    script << "set xlabel 'Time (s)'\n";
    script << "set ylabel 'Pressure (Pa)'\n";
    script << "set title 'Pressure History'\n";
    script << "set grid\n";
    script << "plot '" << full_path << ".dat' using 1:2 with lines lw 2 notitle\n";
    script.close();
    
    system(("gnuplot " + full_path + ".gp 2>/dev/null").c_str());
}

void Visualization::plotProductionHistory(const std::vector<double>& times,
                                         const std::map<std::string, std::vector<double>>& rates,
                                         const std::string& filename) {
    std::string full_path = output_dir + "/" + filename;
    
    // Write data
    std::ofstream data(full_path + ".dat");
    for (size_t i = 0; i < times.size(); ++i) {
        data << times[i];
        for (const auto& pair : rates) {
            if (i < pair.second.size()) {
                data << " " << pair.second[i];
            }
        }
        data << "\n";
    }
    data.close();
    
    // Generate gnuplot script
    std::ofstream script(full_path + ".gp");
    script << "set terminal png size 800,600\n";
    script << "set output '" << full_path << ".png'\n";
    script << "set xlabel 'Time (days)'\n";
    script << "set ylabel 'Rate (m³/day)'\n";
    script << "set title 'Production History'\n";
    script << "set grid\n";
    script << "set key outside right\n";
    
    script << "plot ";
    int col = 2;
    for (const auto& pair : rates) {
        if (col > 2) script << ", \\\n     ";
        script << "'" << full_path << ".dat' using 1:" << col 
               << " with lines lw 2 title '" << pair.first << "'";
        col++;
    }
    script << "\n";
    script.close();
    
    system(("gnuplot " + full_path + ".gp 2>/dev/null").c_str());
}

void Visualization::plotConvergenceHistory(const std::vector<int>& iterations,
                                          const std::vector<double>& residuals,
                                          const std::string& filename) {
    std::string full_path = output_dir + "/" + filename;
    
    std::ofstream data(full_path + ".dat");
    for (size_t i = 0; i < iterations.size(); ++i) {
        data << iterations[i] << " " << residuals[i] << "\n";
    }
    data.close();
    
    std::ofstream script(full_path + ".gp");
    script << "set terminal png size 800,600\n";
    script << "set output '" << full_path << ".png'\n";
    script << "set xlabel 'Iteration'\n";
    script << "set ylabel 'Residual'\n";
    script << "set logscale y\n";
    script << "set title 'Convergence History'\n";
    script << "set grid\n";
    script << "plot '" << full_path << ".dat' using 1:2 with linespoints lw 2 notitle\n";
    script.close();
    
    system(("gnuplot " + full_path + ".gp 2>/dev/null").c_str());
}

void Visualization::plotContour2D(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 const std::vector<double>& values,
                                 const std::string& title,
                                 const std::string& filename) {
    std::string full_path = output_dir + "/" + filename;
    
    // Write data in gnuplot matrix format
    std::ofstream data(full_path + ".dat");
    size_t nx = x.size();
    size_t ny = y.size();
    
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            size_t idx = i + j * nx;
            if (idx < values.size()) {
                data << x[i] << " " << y[j] << " " << values[idx] << "\n";
            }
        }
        data << "\n";
    }
    data.close();
    
    std::ofstream script(full_path + ".gp");
    script << "set terminal png size 800,600\n";
    script << "set output '" << full_path << ".png'\n";
    script << "set xlabel 'X'\n";
    script << "set ylabel 'Y'\n";
    script << "set title '" << title << "'\n";
    script << "set pm3d map\n";
    script << "set palette defined (0 'blue', 0.5 'white', 1 'red')\n";
    script << "splot '" << full_path << ".dat' using 1:2:3 with pm3d notitle\n";
    script.close();
    
    system(("gnuplot " + full_path + ".gp 2>/dev/null").c_str());
}

void Visualization::plot3DField(Vec solution, DM dm, const std::string& field_name,
                               const std::string& filename) {
    // 3D visualization would use ParaView/VTK
    writeVTK(solution, dm, filename);
}

void Visualization::writeSummaryTable(const std::map<std::string, double>& metrics,
                                     const std::string& filename) {
    std::string full_path = output_dir + "/" + filename;
    std::ofstream table(full_path);
    
    table << std::setw(30) << std::left << "Metric" 
          << std::setw(20) << std::right << "Value" << "\n";
    table << std::string(50, '-') << "\n";
    
    for (const auto& pair : metrics) {
        table << std::setw(30) << std::left << pair.first
              << std::setw(20) << std::right << std::scientific 
              << std::setprecision(6) << pair.second << "\n";
    }
    
    table.close();
}

void Visualization::plotScalability(const std::vector<int>& num_procs,
                                   const std::vector<double>& solve_times,
                                   const std::string& filename) {
    std::string full_path = output_dir + "/" + filename;
    
    std::ofstream data(full_path + ".dat");
    for (size_t i = 0; i < num_procs.size(); ++i) {
        double speedup = solve_times[0] / solve_times[i];
        double efficiency = speedup / num_procs[i];
        data << num_procs[i] << " " << speedup << " " 
             << num_procs[i] << " " << efficiency << "\n";
    }
    data.close();
    
    std::ofstream script(full_path + ".gp");
    script << "set terminal png size 1200,500\n";
    script << "set output '" << full_path << ".png'\n";
    script << "set multiplot layout 1,2\n";
    
    script << "set xlabel 'Number of Processors'\n";
    script << "set ylabel 'Speedup'\n";
    script << "set title 'Strong Scaling'\n";
    script << "set grid\n";
    script << "plot '" << full_path << ".dat' using 1:2 with linespoints lw 2 title 'Actual', \\\n";
    script << "     '" << full_path << ".dat' using 1:3 with lines dt 2 lw 2 title 'Ideal'\n";
    
    script << "set ylabel 'Parallel Efficiency'\n";
    script << "set title 'Parallel Efficiency'\n";
    script << "plot '" << full_path << ".dat' using 1:4 with linespoints lw 2 notitle\n";
    
    script << "unset multiplot\n";
    script.close();
    
    system(("gnuplot " + full_path + ".gp 2>/dev/null").c_str());
}

void Visualization::plotConvergenceRate(const std::vector<double>& mesh_sizes,
                                       const std::vector<double>& errors,
                                       double expected_rate,
                                       const std::string& filename) {
    std::string full_path = output_dir + "/" + filename;
    
    std::ofstream data(full_path + ".dat");
    for (size_t i = 0; i < mesh_sizes.size(); ++i) {
        double ideal_error = errors[0] * std::pow(mesh_sizes[i] / mesh_sizes[0], expected_rate);
        data << mesh_sizes[i] << " " << errors[i] << " " << ideal_error << "\n";
    }
    data.close();
    
    std::ofstream script(full_path + ".gp");
    script << "set terminal png size 800,600\n";
    script << "set output '" << full_path << ".png'\n";
    script << "set logscale xy\n";
    script << "set xlabel 'Mesh Size h'\n";
    script << "set ylabel 'Error'\n";
    script << "set title 'Convergence Rate'\n";
    script << "set grid\n";
    script << "plot '" << full_path << ".dat' using 1:2 with linespoints lw 2 title 'Computed', \\\n";
    script << "     '" << full_path << ".dat' using 1:3 with lines dt 2 lw 2 title 'O(h^" 
           << expected_rate << ")'\n";
    script.close();
    
    system(("gnuplot " + full_path + ".gp 2>/dev/null").c_str());
}

void Visualization::generateAnimation(const std::vector<std::string>& vtk_files,
                                     const std::string& output_file) {
    // Would use ffmpeg or ParaView Python scripting
    // Simplified: just list the files
    
    std::string full_path = output_dir + "/" + output_file + ".pvd";
    std::ofstream pvd(full_path);
    
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd << "  <Collection>\n";
    
    for (size_t i = 0; i < vtk_files.size(); ++i) {
        pvd << "    <DataSet timestep=\"" << i << "\" file=\"" 
            << vtk_files[i] << "\"/>\n";
    }
    
    pvd << "  </Collection>\n";
    pvd << "</VTKFile>\n";
    
    pvd.close();
}

void Visualization::generateGnuplotScript(const std::string& script_content,
                                         const std::string& output_file) {
    std::string full_path = output_dir + "/" + output_file;
    std::ofstream script(full_path);
    script << script_content;
    script.close();
    
    system(("gnuplot " + full_path + " 2>/dev/null").c_str());
}

// ============================================================================
// PlotGenerator2D Methods
// ============================================================================

void PlotGenerator2D::drawField(const std::vector<std::vector<double>>& data,
                               const ColorMap& cmap, double vmin, double vmax) {
    if (data.empty() || data[0].empty()) return;
    
    int ny = data.size();
    int nx = data[0].size();
    
    for (int py = 0; py < field_height_; ++py) {
        for (int px = 0; px < field_width_; ++px) {
            int i = px * nx / field_width_;
            int j = (field_height_ - 1 - py) * ny / field_height_;
            
            i = std::min(i, nx - 1);
            j = std::min(j, ny - 1);
            
            double value = (data[j][i] - vmin) / (vmax - vmin + 1e-30);
            Color color = cmap.getColor(value);
            
            int img_y = margin_t_ + py;
            int img_x = margin_l_ + px;
            
            if (img_x < plot_width_ && img_y < plot_height_) {
                image_[img_y][img_x] = color;
            }
        }
    }
    
    // Draw border
    drawRect(margin_l_ - 1, margin_t_ - 1,
            margin_l_ + field_width_, margin_t_ + field_height_,
            Color(0, 0, 0), 2);
}

void PlotGenerator2D::drawColorbar(const ColorMap& cmap, double vmin, double vmax,
                                  const std::string& label) {
    int cbar_x = plot_width_ - margin_r_ + 20;
    int cbar_width = 30;
    int cbar_height = field_height_;
    
    // Draw colorbar gradient
    for (int py = 0; py < cbar_height; ++py) {
        double value = 1.0 - (double)py / cbar_height;
        Color color = cmap.getColor(value);
        
        for (int px = 0; px < cbar_width; ++px) {
            int img_y = margin_t_ + py;
            int img_x = cbar_x + px;
            if (img_x < plot_width_ && img_y < plot_height_) {
                image_[img_y][img_x] = color;
            }
        }
    }
    
    // Draw colorbar border
    drawRect(cbar_x - 1, margin_t_ - 1,
            cbar_x + cbar_width, margin_t_ + cbar_height,
            Color(0, 0, 0), 1);
    
    // Draw tick marks and labels
    Color text_color(0, 0, 0);
    int num_ticks = 5;
    
    for (int i = 0; i <= num_ticks; ++i) {
        int tick_y = margin_t_ + (i * cbar_height) / num_ticks;
        double value = vmax - (vmax - vmin) * i / num_ticks;
        
        // Draw tick mark
        for (int dx = 0; dx < 5; ++dx) {
            int px = cbar_x + cbar_width + dx;
            if (px < plot_width_ && tick_y >= 0 && tick_y < plot_height_) {
                image_[tick_y][px] = text_color;
            }
        }
        
        // Format value
        char value_str[32];
        if (std::abs(value) >= 1e6 || (std::abs(value) < 1e-3 && std::abs(value) > 0)) {
            snprintf(value_str, sizeof(value_str), "%.1e", value);
        } else if (std::abs(value) >= 100) {
            snprintf(value_str, sizeof(value_str), "%.0f", value);
        } else if (std::abs(value) >= 1) {
            snprintf(value_str, sizeof(value_str), "%.1f", value);
        } else {
            snprintf(value_str, sizeof(value_str), "%.3f", value);
        }
        
        // Draw label with proper spacing for larger font
        int label_x = cbar_x + cbar_width + 15;
        int label_y = tick_y - BitmapFont::getCharHeight(3) / 2;
        drawText(label_x, label_y, value_str, text_color);
    }
}

void PlotGenerator2D::drawWell(double x_pos, double y_pos, double Lx, double Ly,
                              bool is_injector, const std::string& label) {
    int img_x = margin_l_ + static_cast<int>((x_pos / Lx) * field_width_);
    int img_y = margin_t_ + field_height_ - static_cast<int>((y_pos / Ly) * field_height_);
    
    Color well_color = is_injector ? Color(0, 200, 0) : Color(200, 0, 0);
    int radius = 6;
    
    // Draw filled circle
    for (int dy = -radius; dy <= radius; ++dy) {
        for (int dx = -radius; dx <= radius; ++dx) {
            if (dx*dx + dy*dy <= radius*radius) {
                int px = img_x + dx;
                int py = img_y + dy;
                if (px >= 0 && px < plot_width_ && py >= 0 && py < plot_height_) {
                    image_[py][px] = well_color;
                }
            }
        }
    }
    
    // Draw black outline
    for (int angle = 0; angle < 360; angle += 15) {
        double rad = angle * M_PI / 180.0;
        int px = img_x + static_cast<int>(radius * cos(rad));
        int py = img_y + static_cast<int>(radius * sin(rad));
        if (px >= 0 && px < plot_width_ && py >= 0 && py < plot_height_) {
            image_[py][px] = Color(0, 0, 0);
        }
    }
    
    // Draw well label
    if (!label.empty()) {
        Color text_color(0, 0, 0);
        int label_x = img_x + radius + 5;
        int label_y = img_y - radius;
        drawText(label_x, label_y, label, text_color);
    }
}

void PlotGenerator2D::drawGrid(int num_x_lines, int num_y_lines) {
    Color grid_color(200, 200, 200);
    
    for (int i = 0; i <= num_x_lines; ++i) {
        int x = margin_l_ + (i * field_width_) / num_x_lines;
        for (int y = margin_t_; y < margin_t_ + field_height_; y += 3) {
            if (x >= 0 && x < plot_width_ && y >= 0 && y < plot_height_) {
                image_[y][x] = grid_color;
            }
        }
    }
    
    for (int j = 0; j <= num_y_lines; ++j) {
        int y = margin_t_ + (j * field_height_) / num_y_lines;
        for (int x = margin_l_; x < margin_l_ + field_width_; x += 3) {
            if (x >= 0 && x < plot_width_ && y >= 0 && y < plot_height_) {
                image_[y][x] = grid_color;
            }
        }
    }
}

void PlotGenerator2D::drawRect(int x0, int y0, int x1, int y1,
                              const Color& color, int thickness) {
    for (int t = 0; t < thickness; ++t) {
        for (int x = x0 - t; x <= x1 + t; ++x) {
            if (x >= 0 && x < plot_width_) {
                if (y0 - t >= 0 && y0 - t < plot_height_)
                    image_[y0 - t][x] = color;
                if (y1 + t >= 0 && y1 + t < plot_height_)
                    image_[y1 + t][x] = color;
            }
        }
        for (int y = y0 - t; y <= y1 + t; ++y) {
            if (y >= 0 && y < plot_height_) {
                if (x0 - t >= 0 && x0 - t < plot_width_)
                    image_[y][x0 - t] = color;
                if (x1 + t >= 0 && x1 + t < plot_width_)
                    image_[y][x1 + t] = color;
            }
        }
    }
}

void PlotGenerator2D::setTitle(const std::string& title) {
    title_ = title;
}

void PlotGenerator2D::setXLabel(const std::string& label) {
    xlabel_ = label;
}

void PlotGenerator2D::setYLabel(const std::string& label) {
    ylabel_ = label;
}

void PlotGenerator2D::drawText(int x, int y, const std::string& text, const Color& color) {
    BitmapFont::drawString(const_cast<std::vector<std::vector<Color>>&>(image_), x, y, text, color, 3);
}

void PlotGenerator2D::writePPM(const std::string& filename) const {
    // First, render title and labels on the image
    Color text_color(0, 0, 0);
    
    // Render title (centered at top)
    if (!title_.empty()) {
        int title_x = (plot_width_ - title_.length() * BitmapFont::getCharWidth(3)) / 2;
        int title_y = 30;
        const_cast<PlotGenerator2D*>(this)->drawText(title_x, title_y, title_, text_color);
    }
    
    // Render X label (centered below field)
    if (!xlabel_.empty()) {
        int xlabel_x = margin_l_ + (field_width_ - xlabel_.length() * BitmapFont::getCharWidth(3)) / 2;
        int xlabel_y = margin_t_ + field_height_ + 40;
        const_cast<PlotGenerator2D*>(this)->drawText(xlabel_x, xlabel_y, xlabel_, text_color);
    }
    
    // Render Y label (at left)
    if (!ylabel_.empty()) {
        int ylabel_x = 10;
        int ylabel_y = margin_t_ + field_height_ / 2 - ylabel_.length() * BitmapFont::getCharWidth(3) / 2;
        const_cast<PlotGenerator2D*>(this)->drawText(ylabel_x, ylabel_y, ylabel_, text_color);
    }
    
    std::ofstream file(filename, std::ios::binary);
    if (!file) return;
    
    file << "P6\n" << plot_width_ << " " << plot_height_ << "\n255\n";
    
    for (int y = 0; y < plot_height_; ++y) {
        for (int x = 0; x < plot_width_; ++x) {
            file.write(reinterpret_cast<const char*>(&image_[y][x].r), 1);
            file.write(reinterpret_cast<const char*>(&image_[y][x].g), 1);
            file.write(reinterpret_cast<const char*>(&image_[y][x].b), 1);
        }
    }
    file.close();
}

bool PlotGenerator2D::convertToPNG(const std::string& ppm_file,
                                  const std::string& png_file) const {
    std::string cmd = "command -v convert > /dev/null 2>&1 && convert " +
                     ppm_file + " " + png_file + " 2>/dev/null && rm -f " + ppm_file;
    return (system(cmd.c_str()) == 0);
}

void PlotGenerator2D::writeImage(const std::string& base_filename) const {
    std::string ppm_file = base_filename + ".ppm";
    std::string png_file = base_filename + ".png";
    
    writePPM(ppm_file);
    
    if (!convertToPNG(ppm_file, png_file)) {
        // PNG conversion failed, keep PPM
    }
}

// ============================================================================
// Visualization::plot2DField Implementation
// ============================================================================

void Visualization::plot2DField(const std::vector<std::vector<double>>& data,
                               double Lx, double Ly,
                               const std::string& title,
                               const std::string& colormap_name,
                               const std::string& filename,
                               double vmin, double vmax) {
    // Create appropriate colormap
    const ColorMap* cmap = nullptr;
    JetColorMap jet;
    ViridisColorMap viridis;
    PressureColorMap pressure;
    SubsidenceColorMap subsidence;
    
    if (colormap_name == "jet") {
        cmap = &jet;
    } else if (colormap_name == "viridis") {
        cmap = &viridis;
    } else if (colormap_name == "pressure") {
        cmap = &pressure;
    } else if (colormap_name == "subsidence") {
        cmap = &subsidence;
    } else {
        cmap = &jet;  // default
    }
    
    // Create plot
    PlotGenerator2D plot(1200, 600);
    plot.setTitle(title);
    plot.drawField(data, *cmap, vmin, vmax);
    plot.drawColorbar(*cmap, vmin, vmax);
    
    // Write output
    std::string full_path = output_dir + "/" + filename;
    plot.writeImage(full_path);
}

} // namespace FSRM
