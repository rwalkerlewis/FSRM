#include "GnuplotViz.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>

namespace FSRM {

GnuplotViz::GnuplotViz(const std::string& output_dir)
    : output_dir_(output_dir) {
    mkdir(output_dir_.c_str(), 0755);
}

GnuplotViz::~GnuplotViz() {}

std::string GnuplotViz::getPalette(const std::string& colormap) {
    if (colormap == "viridis") {
        return "defined ( 0 '#440154', 1 '#482475', 2 '#414487', 3 '#355f8d', 4 '#2a788e', 5 '#21918c', 6 '#22a884', 7 '#44bf70', 8 '#7ad151', 9 '#bddf26', 10 '#fde724' )";
    } else if (colormap == "plasma") {
        return "defined ( 0 '#0d0887', 1 '#41049d', 2 '#6a00a8', 3 '#8f0da4', 4 '#b12a90', 5 '#cc4778', 6 '#e16462', 7 '#f2844b', 8 '#fca636', 9 '#fcce25', 10 '#f0f921' )";
    } else if (colormap == "jet") {
        return "defined ( 0 '#000080', 1 '#0000ff', 2 '#0080ff', 3 '#00ffff', 4 '#80ff80', 5 '#ffff00', 6 '#ff8000', 7 '#ff0000', 8 '#800000' )";
    } else if (colormap == "coolwarm") {
        return "defined ( 0 '#3b4cc0', 1 '#7396f5', 2 '#b0d5f5', 3 '#edd1c2', 4 '#f7a789', 5 '#e36a53', 6 '#b40426' )";
    } else {
        // Default: parula-like
        return "defined ( 0 '#352a87', 1 '#0363e1', 2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e' )";
    }
}

void GnuplotViz::writeDataMatrix(const std::vector<std::vector<double>>& data,
                                 double Lx, double Ly,
                                 const std::string& filename) {
    std::ofstream file(filename);
    if (!file) return;
    
    int ny = data.size();
    int nx = data[0].size();
    double dx = Lx / nx;
    double dy = Ly / ny;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double x = i * dx + dx/2;
            double y = j * dy + dy/2;
            file << x << " " << y << " " << data[j][i] << "\n";
        }
        file << "\n";  // Blank line for gnuplot matrix
    }
    file.close();
}

void GnuplotViz::plot2DField(const std::vector<std::vector<double>>& data,
                            double Lx, double Ly,
                            const std::string& title,
                            const std::string& xlabel,
                            const std::string& ylabel,
                            const std::string& colormap,
                            const std::string& filename,
                            double vmin, double vmax) {
    // Write data file
    std::string data_file = output_dir_ + "/" + filename + ".dat";
    writeDataMatrix(data, Lx, Ly, data_file);
    
    // Auto-detect range if both are zero
    if (vmin == 0.0 && vmax == 0.0) {
        double dmin = 1e30, dmax = -1e30;
        for (const auto& row : data) {
            for (double val : row) {
                dmin = std::min(dmin, val);
                dmax = std::max(dmax, val);
            }
        }
        vmin = dmin;
        vmax = dmax;
        
        // If still the same, add small range
        if (vmin == vmax) {
            if (std::abs(vmin) < 1e-10) {
                vmin = -1e-10;
                vmax = 1e-10;
            } else {
                double delta = std::abs(vmin) * 0.1;
                vmin -= delta;
                vmax += delta;
            }
        }
    }
    
    // Create gnuplot script
    std::string script_file = output_dir_ + "/" + filename + ".gp";
    std::ofstream script(script_file);
    
    script << "set terminal pngcairo size 1400,800 enhanced font 'Arial,14'\n";
    script << "set output '" << filename << ".png'\n\n";
    
    script << "set title '" << title << "' font 'Arial,16'\n";
    script << "set xlabel '" << xlabel << "' font 'Arial,14'\n";
    script << "set ylabel '" << ylabel << "' font 'Arial,14'\n\n";
    
    script << "set pm3d map\n";
    script << "set palette " << getPalette(colormap) << "\n";
    script << "set cbrange [" << vmin << ":" << vmax << "]\n";
    script << "set cblabel offset 2\n\n";
    
    script << "set xrange [0:" << Lx << "]\n";
    script << "set yrange [0:" << Ly << "]\n\n";
    
    script << "set style data pm3d\n";
    script << "set pm3d interpolate 0,0\n\n";
    
    // Plot the field
    script << "splot '" << filename << ".dat' using 1:2:3 notitle\n";
    
    // Add wells if any
    if (!current_wells_.empty()) {
        script << "\nset object circle at ";
        for (size_t i = 0; i < current_wells_.size(); ++i) {
            const auto& w = current_wells_[i];
            if (i > 0) script << ", ";
            script << "first " << w.x << "," << w.y << " radius char 0.5 ";
            script << "fillcolor rgb '" << (w.is_injector ? "green" : "red") << "' ";
            script << "fillstyle solid border rgb 'black'";
        }
        script << "\n";
        
        for (const auto& w : current_wells_) {
            script << "set label '" << w.label << "' at " << w.x << "," << w.y 
                   << " offset char 1,0.5 font 'Arial,12'\n";
        }
    }
    
    script.close();
    
    // Execute gnuplot from output directory
    std::string cmd = "cd " + output_dir_ + " && gnuplot " + filename + ".gp 2>/dev/null";
    system(cmd.c_str());
    
    current_wells_.clear();
}

void GnuplotViz::addWells(const std::vector<WellMarker>& wells) {
    current_wells_ = wells;
}

void GnuplotViz::plotLines(const std::vector<double>& x,
                          const std::map<std::string, std::vector<double>>& y_data,
                          const std::string& title,
                          const std::string& xlabel,
                          const std::string& ylabel,
                          const std::string& filename) {
    // Write data file
    std::string data_file = output_dir_ + "/" + filename + ".dat";
    std::ofstream data(data_file);
    
    for (size_t i = 0; i < x.size(); ++i) {
        data << x[i];
        for (const auto& pair : y_data) {
            if (i < pair.second.size()) {
                data << " " << pair.second[i];
            }
        }
        data << "\n";
    }
    data.close();
    
    // Create gnuplot script
    std::string script_file = output_dir_ + "/" + filename + ".gp";
    std::ofstream script(script_file);
    
    script << "set terminal pngcairo size 1200,700 enhanced font 'Arial,14'\n";
    script << "set output '" << filename << ".png'\n\n";
    
    script << "set title '" << title << "' font 'Arial,16'\n";
    script << "set xlabel '" << xlabel << "'\n";
    script << "set ylabel '" << ylabel << "'\n";
    script << "set grid\n";
    script << "set key outside right\n\n";
    
    script << "plot ";
    int col = 2;
    for (const auto& pair : y_data) {
        if (col > 2) script << ", \\\n     ";
        script << "'" << filename << ".dat' using 1:" << col 
               << " with linespoints lw 2 pt 7 ps 0.5 title '" << pair.first << "'";
        col++;
    }
    script << "\n";
    script.close();
    
    std::string cmd = "cd " + output_dir_ + " && gnuplot " + filename + ".gp 2>/dev/null";
    system(cmd.c_str());
}

void GnuplotViz::plotMultiPanel(const std::vector<std::vector<std::vector<double>>>& data_panels,
                               const std::vector<std::string>& titles,
                               double Lx, double Ly,
                               const std::string& colormap,
                               const std::string& filename,
                               const std::vector<std::pair<double,double>>& ranges) {
    // Write data files for each panel
    std::vector<std::string> data_files;
    for (size_t i = 0; i < data_panels.size(); ++i) {
        std::string df = output_dir_ + "/" + filename + "_panel" + std::to_string(i) + ".dat";
        writeDataMatrix(data_panels[i], Lx, Ly, df);
        data_files.push_back(df);
    }
    
    // Create multiplot script
    std::string script_file = output_dir_ + "/" + filename + ".gp";
    std::ofstream script(script_file);
    
    script << "set terminal pngcairo size 1600,1400 enhanced font 'Arial,12'\n";
    script << "set output '" << filename << ".png'\n\n";
    
    script << "set multiplot layout 2,2 title '" << filename << "' font 'Arial,16'\n\n";
    script << "set pm3d map\n";
    script << "set palette " << getPalette(colormap) << "\n";
    script << "set style data pm3d\n";
    script << "set pm3d interpolate 0,0\n\n";
    
    for (size_t i = 0; i < data_panels.size() && i < 4; ++i) {
        script << "set title '" << titles[i] << "'\n";
        
        // Auto-detect range if both are zero
        double vmin = ranges[i].first;
        double vmax = ranges[i].second;
        if (vmin == 0.0 && vmax == 0.0) {
            double dmin = 1e30, dmax = -1e30;
            for (const auto& row : data_panels[i]) {
                for (double val : row) {
                    dmin = std::min(dmin, val);
                    dmax = std::max(dmax, val);
                }
            }
            vmin = dmin;
            vmax = dmax;
            
            // If still the same, add small range
            if (vmin == vmax) {
                if (std::abs(vmin) < 1e-10) {
                    vmin = -1e-10;
                    vmax = 1e-10;
                } else {
                    double delta = std::abs(vmin) * 0.1;
                    vmin -= delta;
                    vmax += delta;
                }
            }
        }
        
        script << "set cbrange [" << vmin << ":" << vmax << "]\n";
        script << "splot '" << filename << "_panel" << i << ".dat' using 1:2:3 notitle\n\n";
    }
    
    script << "unset multiplot\n";
    script.close();
    
    std::string cmd = "cd " + output_dir_ + " && gnuplot " + filename + ".gp 2>/dev/null";
    system(cmd.c_str());
}

} // namespace FSRM
