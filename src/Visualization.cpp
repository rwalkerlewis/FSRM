#include "Visualization.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace ResSim {

Visualization::Visualization() 
    : format(OutputFormat::VTK), output_dir("output") {}

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

} // namespace ResSim
