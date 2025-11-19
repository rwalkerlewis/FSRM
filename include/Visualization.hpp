#ifndef VISUALIZATION_HPP
#define VISUALIZATION_HPP

#include "ReservoirSim.hpp"
#include <string>
#include <vector>
#include <map>

namespace ResSim {

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
