/*
 * Coulomb Stress Distribution from Slip on Parallel Faults
 * 
 * This example calculates the Coulomb stress change (ΔCFF) distribution
 * induced by slip on one fault ("source fault") onto a parallel fault 
 * ("receiver fault") in a pair of parallel faults.
 * 
 * This is a fundamental calculation in earthquake science to understand:
 * - How slip on one fault can trigger or inhibit slip on nearby faults
 * - Coulomb stress transfer and earthquake triggering
 * - Fault interactions in multi-fault systems
 * 
 * The example provides options to output results as:
 * - Snapshot images (PPM/PNG format)
 * - ASCII data files for external visualization
 * 
 * Usage:
 *   ./coulomb_stress_parallel_faults [options]
 * 
 * Options:
 *   -c <config_file>      Configuration file (optional)
 *   -slip <value>         Source fault slip in meters (default: 1.0)
 *   -sep <value>          Fault separation distance in meters (default: 5000.0)
 *   -output_format <fmt>  Output format: 'snapshot', 'ascii', or 'both' (default: both)
 *   -output_dir <dir>     Output directory (default: output)
 *   -grid_nx <n>          Grid points in x-direction (default: 100)
 *   -grid_nz <n>          Grid points in z-direction (default: 50)
 * 
 * Reference:
 *   - King et al. (1994), JGR: Static stress changes and the triggering of earthquakes
 *   - Okada (1992), BSSA: Internal deformation due to shear and tensile faults
 */

#include "FaultModel.hpp"
#include "Visualization.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>

using namespace FSRM;

// =============================================================================
// Output Format Options
// =============================================================================
enum class OutputFormat {
    SNAPSHOT,       // PNG/PPM image output
    ASCII,          // ASCII text files
    BOTH            // Both snapshot and ASCII
};

// =============================================================================
// Helper: Create output directory
// =============================================================================
void createOutputDirectory(const std::string& dir) {
    #ifdef _WIN32
        _mkdir(dir.c_str());
    #else
        mkdir(dir.c_str(), 0755);
    #endif
}

// =============================================================================
// Helper: Parse command line options
// =============================================================================
struct Options {
    double slip = 1.0;              // Source fault slip (m)
    double separation = 5000.0;     // Fault separation distance (m)
    double fault_length = 10000.0;  // Fault length (m)
    double fault_width = 5000.0;    // Fault width (m) - down-dip
    double strike = 0.0;            // Fault strike (degrees from North)
    double dip = 90.0;              // Fault dip (degrees)
    double friction = 0.4;          // Effective friction coefficient
    double shear_modulus = 30.0e9;  // Shear modulus (Pa)
    double poisson_ratio = 0.25;    // Poisson's ratio
    
    int grid_nx = 100;              // Grid points in along-strike direction
    int grid_nz = 50;               // Grid points in depth direction
    double grid_Lx = 30000.0;       // Grid extent in x (m)
    double grid_Lz = 15000.0;       // Grid extent in z (m)
    
    OutputFormat output_format = OutputFormat::BOTH;
    std::string output_dir = "output";
    std::string config_file = "";
};

Options parseOptions(int argc, char** argv) {
    Options opts;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-c" && i + 1 < argc) {
            opts.config_file = argv[++i];
        } else if (arg == "-slip" && i + 1 < argc) {
            opts.slip = std::stod(argv[++i]);
        } else if (arg == "-sep" && i + 1 < argc) {
            opts.separation = std::stod(argv[++i]);
        } else if (arg == "-length" && i + 1 < argc) {
            opts.fault_length = std::stod(argv[++i]);
        } else if (arg == "-width" && i + 1 < argc) {
            opts.fault_width = std::stod(argv[++i]);
        } else if (arg == "-strike" && i + 1 < argc) {
            opts.strike = std::stod(argv[++i]);
        } else if (arg == "-dip" && i + 1 < argc) {
            opts.dip = std::stod(argv[++i]);
        } else if (arg == "-friction" && i + 1 < argc) {
            opts.friction = std::stod(argv[++i]);
        } else if (arg == "-shear_modulus" && i + 1 < argc) {
            opts.shear_modulus = std::stod(argv[++i]);
        } else if (arg == "-grid_nx" && i + 1 < argc) {
            opts.grid_nx = std::stoi(argv[++i]);
        } else if (arg == "-grid_nz" && i + 1 < argc) {
            opts.grid_nz = std::stoi(argv[++i]);
        } else if (arg == "-grid_Lx" && i + 1 < argc) {
            opts.grid_Lx = std::stod(argv[++i]);
        } else if (arg == "-grid_Lz" && i + 1 < argc) {
            opts.grid_Lz = std::stod(argv[++i]);
        } else if (arg == "-output_format" && i + 1 < argc) {
            std::string fmt = argv[++i];
            if (fmt == "snapshot" || fmt == "SNAPSHOT") {
                opts.output_format = OutputFormat::SNAPSHOT;
            } else if (fmt == "ascii" || fmt == "ASCII") {
                opts.output_format = OutputFormat::ASCII;
            } else {
                opts.output_format = OutputFormat::BOTH;
            }
        } else if (arg == "-output_dir" && i + 1 < argc) {
            opts.output_dir = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Coulomb Stress Distribution from Parallel Faults\n\n";
            std::cout << "Usage: " << argv[0] << " [options]\n\n";
            std::cout << "Options:\n";
            std::cout << "  -c <file>           Configuration file\n";
            std::cout << "  -slip <value>       Source fault slip (m, default: 1.0)\n";
            std::cout << "  -sep <value>        Fault separation (m, default: 5000.0)\n";
            std::cout << "  -length <value>     Fault length (m, default: 10000.0)\n";
            std::cout << "  -width <value>      Fault width (m, default: 5000.0)\n";
            std::cout << "  -strike <value>     Fault strike (deg, default: 0.0)\n";
            std::cout << "  -dip <value>        Fault dip (deg, default: 90.0)\n";
            std::cout << "  -friction <value>   Effective friction (default: 0.4)\n";
            std::cout << "  -shear_modulus <v>  Shear modulus (Pa, default: 30e9)\n";
            std::cout << "  -grid_nx <n>        Grid points in x (default: 100)\n";
            std::cout << "  -grid_nz <n>        Grid points in z (default: 50)\n";
            std::cout << "  -grid_Lx <value>    Grid extent x (m, default: 30000)\n";
            std::cout << "  -grid_Lz <value>    Grid extent z (m, default: 15000)\n";
            std::cout << "  -output_format <f>  'snapshot', 'ascii', or 'both'\n";
            std::cout << "  -output_dir <dir>   Output directory (default: output)\n";
            std::cout << "  -h, --help          Show this help message\n";
            exit(0);
        }
    }
    
    return opts;
}

// =============================================================================
// Compute Okada (1992) stress changes from slip on a rectangular fault
// Simplified implementation for strike-slip faulting
// =============================================================================
struct StressTensor {
    double sxx, syy, szz, sxy, sxz, syz;
};

StressTensor computeOkadaStress(
    double x, double y, double z,           // Observation point
    double fault_x, double fault_y, double fault_depth,  // Fault center
    double strike_rad, double dip_rad,      // Fault orientation
    double length, double width,            // Fault dimensions
    double slip,                            // Slip amount
    double mu, double nu)                   // Elastic parameters
{
    StressTensor stress = {0, 0, 0, 0, 0, 0};
    
    // Transform to fault-centered coordinates
    double cos_s = std::cos(strike_rad);
    double sin_s = std::sin(strike_rad);
    
    double dx = x - fault_x;
    double dy = y - fault_y;
    double dz = z - fault_depth;
    
    // Rotate to fault-aligned coordinates
    double xp = dx * sin_s + dy * cos_s;     // Along-strike
    double yp = -dx * cos_s + dy * sin_s;    // Fault-perpendicular
    double zp = dz;                           // Depth
    
    // Distance from fault plane
    double r = std::sqrt(xp * xp + yp * yp + zp * zp);
    if (r < 1.0) r = 1.0;  // Avoid singularity
    
    // Simplified stress calculation (leading-order behavior)
    // For a vertical strike-slip fault, the stress pattern follows:
    // - Positive CFF in lobes at ~45° from fault tips
    // - Negative CFF ("stress shadow") perpendicular to fault
    
    double area = length * width;
    double M0 = mu * area * slip;  // Seismic moment
    
    // Characteristic stress scale
    double stress_scale = M0 / (4.0 * M_PI * r * r * r);
    
    // Angular factors for strike-slip fault
    double theta = std::atan2(yp, xp);
    double phi = std::atan2(zp, std::sqrt(xp * xp + yp * yp));
    
    // Strike-slip stress pattern (simplified)
    // σ_xy dominates, with characteristic lobes
    double sin2theta = std::sin(2.0 * theta);
    double cos2theta = std::cos(2.0 * theta);
    double cos_phi = std::cos(phi);
    
    // Stress components (approximate Okada)
    stress.sxx = -stress_scale * sin2theta * cos_phi;
    stress.syy = stress_scale * sin2theta * cos_phi;
    stress.sxy = stress_scale * cos2theta * cos_phi;
    stress.szz = 0.0;  // Simplified: neglect vertical stress
    stress.sxz = 0.0;
    stress.syz = 0.0;
    
    // Add depth decay factor
    double depth_factor = std::exp(-std::abs(zp) / (2.0 * width));
    stress.sxx *= depth_factor;
    stress.syy *= depth_factor;
    stress.sxy *= depth_factor;
    
    return stress;
}

// =============================================================================
// Compute Coulomb stress change on receiver fault
// ΔCFF = Δτ + μ' × Δσn
// where μ' is effective friction coefficient
// =============================================================================
double computeCoulombStressChange(
    const StressTensor& stress,
    double receiver_strike_rad,
    double receiver_dip_rad,
    double friction)
{
    // Fault normal and slip direction vectors
    double cos_s = std::cos(receiver_strike_rad);
    double sin_s = std::sin(receiver_strike_rad);
    double cos_d = std::cos(receiver_dip_rad);
    double sin_d = std::sin(receiver_dip_rad);
    
    // Normal vector to receiver fault (for vertical fault: points in y-direction)
    double nx = cos_s * sin_d;
    double ny = -sin_s * sin_d;
    double nz = cos_d;
    
    // Slip direction (strike direction for strike-slip)
    double sx = sin_s;
    double sy = cos_s;
    double sz = 0.0;
    
    // Compute traction on receiver fault: t_i = σ_ij * n_j
    double tx = stress.sxx * nx + stress.sxy * ny + stress.sxz * nz;
    double ty = stress.sxy * nx + stress.syy * ny + stress.syz * nz;
    double tz = stress.sxz * nx + stress.syz * ny + stress.szz * nz;
    
    // Normal stress change (tension positive)
    double delta_sigma_n = tx * nx + ty * ny + tz * nz;
    
    // Shear stress change (in slip direction)
    double delta_tau = tx * sx + ty * sy + tz * sz;
    
    // Coulomb stress change
    // Positive ΔCFF promotes failure, negative inhibits
    double delta_CFF = delta_tau + friction * delta_sigma_n;
    
    return delta_CFF;
}

// =============================================================================
// Write ASCII output file
// =============================================================================
void writeASCIIOutput(
    const std::string& filename,
    const std::vector<std::vector<double>>& data,
    const std::vector<double>& x_coords,
    const std::vector<double>& z_coords,
    const std::string& description)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << filename << std::endl;
        return;
    }
    
    file << "# " << description << "\n";
    file << "# Format: x (m), z (m), ΔCFF (Pa)\n";
    file << "# Columns: " << data[0].size() << ", Rows: " << data.size() << "\n";
    file << "#\n";
    
    // Write header with x coordinates
    file << "# x\\z";
    for (size_t j = 0; j < z_coords.size(); ++j) {
        file << "\t" << std::scientific << std::setprecision(4) << z_coords[j];
    }
    file << "\n";
    
    // Write data
    for (size_t i = 0; i < data.size(); ++i) {
        file << std::scientific << std::setprecision(4) << x_coords[i];
        for (size_t j = 0; j < data[i].size(); ++j) {
            file << "\t" << std::scientific << std::setprecision(6) << data[i][j];
        }
        file << "\n";
    }
    
    file.close();
    std::cout << "  Written: " << filename << std::endl;
}

// =============================================================================
// Write summary file with parameters
// =============================================================================
void writeSummaryFile(
    const std::string& filename,
    const Options& opts,
    double max_cff, double min_cff)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << filename << std::endl;
        return;
    }
    
    file << "# Coulomb Stress Distribution Analysis Summary\n";
    file << "#\n";
    file << "# Source Fault Parameters:\n";
    file << "slip = " << opts.slip << " m\n";
    file << "length = " << opts.fault_length << " m\n";
    file << "width = " << opts.fault_width << " m\n";
    file << "strike = " << opts.strike << " degrees\n";
    file << "dip = " << opts.dip << " degrees\n";
    file << "#\n";
    file << "# Receiver Fault Parameters:\n";
    file << "separation = " << opts.separation << " m (parallel offset)\n";
    file << "friction = " << opts.friction << "\n";
    file << "#\n";
    file << "# Material Properties:\n";
    file << "shear_modulus = " << std::scientific << opts.shear_modulus << " Pa\n";
    file << "poisson_ratio = " << opts.poisson_ratio << "\n";
    file << "#\n";
    file << "# Grid Parameters:\n";
    file << "grid_nx = " << opts.grid_nx << "\n";
    file << "grid_nz = " << opts.grid_nz << "\n";
    file << "grid_Lx = " << opts.grid_Lx << " m\n";
    file << "grid_Lz = " << opts.grid_Lz << " m\n";
    file << "#\n";
    file << "# Results:\n";
    file << "max_delta_CFF = " << std::scientific << max_cff << " Pa\n";
    file << "min_delta_CFF = " << std::scientific << min_cff << " Pa\n";
    file << "#\n";
    file << "# Interpretation:\n";
    file << "# Positive ΔCFF: promotes failure (brought closer to failure)\n";
    file << "# Negative ΔCFF: inhibits failure (stress shadow)\n";
    
    file.close();
    std::cout << "  Written: " << filename << std::endl;
}

// =============================================================================
// Write snapshot image (PPM format with optional PNG conversion)
// =============================================================================
void writeSnapshotImage(
    const std::string& filename,
    const std::vector<std::vector<double>>& data,
    double vmin, double vmax,
    const std::string& title,
    double Lx, double Lz,
    double source_x, double receiver_x)
{
    // Use the built-in visualization class
    int plot_width = 800;
    int plot_height = 500;
    int margin_left = 100;
    int margin_right = 150;
    int margin_top = 80;
    int margin_bottom = 60;
    
    int field_width = plot_width - margin_left - margin_right;
    int field_height = plot_height - margin_top - margin_bottom;
    
    // Create image buffer
    std::vector<std::vector<Color>> image(plot_height, 
        std::vector<Color>(plot_width, Color(255, 255, 255)));
    
    // Create colormap (blue-white-red for Coulomb stress)
    auto getColor = [](double val, double vmin, double vmax) -> Color {
        double t = (val - vmin) / (vmax - vmin);
        t = std::max(0.0, std::min(1.0, t));
        
        // Blue (negative) - White (zero) - Red (positive)
        if (t < 0.5) {
            // Blue to white
            double s = t * 2.0;
            return Color(
                static_cast<unsigned char>(s * 255),
                static_cast<unsigned char>(s * 255),
                255
            );
        } else {
            // White to red
            double s = (t - 0.5) * 2.0;
            return Color(
                255,
                static_cast<unsigned char>((1.0 - s) * 255),
                static_cast<unsigned char>((1.0 - s) * 255)
            );
        }
    };
    
    // Draw field data
    int nx = data.size();
    int nz = data[0].size();
    
    for (int iy = 0; iy < field_height; ++iy) {
        for (int ix = 0; ix < field_width; ++ix) {
            // Map pixel to data indices
            int di = static_cast<int>((ix * nx) / field_width);
            int dj = static_cast<int>((iy * nz) / field_height);
            
            di = std::min(di, nx - 1);
            dj = std::min(dj, nz - 1);
            
            double val = data[di][dj];
            Color c = getColor(val, vmin, vmax);
            
            int px = margin_left + ix;
            int py = margin_top + iy;
            
            if (px >= 0 && px < plot_width && py >= 0 && py < plot_height) {
                image[py][px] = c;
            }
        }
    }
    
    // Draw colorbar
    int cb_x = plot_width - margin_right + 20;
    int cb_width = 20;
    int cb_height = field_height;
    
    for (int iy = 0; iy < cb_height; ++iy) {
        double t = 1.0 - static_cast<double>(iy) / cb_height;
        double val = vmin + t * (vmax - vmin);
        Color c = getColor(val, vmin, vmax);
        
        for (int ix = 0; ix < cb_width; ++ix) {
            int px = cb_x + ix;
            int py = margin_top + iy;
            if (px >= 0 && px < plot_width && py >= 0 && py < plot_height) {
                image[py][px] = c;
            }
        }
    }
    
    // Draw fault locations as vertical lines
    // Source fault (black solid line)
    int src_px = margin_left + static_cast<int>((source_x / Lx + 0.5) * field_width);
    if (src_px >= margin_left && src_px < margin_left + field_width) {
        for (int iy = margin_top; iy < margin_top + field_height; ++iy) {
            image[iy][src_px] = Color(0, 0, 0);
            if (src_px + 1 < plot_width) image[iy][src_px + 1] = Color(0, 0, 0);
        }
    }
    
    // Receiver fault (dashed gray line)
    int rcv_px = margin_left + static_cast<int>((receiver_x / Lx + 0.5) * field_width);
    if (rcv_px >= margin_left && rcv_px < margin_left + field_width) {
        for (int iy = margin_top; iy < margin_top + field_height; iy += 4) {
            if (iy + 2 < margin_top + field_height) {
                image[iy][rcv_px] = Color(100, 100, 100);
                image[iy + 1][rcv_px] = Color(100, 100, 100);
            }
        }
    }
    
    // Write PPM file
    std::string ppm_filename = filename + ".ppm";
    std::ofstream ppm(ppm_filename, std::ios::binary);
    if (ppm.is_open()) {
        ppm << "P6\n" << plot_width << " " << plot_height << "\n255\n";
        for (int y = 0; y < plot_height; ++y) {
            for (int x = 0; x < plot_width; ++x) {
                ppm.put(image[y][x].r);
                ppm.put(image[y][x].g);
                ppm.put(image[y][x].b);
            }
        }
        ppm.close();
        std::cout << "  Written: " << ppm_filename << std::endl;
        
        // Try to convert to PNG
        std::string png_filename = filename + ".png";
        std::string cmd = "convert " + ppm_filename + " " + png_filename + " 2>/dev/null";
        int ret = system(cmd.c_str());
        if (ret == 0) {
            std::cout << "  Converted to: " << png_filename << std::endl;
        }
    }
}

// =============================================================================
// Main function
// =============================================================================
int main(int argc, char** argv) {
    std::cout << "============================================================\n";
    std::cout << "  Coulomb Stress Distribution from Parallel Faults\n";
    std::cout << "============================================================\n\n";
    
    // Parse command line options
    Options opts = parseOptions(argc, argv);
    
    // Create output directory
    createOutputDirectory(opts.output_dir);
    
    // Print configuration
    std::cout << "Configuration:\n";
    std::cout << "  Source fault slip:    " << opts.slip << " m\n";
    std::cout << "  Fault separation:     " << opts.separation << " m\n";
    std::cout << "  Fault length:         " << opts.fault_length << " m\n";
    std::cout << "  Fault width:          " << opts.fault_width << " m\n";
    std::cout << "  Strike:               " << opts.strike << " deg\n";
    std::cout << "  Dip:                  " << opts.dip << " deg\n";
    std::cout << "  Friction coefficient: " << opts.friction << "\n";
    std::cout << "  Shear modulus:        " << std::scientific << opts.shear_modulus << " Pa\n";
    std::cout << "  Grid size:            " << opts.grid_nx << " x " << opts.grid_nz << "\n";
    std::cout << "  Output format:        ";
    switch (opts.output_format) {
        case OutputFormat::SNAPSHOT: std::cout << "Snapshot images\n"; break;
        case OutputFormat::ASCII: std::cout << "ASCII files\n"; break;
        case OutputFormat::BOTH: std::cout << "Both snapshot and ASCII\n"; break;
    }
    std::cout << "  Output directory:     " << opts.output_dir << "\n\n";
    
    // Convert angles to radians
    double strike_rad = opts.strike * M_PI / 180.0;
    double dip_rad = opts.dip * M_PI / 180.0;
    
    // Define fault positions
    // Source fault at center (x = 0)
    double source_x = 0.0;
    double source_y = 0.0;
    double source_depth = opts.fault_width / 2.0;  // Center depth
    
    // Receiver fault offset by separation distance
    double receiver_x = opts.separation;
    
    // Create computation grid on receiver fault plane
    std::vector<double> x_coords(opts.grid_nx);
    std::vector<double> z_coords(opts.grid_nz);
    
    double dx = opts.grid_Lx / (opts.grid_nx - 1);
    double dz = opts.grid_Lz / (opts.grid_nz - 1);
    
    for (int i = 0; i < opts.grid_nx; ++i) {
        x_coords[i] = -opts.grid_Lx / 2.0 + i * dx;  // Centered on fault
    }
    for (int j = 0; j < opts.grid_nz; ++j) {
        z_coords[j] = j * dz;  // Surface to depth
    }
    
    // Compute Coulomb stress change at each grid point
    std::cout << "Computing Coulomb stress distribution...\n";
    
    std::vector<std::vector<double>> delta_CFF(opts.grid_nx, 
        std::vector<double>(opts.grid_nz, 0.0));
    
    double max_cff = -1e30;
    double min_cff = 1e30;
    
    for (int i = 0; i < opts.grid_nx; ++i) {
        for (int j = 0; j < opts.grid_nz; ++j) {
            // Observation point on receiver fault
            double obs_x = x_coords[i];
            double obs_y = receiver_x;  // On receiver fault plane
            double obs_z = z_coords[j];
            
            // Compute stress from source fault slip
            StressTensor stress = computeOkadaStress(
                obs_x, obs_y, obs_z,
                source_x, source_y, source_depth,
                strike_rad, dip_rad,
                opts.fault_length, opts.fault_width,
                opts.slip,
                opts.shear_modulus, opts.poisson_ratio
            );
            
            // Compute Coulomb stress change on receiver
            double cff = computeCoulombStressChange(
                stress, strike_rad, dip_rad, opts.friction
            );
            
            delta_CFF[i][j] = cff;
            
            max_cff = std::max(max_cff, cff);
            min_cff = std::min(min_cff, cff);
        }
    }
    
    std::cout << "  Maximum ΔCFF: " << std::scientific << max_cff << " Pa\n";
    std::cout << "  Minimum ΔCFF: " << std::scientific << min_cff << " Pa\n\n";
    
    // Determine symmetric color range
    double vmax = std::max(std::abs(max_cff), std::abs(min_cff));
    double vmin = -vmax;
    
    // Output results
    std::cout << "Writing output files...\n";
    
    // Write summary file (always)
    writeSummaryFile(
        opts.output_dir + "/coulomb_stress_summary.txt",
        opts, max_cff, min_cff
    );
    
    // Write ASCII output
    if (opts.output_format == OutputFormat::ASCII || 
        opts.output_format == OutputFormat::BOTH) {
        
        writeASCIIOutput(
            opts.output_dir + "/coulomb_stress_delta_CFF.dat",
            delta_CFF, x_coords, z_coords,
            "Coulomb stress change (Pa) on receiver fault"
        );
        
        // Write coordinate files for external plotting
        std::ofstream x_file(opts.output_dir + "/coulomb_stress_x_coords.dat");
        x_file << "# X coordinates (m) along strike\n";
        for (int i = 0; i < opts.grid_nx; ++i) {
            x_file << x_coords[i] << "\n";
        }
        x_file.close();
        std::cout << "  Written: " << opts.output_dir << "/coulomb_stress_x_coords.dat\n";
        
        std::ofstream z_file(opts.output_dir + "/coulomb_stress_z_coords.dat");
        z_file << "# Z coordinates (m) depth\n";
        for (int j = 0; j < opts.grid_nz; ++j) {
            z_file << z_coords[j] << "\n";
        }
        z_file.close();
        std::cout << "  Written: " << opts.output_dir << "/coulomb_stress_z_coords.dat\n";
        
        // Write gnuplot script for visualization
        std::ofstream gp(opts.output_dir + "/plot_coulomb_stress.gp");
        gp << "# Gnuplot script for Coulomb stress visualization\n";
        gp << "set terminal pngcairo size 1000,600 enhanced font 'Arial,12'\n";
        gp << "set output 'coulomb_stress_gnuplot.png'\n";
        gp << "\n";
        gp << "set title 'Coulomb Stress Change on Receiver Fault'\n";
        gp << "set xlabel 'Along-strike distance (m)'\n";
        gp << "set ylabel 'Depth (m)'\n";
        gp << "set cblabel '{/Symbol D}CFF (Pa)'\n";
        gp << "\n";
        gp << "set yrange [" << opts.grid_Lz << ":0]\n";
        gp << "set palette defined (-1 'blue', 0 'white', 1 'red')\n";
        gp << "set cbrange [" << vmin << ":" << vmax << "]\n";
        gp << "\n";
        gp << "# Source fault location (x=0)\n";
        gp << "set arrow from 0,0 to 0," << opts.grid_Lz << " nohead lw 2 lc 'black'\n";
        gp << "\n";
        gp << "plot 'coulomb_stress_delta_CFF.dat' matrix nonuniform with image notitle\n";
        gp.close();
        std::cout << "  Written: " << opts.output_dir << "/plot_coulomb_stress.gp\n";
    }
    
    // Write snapshot image
    if (opts.output_format == OutputFormat::SNAPSHOT || 
        opts.output_format == OutputFormat::BOTH) {
        
        std::stringstream title;
        title << "Coulomb Stress Change (Slip=" << opts.slip << "m, Sep=" 
              << opts.separation << "m)";
        
        writeSnapshotImage(
            opts.output_dir + "/coulomb_stress_snapshot",
            delta_CFF, vmin, vmax,
            title.str(),
            opts.grid_Lx, opts.grid_Lz,
            source_x, receiver_x
        );
    }
    
    std::cout << "\n============================================================\n";
    std::cout << "  Coulomb Stress Analysis Complete\n";
    std::cout << "============================================================\n\n";
    std::cout << "Results:\n";
    std::cout << "  - Maximum positive ΔCFF (promotes failure): " 
              << std::scientific << max_cff << " Pa\n";
    std::cout << "  - Maximum negative ΔCFF (stress shadow): " 
              << std::scientific << min_cff << " Pa\n";
    std::cout << "\n";
    std::cout << "Interpretation:\n";
    std::cout << "  - Positive ΔCFF regions: fault is brought closer to failure\n";
    std::cout << "  - Negative ΔCFF regions: fault is moved away from failure (stress shadow)\n";
    std::cout << "  - Pattern shows characteristic lobes of Coulomb stress transfer\n";
    std::cout << "\n";
    std::cout << "Output files in: " << opts.output_dir << "/\n";
    
    return 0;
}
