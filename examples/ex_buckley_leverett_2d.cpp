#include "Simulator.hpp"
#include "PhysicsKernel.hpp"
#include "Visualization.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sys/stat.h>

using namespace FSRM;

// Buckley-Leverett fractional flow function
double fractionalFlow(double Sw, double mu_w, double mu_o, double krw_max = 1.0, double kro_max = 1.0) {
    // Corey relative permeability model
    double n_w = 2.0;  // Water exponent
    double n_o = 2.0;  // Oil exponent
    double Sw_c = 0.2; // Connate water saturation
    double Sor = 0.2;  // Residual oil saturation
    
    // Normalize saturation
    double Sw_norm = (Sw - Sw_c) / (1.0 - Sw_c - Sor);
    Sw_norm = std::max(0.0, std::min(1.0, Sw_norm));
    
    // Relative permeabilities
    double krw = krw_max * std::pow(Sw_norm, n_w);
    double kro = kro_max * std::pow(1.0 - Sw_norm, n_o);
    
    // Fractional flow
    double lambda_w = krw / mu_w;
    double lambda_o = kro / mu_o;
    
    if (lambda_w + lambda_o < 1e-12) return 0.0;
    
    return lambda_w / (lambda_w + lambda_o);
}

// Simple PPM image writer (Portable Pixmap Format)
void writePPM(const std::string& filename, const std::vector<std::vector<double>>& data,
              int width, int height) {
    std::ofstream file(filename, std::ios::binary);
    
    // PPM header
    file << "P6\n" << width << " " << height << "\n255\n";
    
    // Color map: blue (low) -> white (mid) -> red (high)
    auto getColor = [](double value) -> std::array<unsigned char, 3> {
        // Value should be 0-1, map to color
        value = std::max(0.0, std::min(1.0, value));
        
        unsigned char r, g, b;
        if (value < 0.5) {
            // Blue to white
            double t = value * 2.0;
            r = static_cast<unsigned char>(255 * t);
            g = static_cast<unsigned char>(255 * t);
            b = 255;
        } else {
            // White to red
            double t = (value - 0.5) * 2.0;
            r = 255;
            g = static_cast<unsigned char>(255 * (1.0 - t));
            b = static_cast<unsigned char>(255 * (1.0 - t));
        }
        return {r, g, b};
    };
    
    // Write pixels (flip vertically for correct orientation)
    for (int j = height - 1; j >= 0; --j) {
        for (int i = 0; i < width; ++i) {
            auto color = getColor(data[j][i]);
            file.write(reinterpret_cast<char*>(color.data()), 3);
        }
    }
    file.close();
}

// Convert PPM to PNG using ImageMagick if available
void convertPPMtoPNG(const std::string& ppm_file, const std::string& png_file) {
    std::string cmd = "command -v convert > /dev/null 2>&1 && convert " + 
                     ppm_file + " " + png_file + " 2>/dev/null";
    int result = system(cmd.c_str());
    if (result == 0) {
        // Remove temporary PPM file
        std::remove(ppm_file.c_str());
    }
}

// Generate saturation plot with colorbar and labels
void generateSaturationPlot(const std::string& output_dir, int step,
                           const std::vector<std::vector<double>>& Sw,
                           int nx, int ny, double Lx, double Ly, double time_days) {
    const int plot_width = 1200;
    const int plot_height = 400;
    const int margin_left = 80;
    const int margin_right = 150;
    const int margin_top = 80;
    const int margin_bottom = 60;
    
    int field_width = plot_width - margin_left - margin_right;
    int field_height = plot_height - margin_top - margin_bottom;
    
    // Create image buffer (white background)
    std::vector<std::vector<std::array<unsigned char, 3>>> image(
        plot_height, std::vector<std::array<unsigned char, 3>>(plot_width, {255, 255, 255}));
    
    // Color mapping function
    auto getSaturationColor = [](double Sw) -> std::array<unsigned char, 3> {
        Sw = std::max(0.0, std::min(1.0, Sw));
        
        // Blue (water) to red (oil depleted)
        unsigned char r = static_cast<unsigned char>(255 * Sw);
        unsigned char g = static_cast<unsigned char>(100 * (1.0 - Sw));
        unsigned char b = static_cast<unsigned char>(255 * (1.0 - Sw));
        return {r, g, b};
    };
    
    // Draw saturation field
    for (int py = 0; py < field_height; ++py) {
        for (int px = 0; px < field_width; ++px) {
            int i = px * nx / field_width;
            int j = py * ny / field_height;
            
            int img_y = margin_top + py;
            int img_x = margin_left + px;
            
            image[img_y][img_x] = getSaturationColor(Sw[j][i]);
        }
    }
    
    // Draw colorbar
    int cbar_x = plot_width - margin_right + 20;
    int cbar_width = 30;
    int cbar_height = field_height;
    
    for (int py = 0; py < cbar_height; ++py) {
        double value = 1.0 - (double)py / cbar_height;
        auto color = getSaturationColor(value);
        
        for (int px = 0; px < cbar_width; ++px) {
            int img_y = margin_top + py;
            int img_x = cbar_x + px;
            if (img_x < plot_width) {
                image[img_y][img_x] = color;
            }
        }
    }
    
    // Draw border around field
    auto drawRect = [&](int x0, int y0, int x1, int y1, std::array<unsigned char, 3> color) {
        for (int x = x0; x <= x1; ++x) {
            if (x >= 0 && x < plot_width && y0 >= 0 && y0 < plot_height)
                image[y0][x] = color;
            if (x >= 0 && x < plot_width && y1 >= 0 && y1 < plot_height)
                image[y1][x] = color;
        }
        for (int y = y0; y <= y1; ++y) {
            if (y >= 0 && y < plot_height && x0 >= 0 && x0 < plot_width)
                image[y][x0] = color;
            if (y >= 0 && y < plot_height && x1 >= 0 && x1 < plot_width)
                image[y][x1] = color;
        }
    };
    
    drawRect(margin_left - 1, margin_top - 1, 
             margin_left + field_width, margin_top + field_height, {0, 0, 0});
    drawRect(cbar_x - 1, margin_top - 1,
             cbar_x + cbar_width, margin_top + cbar_height, {0, 0, 0});
    
    // Write to PPM file
    char filename_ppm[512];
    snprintf(filename_ppm, sizeof(filename_ppm), "%s/saturation_%04d.ppm", 
             output_dir.c_str(), step);
    
    std::ofstream file(filename_ppm, std::ios::binary);
    file << "P6\n" << plot_width << " " << plot_height << "\n255\n";
    
    for (int y = 0; y < plot_height; ++y) {
        for (int x = 0; x < plot_width; ++x) {
            file.write(reinterpret_cast<char*>(image[y][x].data()), 3);
        }
    }
    file.close();
    
    // Try to convert to PNG
    char filename_png[512];
    snprintf(filename_png, sizeof(filename_png), "%s/saturation_%04d.png", 
             output_dir.c_str(), step);
    convertPPMtoPNG(filename_ppm, filename_png);
    
    // Also create text overlay file with metadata
    char metadata_file[512];
    snprintf(metadata_file, sizeof(metadata_file), "%s/info_%04d.txt", 
             output_dir.c_str(), step);
    std::ofstream info(metadata_file);
    info << "Step: " << step << "\n";
    info << "Time: " << std::fixed << std::setprecision(2) << time_days << " days\n";
    info << "Domain: " << Lx << " m x " << Ly << " m\n";
    info << "Grid: " << nx << " x " << ny << " cells\n";
    
    double avg_sw = 0.0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            avg_sw += Sw[j][i];
        }
    }
    avg_sw /= (nx * ny);
    info << "Average Sw: " << std::setprecision(4) << avg_sw << "\n";
    info.close();
}

// Write data to CSV for plotting
void writeCSV(const std::string& filename, const std::vector<double>& x, 
              const std::vector<double>& y, const std::vector<std::vector<double>>& data,
              const std::string& data_name) {
    std::ofstream file(filename);
    file << "x,y," << data_name << "\n";
    
    for (size_t j = 0; j < y.size(); ++j) {
        for (size_t i = 0; i < x.size(); ++i) {
            file << x[i] << "," << y[j] << "," << data[j][i] << "\n";
        }
    }
    file.close();
}

// Generate Python plotting script
void generatePlotScript(const std::string& output_dir, int num_steps) {
    std::ofstream script(output_dir + "/plot_results.py");
    
    script << "#!/usr/bin/env python3\n";
    script << "import numpy as np\n";
    script << "import matplotlib.pyplot as plt\n";
    script << "import pandas as pd\n";
    script << "import os\n\n";
    
    script << "# Set style\n";
    script << "plt.style.use('seaborn-v0_8-darkgrid')\n";
    script << "plt.rcParams['figure.figsize'] = (12, 5)\n\n";
    
    script << "output_dir = '" << output_dir << "'\n\n";
    
    script << "# Get list of timestep files\n";
    script << "timesteps = []\n";
    script << "for i in range(" << num_steps + 1 << "):\n";
    script << "    filename = f'saturation_step_{i:04d}.csv'\n";
    script << "    if os.path.exists(os.path.join(output_dir, filename)):\n";
    script << "        timesteps.append(i)\n\n";
    
    script << "print(f'Found {len(timesteps)} timesteps')\n\n";
    
    script << "# Plot each timestep\n";
    script << "for step in timesteps:\n";
    script << "    filename = f'saturation_step_{step:04d}.csv'\n";
    script << "    df = pd.read_csv(os.path.join(output_dir, filename))\n";
    script << "    \n";
    script << "    # Pivot data for 2D plotting\n";
    script << "    x_unique = sorted(df['x'].unique())\n";
    script << "    y_unique = sorted(df['y'].unique())\n";
    script << "    \n";
    script << "    nx, ny = len(x_unique), len(y_unique)\n";
    script << "    saturation = df['saturation'].values.reshape(ny, nx)\n";
    script << "    \n";
    script << "    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))\n";
    script << "    \n";
    script << "    # 2D saturation map\n";
    script << "    im = ax1.imshow(saturation, extent=[min(x_unique), max(x_unique), \n";
    script << "                                         min(y_unique), max(y_unique)],\n";
    script << "                    origin='lower', cmap='RdYlBu_r', vmin=0, vmax=1, aspect='auto')\n";
    script << "    ax1.set_xlabel('X (m)', fontsize=12)\n";
    script << "    ax1.set_ylabel('Y (m)', fontsize=12)\n";
    script << "    ax1.set_title(f'Water Saturation - Step {step}', fontsize=14, fontweight='bold')\n";
    script << "    cbar = plt.colorbar(im, ax=ax1)\n";
    script << "    cbar.set_label('Water Saturation', fontsize=11)\n";
    script << "    ax1.grid(True, alpha=0.3)\n";
    script << "    \n";
    script << "    # 1D profile at mid-height\n";
    script << "    mid_y = ny // 2\n";
    script << "    ax2.plot(x_unique, saturation[mid_y, :], 'b-', linewidth=2, label='Numerical')\n";
    script << "    ax2.set_xlabel('X (m)', fontsize=12)\n";
    script << "    ax2.set_ylabel('Water Saturation', fontsize=12)\n";
    script << "    ax2.set_title('Saturation Profile (y = mid)', fontsize=14, fontweight='bold')\n";
    script << "    ax2.set_ylim(0, 1)\n";
    script << "    ax2.grid(True, alpha=0.3)\n";
    script << "    ax2.legend(fontsize=10)\n";
    script << "    \n";
    script << "    plt.tight_layout()\n";
    script << "    plt.savefig(os.path.join(output_dir, f'saturation_{step:04d}.png'), dpi=150, bbox_inches='tight')\n";
    script << "    plt.close()\n";
    script << "    \n";
    script << "    print(f'Generated plot for step {step}')\n\n";
    
    script << "# Create animation frames summary plot\n";
    script << "if len(timesteps) >= 4:\n";
    script << "    fig, axes = plt.subplots(2, 2, figsize=(12, 10))\n";
    script << "    axes = axes.flatten()\n";
    script << "    \n";
    script << "    plot_steps = [timesteps[0], timesteps[len(timesteps)//3], \n";
    script << "                  timesteps[2*len(timesteps)//3], timesteps[-1]]\n";
    script << "    \n";
    script << "    for idx, step in enumerate(plot_steps):\n";
    script << "        filename = f'saturation_step_{step:04d}.csv'\n";
    script << "        df = pd.read_csv(os.path.join(output_dir, filename))\n";
    script << "        \n";
    script << "        x_unique = sorted(df['x'].unique())\n";
    script << "        y_unique = sorted(df['y'].unique())\n";
    script << "        nx, ny = len(x_unique), len(y_unique)\n";
    script << "        saturation = df['saturation'].values.reshape(ny, nx)\n";
    script << "        \n";
    script << "        im = axes[idx].imshow(saturation, extent=[min(x_unique), max(x_unique),\n";
    script << "                                                   min(y_unique), max(y_unique)],\n";
    script << "                              origin='lower', cmap='RdYlBu_r', vmin=0, vmax=1, aspect='auto')\n";
    script << "        axes[idx].set_title(f'Step {step}', fontsize=12, fontweight='bold')\n";
    script << "        axes[idx].set_xlabel('X (m)')\n";
    script << "        axes[idx].set_ylabel('Y (m)')\n";
    script << "        plt.colorbar(im, ax=axes[idx], label='Sw')\n";
    script << "    \n";
    script << "    plt.suptitle('Buckley-Leverett 2D Waterflooding Evolution', fontsize=16, fontweight='bold')\n";
    script << "    plt.tight_layout()\n";
    script << "    plt.savefig(os.path.join(output_dir, 'evolution_summary.png'), dpi=150, bbox_inches='tight')\n";
    script << "    plt.close()\n";
    script << "    print('Generated evolution summary plot')\n\n";
    
    script << "print('All plots generated successfully!')\n";
    script << "print(f'Plots saved to: {output_dir}')\n";
    
    script.close();
    
    // Make script executable
    chmod((output_dir + "/plot_results.py").c_str(), 0755);
}

int main(int argc, char** argv) {
    PetscErrorCode ierr;
    
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "========================================\n";
        std::cout << "  Buckley-Leverett 2D Waterflooding\n";
        std::cout << "========================================\n";
        std::cout << "\n";
    }
    
    // Create output directory
    std::string output_dir = "output_buckley_leverett";
    if (rank == 0) {
        mkdir(output_dir.c_str(), 0755);
        std::cout << "Output directory: " << output_dir << "\n\n";
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    
    // Domain parameters
    const double Lx = 100.0;  // Length in x (m)
    const double Ly = 20.0;   // Width in y (m)
    const int nx = 100;       // Grid points in x
    const int ny = 20;        // Grid points in y
    const double dx = Lx / nx;
    const double dy = Ly / ny;
    
    // Physical parameters
    const double phi = 0.25;        // Porosity
    const double k = 100e-15;       // Permeability (100 mD)
    const double mu_w = 0.5e-3;     // Water viscosity (Pa·s)
    const double mu_o = 2.0e-3;     // Oil viscosity (Pa·s)
    const double injection_rate = 1e-5; // m/s (Darcy velocity)
    
    // Time stepping
    const double total_time = 3600.0 * 24 * 30;  // 30 days
    const int num_steps = 50;
    const double dt = total_time / num_steps;
    
    if (rank == 0) {
        std::cout << "Domain: " << Lx << " m × " << Ly << " m\n";
        std::cout << "Grid: " << nx << " × " << ny << " cells\n";
        std::cout << "Porosity: " << phi << "\n";
        std::cout << "Permeability: " << k * 1e15 << " mD\n";
        std::cout << "Viscosity ratio (μ_o/μ_w): " << mu_o / mu_w << "\n";
        std::cout << "Time steps: " << num_steps << "\n";
        std::cout << "Total time: " << total_time / (3600 * 24) << " days\n";
        std::cout << "\n";
    }
    
    // Initialize saturation field
    std::vector<std::vector<double>> Sw(ny, std::vector<double>(nx, 0.2));  // Initial water saturation (connate)
    std::vector<std::vector<double>> Sw_new(ny, std::vector<double>(nx, 0.2));
    
    // Injection boundary: left side at x=0
    for (int j = 0; j < ny; ++j) {
        Sw[j][0] = 0.8;  // Injected water saturation
    }
    
    // Grid coordinates for output
    std::vector<double> x_coords(nx), y_coords(ny);
    for (int i = 0; i < nx; ++i) x_coords[i] = (i + 0.5) * dx;
    for (int j = 0; j < ny; ++j) y_coords[j] = (j + 0.5) * dy;
    
    if (rank == 0) {
        std::cout << "Starting simulation...\n\n";
        
        // Write initial condition
        writeCSV(output_dir + "/saturation_step_0000.csv", x_coords, y_coords, Sw, "saturation");
        std::cout << "Step 0 / " << num_steps << " - Time: 0.00 days\n";
    }
    
    // Time stepping loop
    for (int step = 1; step <= num_steps; ++step) {
        double current_time = step * dt;
        
        // Upwind finite volume scheme for saturation equation
        // ∂(φ*Sw)/∂t + ∂(u*fw)/∂x = 0
        
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                
                if (i == 0) {
                    // Injection boundary
                    Sw_new[j][i] = 0.8;
                    continue;
                }
                
                if (i == nx - 1) {
                    // Production boundary (outflow)
                    Sw_new[j][i] = Sw[j][i-1];
                    continue;
                }
                
                // Compute fractional flow at cell faces
                double fw_left = fractionalFlow(Sw[j][i-1], mu_w, mu_o);
                double fw_right = fractionalFlow(Sw[j][i], mu_w, mu_o);
                
                // Upwind flux (assuming flow from left to right)
                double flux_left = injection_rate * fw_left;
                double flux_right = injection_rate * fw_right;
                
                // Update saturation using explicit Euler
                double dSw = -(dt / (phi * dx)) * (flux_right - flux_left);
                
                Sw_new[j][i] = Sw[j][i] + dSw;
                
                // Clamp saturation to physical bounds
                Sw_new[j][i] = std::max(0.2, std::min(0.8, Sw_new[j][i]));
            }
        }
        
        // Update solution
        Sw = Sw_new;
        
        // Output every few steps
        if (rank == 0 && (step % std::max(1, num_steps / 20) == 0 || step == num_steps)) {
            char filename[256];
            snprintf(filename, sizeof(filename), "%s/saturation_step_%04d.csv", 
                    output_dir.c_str(), step);
            writeCSV(filename, x_coords, y_coords, Sw, "saturation");
            
            // Calculate average saturation
            double avg_sw = 0.0;
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    avg_sw += Sw[j][i];
                }
            }
            avg_sw /= (nx * ny);
            
            std::cout << "Step " << step << " / " << num_steps 
                     << " - Time: " << std::fixed << std::setprecision(2) 
                     << current_time / (3600 * 24) << " days"
                     << " - Avg Sw: " << std::setprecision(4) << avg_sw << "\n";
        }
    }
    
    if (rank == 0) {
        std::cout << "\nSimulation completed!\n";
        std::cout << "Generating plots...\n\n";
        
        // Generate plots for saved timesteps
        std::vector<int> plot_steps;
        for (int s = 0; s <= num_steps; ++s) {
            if (s % std::max(1, num_steps / 20) == 0 || s == num_steps) {
                plot_steps.push_back(s);
            }
        }
        
        for (size_t idx = 0; idx < plot_steps.size(); ++idx) {
            int s = plot_steps[idx];
            double t = s * dt;
            
            // Re-read the data for this step
            char csv_file[512];
            snprintf(csv_file, sizeof(csv_file), "%s/saturation_step_%04d.csv", 
                    output_dir.c_str(), s);
            
            std::ifstream csv(csv_file);
            if (csv.is_open()) {
                std::string header;
                std::getline(csv, header);
                
                std::vector<std::vector<double>> plot_data(ny, std::vector<double>(nx));
                
                for (int j = 0; j < ny; ++j) {
                    for (int i = 0; i < nx; ++i) {
                        double x_val, y_val, sw_val;
                        char comma;
                        csv >> x_val >> comma >> y_val >> comma >> sw_val;
                        plot_data[j][i] = sw_val;
                    }
                }
                csv.close();
                
                generateSaturationPlot(output_dir, s, plot_data, nx, ny, Lx, Ly, t / (3600 * 24));
                
                std::cout << "  Generated plot " << (idx + 1) << "/" << plot_steps.size() 
                         << " (step " << s << ")\n";
            }
        }
        
        std::cout << "\nPlots generated successfully!\n";
        std::cout << "Output files in: " << output_dir << "/\n";
        std::cout << "  - saturation_XXXX.ppm (raw image format)\n";
        std::cout << "  - saturation_XXXX.png (if ImageMagick 'convert' is available)\n";
        std::cout << "  - info_XXXX.txt (simulation metadata)\n";
        std::cout << "\n";
        
        // Also generate Python script as backup
        generatePlotScript(output_dir, num_steps);
        std::cout << "Python plotting script also generated: " << output_dir << "/plot_results.py\n";
        std::cout << "(requires numpy, matplotlib, pandas if you want to use it)\n";
        std::cout << "\n";
    }
    
    ierr = PetscFinalize();
    return ierr;
}
