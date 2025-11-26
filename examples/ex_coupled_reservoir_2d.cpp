#include "Simulator.hpp"
#include "PhysicsKernel.hpp"
#include "Visualization.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sys/stat.h>

using namespace FSRM;

// ============================================================================
// PPM Image Generation with Advanced Features
// ============================================================================

struct Color {
    unsigned char r, g, b;
    
    Color(unsigned char red = 0, unsigned char green = 0, unsigned char blue = 0)
        : r(red), g(green), b(blue) {}
    
    static Color lerp(const Color& c1, const Color& c2, double t) {
        t = std::max(0.0, std::min(1.0, t));
        return Color(
            static_cast<unsigned char>(c1.r + t * (c2.r - c1.r)),
            static_cast<unsigned char>(c1.g + t * (c2.g - c1.g)),
            static_cast<unsigned char>(c1.b + t * (c2.b - c1.b))
        );
    }
};

class ColorMap {
public:
    virtual Color getColor(double value) const = 0;
    virtual ~ColorMap() = default;
};

class JetColorMap : public ColorMap {
public:
    Color getColor(double value) const override {
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
};

class ViridisColorMap : public ColorMap {
public:
    Color getColor(double value) const override {
        value = std::max(0.0, std::min(1.0, value));
        
        // Viridis-inspired color palette
        const std::vector<Color> palette = {
            Color(68, 1, 84),      // Dark purple
            Color(59, 82, 139),    // Blue-purple
            Color(33, 145, 140),   // Teal
            Color(94, 201, 98),    // Green
            Color(253, 231, 37)    // Yellow
        };
        
        int n = palette.size() - 1;
        double scaled = value * n;
        int idx = static_cast<int>(scaled);
        idx = std::max(0, std::min(n - 1, idx));
        
        double t = scaled - idx;
        return Color::lerp(palette[idx], palette[idx + 1], t);
    }
};

class PressureColorMap : public ColorMap {
public:
    Color getColor(double value) const override {
        value = std::max(0.0, std::min(1.0, value));
        
        // Blue (low) -> cyan -> green -> yellow -> red (high)
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
};

class SubsidenceColorMap : public ColorMap {
public:
    Color getColor(double value) const override {
        value = std::max(0.0, std::min(1.0, value));
        
        // White (no subsidence) -> yellow -> orange -> red -> dark red (max subsidence)
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
};

// ============================================================================
// Enhanced Plot Generator
// ============================================================================

class PlotGenerator {
public:
    PlotGenerator(int width, int height, int margin_left, int margin_right,
                  int margin_top, int margin_bottom)
        : plot_width(width), plot_height(height),
          margin_l(margin_left), margin_r(margin_right),
          margin_t(margin_top), margin_b(margin_bottom) {
        
        field_width = plot_width - margin_l - margin_r;
        field_height = plot_height - margin_t - margin_b;
        
        // Initialize with white background
        image.resize(plot_height, std::vector<Color>(plot_width, Color(255, 255, 255)));
    }
    
    void drawField(const std::vector<std::vector<double>>& data,
                   const ColorMap& cmap, double vmin, double vmax) {
        int ny = data.size();
        int nx = data[0].size();
        
        for (int py = 0; py < field_height; ++py) {
            for (int px = 0; px < field_width; ++px) {
                int i = px * nx / field_width;
                int j = py * ny / field_height;
                
                double value = (data[j][i] - vmin) / (vmax - vmin);
                Color color = cmap.getColor(value);
                
                int img_y = margin_t + (field_height - 1 - py);  // Flip vertically
                int img_x = margin_l + px;
                
                image[img_y][img_x] = color;
            }
        }
    }
    
    void drawColorbar(const ColorMap& cmap, double vmin, double vmax,
                     int cbar_x, int cbar_width, const std::string& label = "") {
        int cbar_height = field_height;
        
        // Draw colorbar gradient
        for (int py = 0; py < cbar_height; ++py) {
            double value = 1.0 - (double)py / cbar_height;
            Color color = cmap.getColor(value);
            
            for (int px = 0; px < cbar_width; ++px) {
                int img_y = margin_t + py;
                int img_x = cbar_x + px;
                if (img_x < plot_width && img_y < plot_height) {
                    image[img_y][img_x] = color;
                }
            }
        }
        
        // Draw colorbar border
        drawRect(cbar_x - 1, margin_t - 1,
                cbar_x + cbar_width, margin_t + cbar_height, Color(0, 0, 0), 1);
    }
    
    void drawWell(double x, double y, double Lx, double Ly, bool is_injector,
                  const std::string& label = "") {
        // Convert to image coordinates
        int img_x = margin_l + static_cast<int>((x / Lx) * field_width);
        int img_y = margin_t + field_height - static_cast<int>((y / Ly) * field_height);
        
        Color well_color = is_injector ? Color(0, 200, 0) : Color(200, 0, 0);
        int radius = 8;
        
        // Draw filled circle for well
        for (int dy = -radius; dy <= radius; ++dy) {
            for (int dx = -radius; dx <= radius; ++dx) {
                if (dx*dx + dy*dy <= radius*radius) {
                    int px = img_x + dx;
                    int py = img_y + dy;
                    if (px >= 0 && px < plot_width && py >= 0 && py < plot_height) {
                        image[py][px] = well_color;
                    }
                }
            }
        }
        
        // Draw black outline
        for (int angle = 0; angle < 360; angle += 10) {
            double rad = angle * M_PI / 180.0;
            int px = img_x + static_cast<int>(radius * cos(rad));
            int py = img_y + static_cast<int>(radius * sin(rad));
            if (px >= 0 && px < plot_width && py >= 0 && py < plot_height) {
                image[py][px] = Color(0, 0, 0);
            }
        }
    }
    
    void drawRect(int x0, int y0, int x1, int y1, const Color& color, int thickness = 1) {
        for (int t = 0; t < thickness; ++t) {
            // Horizontal lines
            for (int x = x0 - t; x <= x1 + t; ++x) {
                if (x >= 0 && x < plot_width) {
                    if (y0 - t >= 0 && y0 - t < plot_height)
                        image[y0 - t][x] = color;
                    if (y1 + t >= 0 && y1 + t < plot_height)
                        image[y1 + t][x] = color;
                }
            }
            // Vertical lines
            for (int y = y0 - t; y <= y1 + t; ++y) {
                if (y >= 0 && y < plot_height) {
                    if (x0 - t >= 0 && x0 - t < plot_width)
                        image[y][x0 - t] = color;
                    if (x1 + t >= 0 && x1 + t < plot_width)
                        image[y][x1 + t] = color;
                }
            }
        }
    }
    
    void drawGrid(double Lx, double Ly, int num_x_lines, int num_y_lines) {
        Color grid_color(200, 200, 200);
        
        // Vertical grid lines
        for (int i = 0; i <= num_x_lines; ++i) {
            int x = margin_l + (i * field_width) / num_x_lines;
            for (int y = margin_t; y < margin_t + field_height; y += 2) {
                if (x >= 0 && x < plot_width && y >= 0 && y < plot_height) {
                    image[y][x] = grid_color;
                }
            }
        }
        
        // Horizontal grid lines
        for (int j = 0; j <= num_y_lines; ++j) {
            int y = margin_t + (j * field_height) / num_y_lines;
            for (int x = margin_l; x < margin_l + field_width; x += 2) {
                if (x >= 0 && x < plot_width && y >= 0 && y < plot_height) {
                    image[y][x] = grid_color;
                }
            }
        }
    }
    
    void drawTitle(const std::string& title) {
        // Title will be added in metadata file
        current_title = title;
    }
    
    void writePPM(const std::string& filename) const {
        std::ofstream file(filename, std::ios::binary);
        file << "P6\n" << plot_width << " " << plot_height << "\n255\n";
        
        for (int y = 0; y < plot_height; ++y) {
            for (int x = 0; x < plot_width; ++x) {
                file.write(reinterpret_cast<const char*>(&image[y][x].r), 1);
                file.write(reinterpret_cast<const char*>(&image[y][x].g), 1);
                file.write(reinterpret_cast<const char*>(&image[y][x].b), 1);
            }
        }
        file.close();
    }
    
    bool convertToPNG(const std::string& ppm_file, const std::string& png_file) const {
        std::string cmd = "command -v convert > /dev/null 2>&1 && convert " + 
                         ppm_file + " " + png_file + " 2>/dev/null && rm -f " + ppm_file;
        return (system(cmd.c_str()) == 0);
    }
    
private:
    int plot_width, plot_height;
    int margin_l, margin_r, margin_t, margin_b;
    int field_width, field_height;
    std::vector<std::vector<Color>> image;
    std::string current_title;
};

// ============================================================================
// Simulation Data Structures
// ============================================================================

struct Well {
    double x, y;           // Position
    bool is_injector;
    double rate;           // m³/s
    std::string name;
    
    Well(double x_pos, double y_pos, bool inj, double r, const std::string& n)
        : x(x_pos), y(y_pos), is_injector(inj), rate(r), name(n) {}
};

struct SimulationState {
    std::vector<std::vector<double>> pressure;
    std::vector<std::vector<double>> saturation;
    std::vector<std::vector<double>> displacement_x;
    std::vector<std::vector<double>> displacement_y;
    std::vector<std::vector<double>> subsidence;
    std::vector<std::vector<double>> porosity;
    std::vector<std::vector<double>> permeability;
    
    int nx, ny;
    double time;
    
    SimulationState(int nx_, int ny_) : nx(nx_), ny(ny_), time(0.0) {
        pressure.resize(ny, std::vector<double>(nx, 0.0));
        saturation.resize(ny, std::vector<double>(nx, 0.0));
        displacement_x.resize(ny, std::vector<double>(nx, 0.0));
        displacement_y.resize(ny, std::vector<double>(nx, 0.0));
        subsidence.resize(ny, std::vector<double>(nx, 0.0));
        porosity.resize(ny, std::vector<double>(nx, 0.0));
        permeability.resize(ny, std::vector<double>(nx, 0.0));
    }
};

// ============================================================================
// Physics Models
// ============================================================================

class CoupledReservoirSimulator {
public:
    CoupledReservoirSimulator(int nx, int ny, double Lx, double Ly)
        : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly),
          dx_(Lx / nx), dy_(Ly / ny) {
        
        state_ = std::make_unique<SimulationState>(nx, ny);
        
        // Initialize with default properties
        phi0_ = 0.2;
        k0_ = 100e-15;  // 100 mD
        E_ = 10e9;      // 10 GPa
        nu_ = 0.25;
        rho_s_ = 2500.0;
        rho_f_ = 1000.0;
        mu_ = 1e-3;
        alpha_ = 0.7;   // Biot coefficient
        ct_ = 1e-9;     // Total compressibility
        
        P_initial_ = 30e6;  // 30 MPa
        
        initializeState();
    }
    
    void addWell(const Well& well) {
        wells_.push_back(well);
    }
    
    void setRockProperties(double phi, double k, double E, double nu) {
        phi0_ = phi;
        k0_ = k;
        E_ = E;
        nu_ = nu;
    }
    
    void setFluidProperties(double rho, double mu) {
        rho_f_ = rho;
        mu_ = mu;
    }
    
    void initializeState() {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                state_->pressure[j][i] = P_initial_;
                state_->saturation[j][i] = 0.2;  // Connate water
                state_->porosity[j][i] = phi0_;
                state_->permeability[j][i] = k0_;
            }
        }
    }
    
    void step(double dt) {
        // Solve flow equation
        solveFlow(dt);
        
        // Solve geomechanics
        solveGeomechanics(dt);
        
        // Update porosity and permeability from geomechanics
        updatePoroPermeFromGeomechanics();
        
        state_->time += dt;
    }
    
    const SimulationState& getState() const { return *state_; }
    const std::vector<Well>& getWells() const { return wells_; }
    
private:
    void solveFlow(double dt) {
        auto P_new = state_->pressure;
        auto Sw_new = state_->saturation;
        
        // Simple finite difference for demonstration
        for (int j = 1; j < ny_ - 1; ++j) {
            for (int i = 1; i < nx_ - 1; ++i) {
                double P = state_->pressure[j][i];
                double Px_plus = state_->pressure[j][i+1];
                double Px_minus = state_->pressure[j][i-1];
                double Py_plus = state_->pressure[j+1][i];
                double Py_minus = state_->pressure[j-1][i];
                
                double phi = state_->porosity[j][i];
                double k = state_->permeability[j][i];
                
                // Darcy flux
                double qx = -(k / mu_) * (Px_plus - Px_minus) / (2 * dx_);
                double qy = -(k / mu_) * (Py_plus - Py_minus) / (2 * dy_);
                
                // Divergence
                double div_q = (qx / dx_) + (qy / dy_);
                
                // Pressure update (compressible flow)
                double dP = -(dt / (phi * ct_)) * div_q;
                P_new[j][i] = P + dP;
                
                // Saturation update (simplified two-phase)
                double fw = fractionalFlow(state_->saturation[j][i]);
                double flux_x = qx * fw;
                
                double Sw_upwind = (qx > 0) ? state_->saturation[j][i-1] : state_->saturation[j][i+1];
                double dSw = -(dt / phi) * flux_x / dx_;
                
                Sw_new[j][i] = state_->saturation[j][i] + dSw;
                Sw_new[j][i] = std::max(0.2, std::min(0.8, Sw_new[j][i]));
            }
        }
        
        // Apply well conditions
        applyWellConditions(P_new, Sw_new);
        
        state_->pressure = P_new;
        state_->saturation = Sw_new;
    }
    
    void solveGeomechanics(double dt) {
        // Compute subsidence from pressure depletion
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                double dP = state_->pressure[j][i] - P_initial_;
                
                // Uniaxial compaction model
                double K = E_ / (3 * (1 - 2 * nu_));  // Bulk modulus
                double volumetric_strain = alpha_ * dP / K;
                
                // Vertical displacement (subsidence)
                double cell_height = Ly_ / ny_;
                double dz = volumetric_strain * cell_height;
                
                state_->subsidence[j][i] = -dz;  // Negative is downward
                
                // Horizontal displacement (small for uniaxial)
                double lateral_strain = nu_ / (1 - nu_) * volumetric_strain;
                state_->displacement_x[j][i] = lateral_strain * (i * dx_ - Lx_ / 2);
                state_->displacement_y[j][i] = lateral_strain * (j * dy_ - Ly_ / 2);
            }
        }
    }
    
    void updatePoroPermeFromGeomechanics() {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                double dP = state_->pressure[j][i] - P_initial_;
                
                // Porosity change from pressure
                double K = E_ / (3 * (1 - 2 * nu_));
                double volumetric_strain = alpha_ * dP / K;
                state_->porosity[j][i] = phi0_ * (1 + volumetric_strain);
                
                // Permeability change (Kozeny-Carman type)
                double phi = state_->porosity[j][i];
                double phi_ratio = phi / phi0_;
                state_->permeability[j][i] = k0_ * std::pow(phi_ratio, 3) * 
                                            std::pow((1 - phi0_) / (1 - phi), 2);
            }
        }
    }
    
    void applyWellConditions(std::vector<std::vector<double>>& P,
                            std::vector<std::vector<double>>& Sw) {
        for (const auto& well : wells_) {
            int i = static_cast<int>(well.x / dx_);
            int j = static_cast<int>(well.y / dy_);
            
            i = std::max(0, std::min(nx_ - 1, i));
            j = std::max(0, std::min(ny_ - 1, j));
            
            if (well.is_injector) {
                P[j][i] = P_initial_ + 5e6;  // 5 MPa overpressure
                Sw[j][i] = 0.8;  // Inject water
            } else {
                P[j][i] = P_initial_ - 10e6;  // 10 MPa drawdown
                // Production well maintains current saturation
            }
        }
    }
    
    double fractionalFlow(double Sw) const {
        double Sw_norm = (Sw - 0.2) / 0.6;
        Sw_norm = std::max(0.0, std::min(1.0, Sw_norm));
        
        double krw = std::pow(Sw_norm, 2.0);
        double kro = std::pow(1.0 - Sw_norm, 2.0);
        
        double mu_w = mu_;
        double mu_o = 2.0 * mu_;
        
        double lambda_w = krw / mu_w;
        double lambda_o = kro / mu_o;
        
        if (lambda_w + lambda_o < 1e-12) return 0.0;
        return lambda_w / (lambda_w + lambda_o);
    }
    
    int nx_, ny_;
    double Lx_, Ly_;
    double dx_, dy_;
    
    double phi0_, k0_;
    double E_, nu_;
    double rho_s_, rho_f_;
    double mu_;
    double alpha_, ct_;
    double P_initial_;
    
    std::unique_ptr<SimulationState> state_;
    std::vector<Well> wells_;
};

// ============================================================================
// Visualization
// ============================================================================

void generateCombinedPlot(const std::string& output_dir, int step,
                         const SimulationState& state,
                         const std::vector<Well>& wells,
                         double Lx, double Ly) {
    const int plot_width = 1600;
    const int plot_height = 1200;
    const int margin = 100;
    const int margin_right = 200;
    
    // Create 2x2 subplot grid
    const int subplot_width = (plot_width - 3 * margin - margin_right) / 2;
    const int subplot_height = (plot_height - 3 * margin) / 2;
    
    std::vector<std::vector<Color>> image(plot_height,
        std::vector<Color>(plot_width, Color(240, 240, 240)));
    
    // Define subplots
    struct Subplot {
        int x0, y0, width, height;
        std::string title;
    };
    
    std::vector<Subplot> subplots = {
        {margin, margin, subplot_width, subplot_height, "Pressure (MPa)"},
        {2 * margin + subplot_width, margin, subplot_width, subplot_height, "Water Saturation"},
        {margin, 2 * margin + subplot_height, subplot_width, subplot_height, "Subsidence (m)"},
        {2 * margin + subplot_width, 2 * margin + subplot_height, subplot_width, subplot_height, "Permeability (mD)"}
    };
    
    // Color maps
    PressureColorMap pressure_cmap;
    JetColorMap saturation_cmap;
    SubsidenceColorMap subsidence_cmap;
    ViridisColorMap perm_cmap;
    
    // Get value ranges
    double P_min = 30e6;  // Initial pressure
    double P_max = 30e6;
    double sub_min = 0.0;
    double sub_max = 0.0;
    double k_min = 100e-15;
    double k_max = 100e-15;
    
    for (int j = 0; j < state.ny; ++j) {
        for (int i = 0; i < state.nx; ++i) {
            P_min = std::min(P_min, state.pressure[j][i]);
            P_max = std::max(P_max, state.pressure[j][i]);
            sub_min = std::min(sub_min, state.subsidence[j][i]);
            sub_max = std::max(sub_max, state.subsidence[j][i]);
            k_min = std::min(k_min, state.permeability[j][i]);
            k_max = std::max(k_max, state.permeability[j][i]);
        }
    }
    
    // Draw each subplot
    auto drawSubplot = [&](int subplot_idx, const std::vector<std::vector<double>>& data,
                          const ColorMap& cmap, double vmin, double vmax) {
        const auto& sp = subplots[subplot_idx];
        
        for (int py = 0; py < sp.height; ++py) {
            for (int px = 0; px < sp.width; ++px) {
                int i = px * state.nx / sp.width;
                int j = (sp.height - 1 - py) * state.ny / sp.height;
                
                double value = (data[j][i] - vmin) / (vmax - vmin);
                Color color = cmap.getColor(value);
                
                image[sp.y0 + py][sp.x0 + px] = color;
            }
        }
        
        // Draw wells
        for (const auto& well : wells) {
            int well_x = sp.x0 + static_cast<int>((well.x / Lx) * sp.width);
            int well_y = sp.y0 + sp.height - static_cast<int>((well.y / Ly) * sp.height);
            
            Color well_color = well.is_injector ? Color(0, 255, 0) : Color(255, 0, 0);
            int radius = 6;
            
            for (int dy = -radius; dy <= radius; ++dy) {
                for (int dx = -radius; dx <= radius; ++dx) {
                    if (dx*dx + dy*dy <= radius*radius) {
                        int px = well_x + dx;
                        int py = well_y + dy;
                        if (px >= sp.x0 && px < sp.x0 + sp.width &&
                            py >= sp.y0 && py < sp.y0 + sp.height) {
                            image[py][px] = well_color;
                        }
                    }
                }
            }
        }
        
        // Draw border
        for (int x = sp.x0; x < sp.x0 + sp.width; ++x) {
            image[sp.y0][x] = Color(0, 0, 0);
            image[sp.y0 + sp.height - 1][x] = Color(0, 0, 0);
        }
        for (int y = sp.y0; y < sp.y0 + sp.height; ++y) {
            image[y][sp.x0] = Color(0, 0, 0);
            image[y][sp.x0 + sp.width - 1] = Color(0, 0, 0);
        }
    };
    
    // Draw all subplots
    drawSubplot(0, state.pressure, pressure_cmap, P_min, P_max);
    drawSubplot(1, state.saturation, saturation_cmap, 0.0, 1.0);
    drawSubplot(2, state.subsidence, subsidence_cmap, sub_min, sub_max);
    
    // Convert permeability to mD for plotting
    std::vector<std::vector<double>> k_mD = state.permeability;
    for (int j = 0; j < state.ny; ++j) {
        for (int i = 0; i < state.nx; ++i) {
            k_mD[j][i] *= 1e15;  // Convert to mD
        }
    }
    drawSubplot(3, k_mD, perm_cmap, k_min * 1e15, k_max * 1e15);
    
    // Write to file
    char filename_ppm[512];
    snprintf(filename_ppm, sizeof(filename_ppm), 
            "%s/coupled_sim_%04d.ppm", output_dir.c_str(), step);
    
    std::ofstream file(filename_ppm, std::ios::binary);
    file << "P6\n" << plot_width << " " << plot_height << "\n255\n";
    
    for (int y = 0; y < plot_height; ++y) {
        for (int x = 0; x < plot_width; ++x) {
            file.write(reinterpret_cast<const char*>(&image[y][x].r), 1);
            file.write(reinterpret_cast<const char*>(&image[y][x].g), 1);
            file.write(reinterpret_cast<const char*>(&image[y][x].b), 1);
        }
    }
    file.close();
    
    // Convert to PNG
    char filename_png[512];
    snprintf(filename_png, sizeof(filename_png),
            "%s/coupled_sim_%04d.png", output_dir.c_str(), step);
    
    std::string cmd = "command -v convert > /dev/null 2>&1 && convert " + 
                     std::string(filename_ppm) + " " + filename_png + 
                     " 2>/dev/null && rm -f " + filename_ppm;
    system(cmd.c_str());
    
    // Write metadata
    char metadata_file[512];
    snprintf(metadata_file, sizeof(metadata_file),
            "%s/info_%04d.txt", output_dir.c_str(), step);
    
    std::ofstream info(metadata_file);
    info << "═══════════════════════════════════════════════════\n";
    info << "  Coupled Reservoir Simulation - Step " << step << "\n";
    info << "═══════════════════════════════════════════════════\n\n";
    info << "Time: " << std::fixed << std::setprecision(2) 
         << state.time / (3600 * 24) << " days\n\n";
    info << "Pressure:\n";
    info << "  Min: " << std::setprecision(2) << P_min / 1e6 << " MPa\n";
    info << "  Max: " << std::setprecision(2) << P_max / 1e6 << " MPa\n\n";
    info << "Subsidence:\n";
    info << "  Max: " << std::setprecision(4) << std::abs(sub_min) << " m\n\n";
    info << "Permeability:\n";
    info << "  Min: " << std::setprecision(1) << k_min * 1e15 << " mD\n";
    info << "  Max: " << std::setprecision(1) << k_max * 1e15 << " mD\n\n";
    info << "Wells:\n";
    for (const auto& well : wells) {
        info << "  " << well.name << " (" 
             << (well.is_injector ? "Injector" : "Producer") << ")\n";
        info << "    Position: (" << well.x << ", " << well.y << ") m\n";
    }
    info.close();
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    if (rank == 0) {
        std::cout << "\n";
        std::cout << "═══════════════════════════════════════════════════\n";
        std::cout << "  2D Coupled Reservoir Simulation\n";
        std::cout << "  Flow + Wells + Geomechanics\n";
        std::cout << "═══════════════════════════════════════════════════\n\n";
    }
    
    // Create output directory
    std::string output_dir = "output_coupled_2d";
    if (rank == 0) {
        mkdir(output_dir.c_str(), 0755);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    
    // Domain setup
    const double Lx = 500.0;   // 500 m
    const double Ly = 200.0;   // 200 m
    const int nx = 150;
    const int ny = 60;
    
    // Time setup
    const double total_time = 3600.0 * 24 * 365;  // 1 year
    const int num_steps = 100;
    const double dt = total_time / num_steps;
    
    if (rank == 0) {
        std::cout << "Domain: " << Lx << " m × " << Ly << " m\n";
        std::cout << "Grid: " << nx << " × " << ny << " cells\n";
        std::cout << "Time steps: " << num_steps << "\n";
        std::cout << "Total time: " << total_time / (3600 * 24 * 365) << " years\n\n";
    }
    
    // Create simulator
    CoupledReservoirSimulator sim(nx, ny, Lx, Ly);
    
    // Set properties
    sim.setRockProperties(0.2, 100e-15, 10e9, 0.25);
    sim.setFluidProperties(1000.0, 1e-3);
    
    // Add wells
    sim.addWell(Well(100.0, 100.0, true, 0.01, "INJ-1"));   // Injector
    sim.addWell(Well(400.0, 50.0, false, -0.01, "PROD-1"));  // Producer north
    sim.addWell(Well(400.0, 150.0, false, -0.01, "PROD-2")); // Producer south
    
    if (rank == 0) {
        std::cout << "Wells:\n";
        for (const auto& well : sim.getWells()) {
            std::cout << "  " << well.name << " (" 
                     << (well.is_injector ? "Injector" : "Producer") << ")\n";
            std::cout << "    Position: (" << well.x << ", " << well.y << ") m\n";
        }
        std::cout << "\nStarting simulation...\n\n";
    }
    
    // Time loop
    int plot_interval = num_steps / 20;
    
    for (int step = 0; step <= num_steps; ++step) {
        if (step > 0) {
            sim.step(dt);
        }
        
        // Output
        if (rank == 0 && (step % plot_interval == 0 || step == num_steps)) {
            generateCombinedPlot(output_dir, step, sim.getState(), 
                               sim.getWells(), Lx, Ly);
            
            double time_days = sim.getState().time / (3600 * 24);
            std::cout << "Step " << std::setw(4) << step << " / " << num_steps
                     << "  |  Time: " << std::fixed << std::setprecision(1)
                     << std::setw(8) << time_days << " days  |  Plot generated\n";
        }
    }
    
    if (rank == 0) {
        std::cout << "\n═══════════════════════════════════════════════════\n";
        std::cout << "  Simulation Complete!\n";
        std::cout << "═══════════════════════════════════════════════════\n\n";
        std::cout << "Output directory: " << output_dir << "/\n";
        std::cout << "  • coupled_sim_XXXX.png - 2×2 subplot visualization\n";
        std::cout << "  • info_XXXX.txt - Detailed statistics\n\n";
        std::cout << "Plots show:\n";
        std::cout << "  [Top Left]     Pressure distribution (MPa)\n";
        std::cout << "  [Top Right]    Water saturation\n";
        std::cout << "  [Bottom Left]  Surface subsidence (m)\n";
        std::cout << "  [Bottom Right] Permeability changes (mD)\n\n";
        std::cout << "Wells indicated by colored circles:\n";
        std::cout << "  • Green = Injector\n";
        std::cout << "  • Red = Producer\n\n";
    }
    
    ierr = PetscFinalize();
    return ierr;
}
