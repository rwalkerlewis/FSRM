#ifndef RELATIVE_PERMEABILITY_HPP
#define RELATIVE_PERMEABILITY_HPP

/**
 * @file RelativePermeability.hpp
 * @brief ECLIPSE-compatible relative permeability models
 * 
 * Comprehensive relative permeability functionality:
 * - Two-phase models (Corey, Brooks-Corey, Van Genuchten)
 * - Three-phase models (Stone I, Stone II, Baker)
 * - Hysteresis (Killough, Carlson)
 * - End-point scaling (ENDSCALE)
 * - Tabular input (SWOF, SGOF, SOF3, SWFN, SGFN)
 * - Capillary pressure models
 * 
 * ECLIPSE Keywords: SWOF, SGOF, SLGOF, SOF2, SOF3, SWFN, SGFN, 
 *                   EHYSTR, ENDSCALE, SCALECRS
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <array>
#include <cmath>
#include <algorithm>
#include <functional>

namespace FSRM {

/**
 * @brief Relative permeability model type
 */
enum class RelPermModelType {
    COREY,              ///< Power-law (Corey) model
    BROOKS_COREY,       ///< Brooks-Corey model
    VAN_GENUCHTEN,      ///< Van Genuchten model
    LET,                ///< LET model (more flexible)
    TABULAR,            ///< Tabular from SWOF/SGOF etc.
    MODIFIED_COREY      ///< Modified Corey with curvature
};

/**
 * @brief Three-phase relative permeability model
 */
enum class ThreePhaseModel {
    STONE_I,            ///< Stone I method
    STONE_II,           ///< Stone II method
    BAKER,              ///< Baker linear interpolation
    ODD,                ///< Eclipse ODD model
    SEGREGATED          ///< Segregated flow model
};

/**
 * @brief Hysteresis model type
 */
enum class HysteresisModel {
    NONE,               ///< No hysteresis
    KILLOUGH,           ///< Killough model
    CARLSON,            ///< Carlson parallel curve model
    JARGON,             ///< Jargon model
    SIMPLE              ///< Simple scanning curve
};

/**
 * @brief Capillary pressure model type
 */
enum class CapillaryModel {
    NONE,               ///< No capillary pressure
    LINEAR,             ///< Linear Pc(S)
    BROOKS_COREY,       ///< Brooks-Corey Pc
    VAN_GENUCHTEN,      ///< Van Genuchten Pc
    LEVERETT_J,         ///< Leverett J-function
    TABULAR             ///< Tabular Pc
};

/**
 * @brief End-point data for a phase
 */
struct EndPoints {
    double S_connate = 0.0;     ///< Connate/irreducible saturation
    double S_critical = 0.0;    ///< Critical saturation
    double S_max = 1.0;         ///< Maximum saturation
    double kr_max = 1.0;        ///< Maximum relative permeability
    double kr_connate = 0.0;    ///< Kr at connate saturation
    double Pc_max = 0.0;        ///< Maximum capillary pressure
    double Pc_min = 0.0;        ///< Minimum capillary pressure
};

/**
 * @brief End-point scaling parameters
 */
struct EndPointScaling {
    // Water end-points
    double SWL = 0.0;           ///< Connate water saturation
    double SWCR = 0.0;          ///< Critical water saturation
    double SWU = 1.0;           ///< Maximum water saturation
    double SOWCR = 0.0;         ///< Critical oil saturation (water)
    double KRW_max = 1.0;       ///< Maximum water relative perm
    double KRO_OW_max = 1.0;    ///< Max oil rel perm at connate water
    
    // Gas end-points  
    double SGL = 0.0;           ///< Connate gas saturation
    double SGCR = 0.0;          ///< Critical gas saturation
    double SGU = 1.0;           ///< Maximum gas saturation
    double SOGCR = 0.0;         ///< Critical oil saturation (gas)
    double KRG_max = 1.0;       ///< Maximum gas relative perm
    double KRO_OG_max = 1.0;    ///< Max oil rel perm at connate gas
    
    // Capillary pressure scaling
    double PCW = 0.0;           ///< Water-oil Pc scaling
    double PCG = 0.0;           ///< Gas-oil Pc scaling
    
    bool defined = false;
};

/**
 * @brief Hysteresis state for tracking saturation history
 */
struct HysteresisState {
    double S_max_historical = 0.0;  ///< Maximum historical saturation
    double S_min_historical = 1.0;  ///< Minimum historical saturation
    bool increasing = true;          ///< Current direction of change
    double trapped_saturation = 0.0; ///< Trapped non-wetting saturation
    
    void update(double S_new) {
        if (S_new > S_max_historical) {
            S_max_historical = S_new;
            increasing = true;
        } else if (S_new < S_min_historical) {
            S_min_historical = S_new;
            increasing = false;
        }
    }
};

/**
 * @brief Two-phase relative permeability function base class
 */
class TwoPhaseRelPerm {
public:
    TwoPhaseRelPerm();
    virtual ~TwoPhaseRelPerm() = default;
    
    /**
     * @brief Calculate wetting phase relative permeability
     * @param Sw Wetting phase saturation
     * @return kr_w
     */
    virtual double kr_wetting(double Sw) const = 0;
    
    /**
     * @brief Calculate non-wetting phase relative permeability
     * @param Sw Wetting phase saturation (kr_nw is function of 1-Sw)
     * @return kr_nw
     */
    virtual double kr_nonwetting(double Sw) const = 0;
    
    /**
     * @brief Calculate derivatives for Jacobian
     * @param Sw Wetting phase saturation
     * @return {dkr_w/dSw, dkr_nw/dSw}
     */
    virtual std::pair<double, double> derivatives(double Sw) const;
    
    // End-point accessors
    void setEndPoints(const EndPoints& wet_ep, const EndPoints& nonwet_ep);
    const EndPoints& getWettingEndPoints() const { return ep_wetting_; }
    const EndPoints& getNonWettingEndPoints() const { return ep_nonwetting_; }
    
    // Normalized saturation
    double normalizedSaturation(double Sw) const;
    
protected:
    EndPoints ep_wetting_;
    EndPoints ep_nonwetting_;
};

/**
 * @brief Corey (power-law) relative permeability model
 */
class CoreyRelPerm : public TwoPhaseRelPerm {
public:
    CoreyRelPerm();
    
    double kr_wetting(double Sw) const override;
    double kr_nonwetting(double Sw) const override;
    
    // Set Corey exponents
    void setWettingExponent(double n) { n_wetting_ = n; }
    void setNonWettingExponent(double n) { n_nonwetting_ = n; }
    
    double getWettingExponent() const { return n_wetting_; }
    double getNonWettingExponent() const { return n_nonwetting_; }
    
private:
    double n_wetting_ = 2.0;
    double n_nonwetting_ = 2.0;
};

/**
 * @brief Brooks-Corey relative permeability model
 */
class BrooksCoreyRelPerm : public TwoPhaseRelPerm {
public:
    BrooksCoreyRelPerm();
    
    double kr_wetting(double Sw) const override;
    double kr_nonwetting(double Sw) const override;
    
    void setLambda(double lambda) { lambda_ = lambda; }
    double getLambda() const { return lambda_; }
    
private:
    double lambda_ = 2.0;  // Pore size distribution index
};

/**
 * @brief Van Genuchten relative permeability model
 */
class VanGenuchtenRelPerm : public TwoPhaseRelPerm {
public:
    VanGenuchtenRelPerm();
    
    double kr_wetting(double Sw) const override;
    double kr_nonwetting(double Sw) const override;
    
    void setParameters(double m, double n);
    void setParameterM(double m) { m_ = m; n_ = 1.0 / (1.0 - m); }
    
private:
    double m_ = 0.5;
    double n_ = 2.0;
};

/**
 * @brief LET model for more flexible curve shapes
 * 
 * kr = kr_max * S^L / (S^L + E*(1-S)^T)
 */
class LETRelPerm : public TwoPhaseRelPerm {
public:
    LETRelPerm();
    
    double kr_wetting(double Sw) const override;
    double kr_nonwetting(double Sw) const override;
    
    void setWettingLET(double L, double E, double T);
    void setNonWettingLET(double L, double E, double T);
    
private:
    double L_wet_ = 2.0, E_wet_ = 1.0, T_wet_ = 2.0;
    double L_nw_ = 2.0, E_nw_ = 1.0, T_nw_ = 2.0;
};

/**
 * @brief Tabular relative permeability from SWOF/SGOF
 */
class TabularRelPerm : public TwoPhaseRelPerm {
public:
    TabularRelPerm();
    
    double kr_wetting(double Sw) const override;
    double kr_nonwetting(double Sw) const override;
    
    // Set table data
    void setTable(const std::vector<double>& S,
                  const std::vector<double>& kr_w,
                  const std::vector<double>& kr_nw,
                  const std::vector<double>& Pc = {});
    
    // Get capillary pressure (if defined)
    double getPc(double Sw) const;
    bool hasPcData() const { return !Pc_.empty(); }
    
    // Parse from ECLIPSE format
    void parseFromSWOF(const std::vector<std::string>& data);
    void parseFromSGOF(const std::vector<std::string>& data);
    
private:
    std::vector<double> S_;
    std::vector<double> kr_w_;
    std::vector<double> kr_nw_;
    std::vector<double> Pc_;
    
    double interpolate(double S, const std::vector<double>& values) const;
};

/**
 * @brief Hysteresis wrapper for two-phase relative permeability
 */
class HysteresisRelPerm : public TwoPhaseRelPerm {
public:
    HysteresisRelPerm(std::shared_ptr<TwoPhaseRelPerm> drainage,
                      std::shared_ptr<TwoPhaseRelPerm> imbibition = nullptr,
                      HysteresisModel model = HysteresisModel::KILLOUGH);
    
    double kr_wetting(double Sw) const override;
    double kr_nonwetting(double Sw) const override;
    
    // Update hysteresis state (call each timestep)
    void updateState(double Sw);
    
    // Killough parameters
    void setKilloughCurvature(double C) { killough_curvature_ = C; }
    void setLandParameter(double C) { land_parameter_ = C; }
    
    // Get trapped saturation
    double getTrappedSaturation() const;
    
private:
    std::shared_ptr<TwoPhaseRelPerm> drainage_;
    std::shared_ptr<TwoPhaseRelPerm> imbibition_;
    HysteresisModel model_;
    mutable HysteresisState state_;
    
    double killough_curvature_ = 0.1;
    double land_parameter_ = 2.0;  // Land trapping parameter
    
    // Killough scanning curve
    double killoughScanning(double S, bool wetting) const;
    
    // Calculate trapped gas using Land equation
    double landTrapping(double S_max) const;
};

/**
 * @brief Three-phase relative permeability model
 */
class ThreePhaseRelPerm {
public:
    ThreePhaseRelPerm();
    
    // Set two-phase tables
    void setWaterOilTable(std::shared_ptr<TwoPhaseRelPerm> table);
    void setGasOilTable(std::shared_ptr<TwoPhaseRelPerm> table);
    void setGasWaterTable(std::shared_ptr<TwoPhaseRelPerm> table);
    
    void setThreePhaseModel(ThreePhaseModel model) { model_ = model; }
    
    /**
     * @brief Calculate three-phase relative permeabilities
     * @param Sw Water saturation
     * @param So Oil saturation
     * @param Sg Gas saturation (= 1 - Sw - So)
     * @return {kr_w, kr_o, kr_g}
     */
    std::array<double, 3> calculate(double Sw, double So, double Sg) const;
    
    /**
     * @brief Calculate with derivatives
     * @return {{kr_w, kr_o, kr_g}, {dkr_w/dSw, dkr_o/dSw, ...}}
     */
    std::pair<std::array<double, 3>, std::array<double, 6>> 
    calculateWithDerivatives(double Sw, double So, double Sg) const;
    
    // Individual phase accessors
    double kr_water(double Sw, double So, double Sg) const;
    double kr_oil(double Sw, double So, double Sg) const;
    double kr_gas(double Sw, double So, double Sg) const;
    
private:
    std::shared_ptr<TwoPhaseRelPerm> water_oil_;
    std::shared_ptr<TwoPhaseRelPerm> gas_oil_;
    std::shared_ptr<TwoPhaseRelPerm> gas_water_;
    
    ThreePhaseModel model_ = ThreePhaseModel::STONE_I;
    
    // Stone I method
    double stoneI_oil(double Sw, double So, double Sg) const;
    
    // Stone II method  
    double stoneII_oil(double Sw, double So, double Sg) const;
    
    // Baker linear interpolation
    double baker_oil(double Sw, double So, double Sg) const;
};

/**
 * @brief Capillary pressure model
 */
class CapillaryPressure {
public:
    CapillaryPressure();
    
    void setModel(CapillaryModel model) { model_ = model; }
    void setEntryPressure(double Pe) { entry_pressure_ = Pe; }
    void setParameters(double lambda_or_m);
    
    // Calculate Pc(Sw) for water-oil
    double Pc_ow(double Sw) const;
    
    // Calculate Pc(Sg) for gas-oil
    double Pc_go(double Sg) const;
    
    // Calculate Pc(Sw) for gas-water
    double Pc_gw(double Sw) const;
    
    // Set tabular data
    void setTabularPc(const std::vector<double>& S,
                      const std::vector<double>& Pc);
    
    // Leverett J-function
    void setJFunction(double k, double phi, double sigma, double theta);
    double JFunction(double Sw) const;
    double PcFromJ(double Sw, double k, double phi) const;
    
private:
    CapillaryModel model_ = CapillaryModel::BROOKS_COREY;
    double entry_pressure_ = 1e4;  // Pa
    double lambda_ = 2.0;          // Pore size distribution
    
    // Van Genuchten
    double m_vg_ = 0.5;
    double alpha_vg_ = 1e-4;       // 1/Pa
    
    // Tabular
    std::vector<double> S_table_;
    std::vector<double> Pc_table_;
    
    // Leverett J-function parameters
    double sigma_ = 0.03;          // Interfacial tension (N/m)
    double theta_ = 0.0;           // Contact angle (radians)
};

/**
 * @brief End-point scaler for relative permeability
 * 
 * Applies cell-by-cell end-point scaling as in ECLIPSE ENDSCALE
 */
class EndPointScaler {
public:
    EndPointScaler();
    
    void setReferenceEndPoints(const EndPointScaling& ref);
    void setCellEndPoints(int cell_id, const EndPointScaling& cell_ep);
    
    /**
     * @brief Scale relative permeability for a cell
     * 
     * @param cell_id Cell index
     * @param Sw Water saturation
     * @param kr_ref Reference kr values from table
     * @return Scaled kr values
     */
    std::array<double, 3> scaleKr(int cell_id, double Sw, double So, double Sg,
                                  const std::array<double, 3>& kr_ref) const;
    
    /**
     * @brief Scale saturation for lookup
     * 
     * Maps cell saturation to table saturation for interpolation.
     */
    double scaleWaterSaturation(int cell_id, double Sw) const;
    double scaleOilSaturation(int cell_id, double So) const;
    double scaleGasSaturation(int cell_id, double Sg) const;
    
    // Parse from ECLIPSE format
    void parseSCALECRS(const std::vector<std::string>& data);
    void parseENDSCALE(const std::vector<std::string>& data);
    void parseSWATINIT(const std::vector<double>& swatinit, int ncells);
    
private:
    EndPointScaling reference_;
    std::map<int, EndPointScaling> cell_endpoints_;
    
    // Vertical scaling mode
    bool use_three_point_ = true;
    bool use_two_point_ = false;
};

/**
 * @brief Relative permeability region manager
 * 
 * Handles multiple SATNUM regions with different rel perm tables
 */
class RelPermManager {
public:
    RelPermManager();
    
    // Add regions
    void addRegion(int satnum, std::shared_ptr<ThreePhaseRelPerm> relperm);
    void addCapillaryPressure(int satnum, std::shared_ptr<CapillaryPressure> pc);
    
    // Set default region
    void setDefaultRegion(int satnum) { default_region_ = satnum; }
    
    // Set SATNUM array
    void setSATNUM(const std::vector<int>& satnum) { satnum_ = satnum; }
    
    // Enable hysteresis
    void enableHysteresis(HysteresisModel model);
    void setHysteresisParameters(double curvature, double land_param);
    
    // Enable end-point scaling
    void enableEndPointScaling();
    void setEndPointScaling(int cell, const EndPointScaling& ep);
    
    /**
     * @brief Get relative permeabilities for a cell
     * 
     * @param cell Cell index
     * @param Sw Water saturation
     * @param So Oil saturation
     * @param Sg Gas saturation
     * @return {kr_w, kr_o, kr_g}
     */
    std::array<double, 3> getRelPerm(int cell, double Sw, double So, double Sg) const;
    
    /**
     * @brief Get capillary pressures for a cell
     * 
     * @param cell Cell index
     * @param Sw Water saturation
     * @param Sg Gas saturation
     * @return {Pc_ow, Pc_og}
     */
    std::pair<double, double> getCapillaryPressure(int cell, double Sw, double Sg) const;
    
    /**
     * @brief Update hysteresis state for all cells
     * 
     * @param saturations Water saturations for all cells
     */
    void updateHysteresis(const std::vector<double>& saturations);
    
    // Parse ECLIPSE keywords
    void parseSWOF(const std::vector<std::string>& data, int region);
    void parseSGOF(const std::vector<std::string>& data, int region);
    void parseSOF3(const std::vector<std::string>& data, int region);
    void parseSWFN(const std::vector<std::string>& data, int region);
    void parseSGFN(const std::vector<std::string>& data, int region);
    void parseSLGOF(const std::vector<std::string>& data, int region);
    void parseEHYSTR(const std::vector<std::string>& data);
    
    // Output
    void printRegionSummary(int region) const;
    void writeRelPermTable(std::ostream& os, int region) const;
    
private:
    std::map<int, std::shared_ptr<ThreePhaseRelPerm>> regions_;
    std::map<int, std::shared_ptr<CapillaryPressure>> cap_pressure_;
    std::vector<int> satnum_;
    int default_region_ = 1;
    
    // Hysteresis
    bool use_hysteresis_ = false;
    HysteresisModel hysteresis_model_ = HysteresisModel::KILLOUGH;
    std::map<int, HysteresisState> hysteresis_states_;
    double killough_curvature_ = 0.1;
    double land_parameter_ = 2.0;
    
    // End-point scaling
    bool use_scaling_ = false;
    std::unique_ptr<EndPointScaler> scaler_;
    
    int getRegion(int cell) const;
};

/**
 * @brief Factory to create relative permeability model from configuration
 */
std::shared_ptr<TwoPhaseRelPerm> createTwoPhaseRelPerm(
    RelPermModelType type,
    const std::map<std::string, double>& params);

std::shared_ptr<ThreePhaseRelPerm> createThreePhaseRelPerm(
    ThreePhaseModel model,
    std::shared_ptr<TwoPhaseRelPerm> water_oil,
    std::shared_ptr<TwoPhaseRelPerm> gas_oil);

// Parse helper functions
RelPermModelType parseRelPermModelType(const std::string& str);
ThreePhaseModel parseThreePhaseModel(const std::string& str);
HysteresisModel parseHysteresisModel(const std::string& str);
CapillaryModel parseCapillaryModel(const std::string& str);

} // namespace FSRM

#endif // RELATIVE_PERMEABILITY_HPP
