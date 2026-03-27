#ifndef CAPROCK_INTEGRITY_HPP
#define CAPROCK_INTEGRITY_HPP

#include <array>
#include <cmath>

namespace FSRM {

class CaprockIntegrity {
public:
    struct IntegrityResult {
        double mohr_coulomb_safety_factor;
        double max_delta_CFS;
        bool failure_predicted;
        double critical_pressure;
        double estimated_magnitude;
    };

    static double mohrCoulombFOS(double sigma_n, double tau, double cohesion,
                                 double friction_angle);
    static double caprockPermeability(double k0, double sigma_eff,
                                      double sigma_eff_ref, double gamma,
                                      double enhancement_factor, bool failed);
    static double deltaCFS(double delta_tau, double delta_sigma_n, double friction,
                           double delta_pore_pressure);
    static double estimateMagnitude(double fault_area);

    IntegrityResult evaluate(const std::array<double, 6>& stress, double pore_pressure,
                             double cohesion, double friction_angle,
                             const std::array<double, 3>& fault_normal) const;

    void setCaprockProperties(double k0, double gamma, double enhancement);
    void setMohrCoulombParams(double cohesion, double friction_angle);

private:
    double k0_{1e-20};
    double gamma_{1e-7};
    double enhancement_{1000.0};
    double cohesion_{5e6};
    double friction_angle_{30.0};
};

} // namespace FSRM

#endif // CAPROCK_INTEGRITY_HPP
