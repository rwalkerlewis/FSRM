/**
 * @file ThermalPressurization.cpp
 * @brief Implementation of Thermal Pressurization model
 */

#include "ThermalPressurization.hpp"
#include <cmath>

namespace FSRM {

ThermalPressurization::ThermalPressurization() : enabled(false) {
}

void ThermalPressurization::setMaterialProperties(const MaterialProperties& properties) {
    props = properties;
}

void ThermalPressurization::initialize(State& state) const {
    state.temperature = props.initial_temperature;
    state.pore_pressure = props.initial_pore_pressure;
    state.cumulative_slip = 0.0;
    state.slip_rate = 0.0;
    state.shear_stress = 0.0;
    state.normal_stress = 0.0;
}

void ThermalPressurization::update(State& state, double dt) {
    if (!enabled || dt <= 0.0) {
        return;
    }
    
    // Compute frictional heating rate
    double Q = getFrictionalHeatingRate(state);
    
    // Temperature rise from shear heating (adiabatic approximation for short time)
    double dT = Q * dt / props.specific_heat;
    
    // Pore pressure rise from thermal pressurization
    // Λ = (∂p/∂T)_V, the thermal pressurization coefficient
    double dp = props.Lambda * dT;
    
    // Apply diffusion corrections (simplified)
    double t_diff_th = props.half_width * props.half_width / props.thermal_diffusivity;
    double t_diff_hy = props.half_width * props.half_width / props.hydraulic_diffusivity;
    
    // Thermal diffusion reduces temperature rise
    if (dt < t_diff_th) {
        // Short time: use full adiabatic approximation
    } else {
        // Longer time: apply diffusion correction
        dT *= std::exp(-dt / t_diff_th);
    }
    
    // Hydraulic diffusion reduces pore pressure
    if (dt < t_diff_hy) {
        // Short time: undrained
    } else {
        // Longer time: apply diffusion
        dp *= std::exp(-dt / t_diff_hy);
    }
    
    // Update state
    state.temperature += dT;
    state.pore_pressure += dp;
    state.cumulative_slip += state.slip_rate * dt;
}

double ThermalPressurization::getEffectiveNormalStressChange(const State& state) const {
    // Effective normal stress change = -Δp (pore pressure reduces effective stress)
    return -state.pore_pressure + props.initial_pore_pressure;
}

double ThermalPressurization::getFrictionalHeatingRate(const State& state) const {
    // Q = τ * V / w (heating rate per unit volume in shear zone)
    // where τ is shear stress, V is slip rate, w is shear zone half-width
    if (props.half_width <= 0.0) {
        return 0.0;
    }
    return state.shear_stress * state.slip_rate / (2.0 * props.half_width);
}

ThermalPressurization::AnalyticalSolution ThermalPressurization::getAnalyticalSolution(
    double slip, double slip_rate, double shear_stress, double time) const {
    
    (void)slip;  // Unused but part of interface
    
    AnalyticalSolution sol;
    sol.temperature = props.initial_temperature;
    sol.pore_pressure = props.initial_pore_pressure;
    
    if (time <= 0.0 || slip_rate <= 0.0) {
        return sol;
    }
    
    // Frictional heating rate
    double Q = shear_stress * slip_rate / (2.0 * props.half_width);
    
    // Characteristic diffusion lengths
    double L_th = std::sqrt(props.thermal_diffusivity * time);
    double L_hy = std::sqrt(props.hydraulic_diffusivity * time);
    
    // Temperature rise (using half-space solution for w << L_th)
    double dT;
    if (L_th > props.half_width) {
        // Slip-on-plane approximation
        dT = Q * time / props.specific_heat;
        dT *= props.half_width / L_th;  // Diffusion correction
    } else {
        // Adiabatic limit
        dT = Q * time / props.specific_heat;
    }
    
    // Pore pressure from thermal pressurization
    double dp = props.Lambda * dT;
    
    // Apply hydraulic diffusion correction
    if (L_hy > props.half_width) {
        dp *= props.half_width / L_hy;
    }
    
    sol.temperature = props.initial_temperature + dT;
    sol.pore_pressure = props.initial_pore_pressure + dp;
    
    return sol;
}

} // namespace FSRM
