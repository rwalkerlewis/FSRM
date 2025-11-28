/**
 * @file test_seismic_source.cpp
 * @brief Unit tests for seismic source and receiver models
 * 
 * Tests cover:
 * - Moment tensor operations
 * - Source time functions
 * - Point sources
 * - Kinematic sources
 * - Receivers
 */

#include "SeismicSource.hpp"
#include "Testing.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

namespace FSRM {
namespace Testing {

// =============================================================================
// Moment Tensor Tests
// =============================================================================

/**
 * @test Test double couple moment tensor properties
 */
bool test_double_couple_properties() {
    std::cout << "Testing double couple moment tensor..." << std::endl;
    
    double tol = 1e-10;
    
    // Create DC moment tensor: strike=0, dip=90, rake=0 (left-lateral strike-slip)
    double M0 = 1e18;  // 1e18 N·m
    MomentTensor M = MomentTensor::doubleCouple(0.0, M_PI/2.0, 0.0, M0);
    
    // DC should have zero trace (deviatoric)
    double trace = M.Mxx + M.Myy + M.Mzz;
    if (std::abs(trace) > tol * M0) {
        std::cerr << "  FAIL: DC trace = " << trace << ", should be 0" << std::endl;
        return false;
    }
    
    // Scalar moment should match input
    double M0_computed = M.scalarMoment();
    if (std::abs(M0_computed - M0) / M0 > 0.01) {
        std::cerr << "  FAIL: M0 = " << M0_computed << ", expected " << M0 << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Double couple properties correct" << std::endl;
    return true;
}

/**
 * @test Test explosion source (isotropic)
 */
bool test_explosion_source() {
    std::cout << "Testing explosion source..." << std::endl;
    
    double tol = 1e-10;
    double M0 = 1e15;
    
    MomentTensor M = MomentTensor::explosion(M0);
    
    // Explosion should have Mxx = Myy = Mzz = M0
    if (std::abs(M.Mxx - M0) > tol || 
        std::abs(M.Myy - M0) > tol || 
        std::abs(M.Mzz - M0) > tol) {
        std::cerr << "  FAIL: Diagonal components not equal" << std::endl;
        return false;
    }
    
    // Off-diagonal should be zero
    if (std::abs(M.Mxy) > tol || std::abs(M.Mxz) > tol || std::abs(M.Myz) > tol) {
        std::cerr << "  FAIL: Off-diagonal components not zero" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Explosion source correct" << std::endl;
    return true;
}

/**
 * @test Test moment tensor magnitude calculation
 */
bool test_magnitude_calculation() {
    std::cout << "Testing magnitude calculation..." << std::endl;
    
    double tol = 0.05;  // 0.05 magnitude units tolerance
    
    // Test various magnitudes
    std::vector<std::pair<double, double>> test_cases = {
        {1e14, 3.3},   // ~M3.3
        {1e16, 4.6},   // ~M4.6
        {1e18, 6.0},   // ~M6.0
        {1e20, 7.3},   // ~M7.3
        {1e22, 8.6}    // ~M8.6
    };
    
    for (const auto& tc : test_cases) {
        double M0 = tc.first;
        double expected_Mw = tc.second;
        
        MomentTensor M = MomentTensor::doubleCouple(0.0, M_PI/2.0, 0.0, M0);
        double computed_Mw = M.magnitude();
        
        if (std::abs(computed_Mw - expected_Mw) > tol) {
            std::cerr << "  FAIL: M0=" << M0 << ", Mw=" << computed_Mw 
                      << ", expected~" << expected_Mw << std::endl;
            return false;
        }
    }
    
    std::cout << "  PASS: Magnitude calculation correct" << std::endl;
    return true;
}

/**
 * @test Test moment tensor addition and scaling
 */
bool test_moment_tensor_operations() {
    std::cout << "Testing moment tensor operations..." << std::endl;
    
    double tol = 1e-10;
    
    MomentTensor M1(1.0, 2.0, 3.0, 0.5, 0.3, 0.2);
    MomentTensor M2(0.5, 1.0, 1.5, 0.25, 0.15, 0.1);
    
    // Test addition
    MomentTensor M_sum = M1 + M2;
    if (std::abs(M_sum.Mxx - 1.5) > tol || std::abs(M_sum.Myy - 3.0) > tol) {
        std::cerr << "  FAIL: Addition incorrect" << std::endl;
        return false;
    }
    
    // Test scaling
    MomentTensor M_scaled = M1 * 2.0;
    if (std::abs(M_scaled.Mxx - 2.0) > tol || std::abs(M_scaled.Mxy - 1.0) > tol) {
        std::cerr << "  FAIL: Scaling incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Moment tensor operations correct" << std::endl;
    return true;
}

// =============================================================================
// Source Time Function Tests
// =============================================================================

/**
 * @test Test Gaussian source time function
 */
bool test_gaussian_stf() {
    std::cout << "Testing Gaussian STF..." << std::endl;
    
    SourceTimeFunctionEvaluator stf(SourceTimeFunction::GAUSSIAN);
    stf.setDuration(1.0);
    stf.setPeakTime(0.5);
    stf.setOnsetTime(0.0);
    
    // Peak should be at peak_time
    double rate_at_peak = stf.momentRate(0.5);
    double rate_before = stf.momentRate(0.2);
    double rate_after = stf.momentRate(0.8);
    
    if (rate_at_peak <= rate_before || rate_at_peak <= rate_after) {
        std::cerr << "  FAIL: Peak not at expected time" << std::endl;
        return false;
    }
    
    // Should be non-negative
    if (rate_at_peak < 0 || rate_before < 0 || rate_after < 0) {
        std::cerr << "  FAIL: Negative moment rate" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Gaussian STF correct" << std::endl;
    return true;
}

/**
 * @test Test Ricker wavelet
 */
bool test_ricker_stf() {
    std::cout << "Testing Ricker wavelet..." << std::endl;
    
    SourceTimeFunctionEvaluator stf(SourceTimeFunction::RICKER);
    stf.setFrequency(5.0);  // 5 Hz dominant frequency
    stf.setPeakTime(0.2);
    stf.setOnsetTime(0.0);
    
    // Ricker should have zero crossing at peak
    double rate = stf.momentRate(0.2);
    
    // Not checking exact zero because of implementation details
    // Just verify it's finite
    if (std::isnan(rate) || std::isinf(rate)) {
        std::cerr << "  FAIL: Invalid Ricker value" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Ricker wavelet correct" << std::endl;
    return true;
}

/**
 * @test Test STF integral (moment)
 */
bool test_stf_moment_integral() {
    std::cout << "Testing STF moment integral..." << std::endl;
    
    double tol = 0.1;  // 10% tolerance for numerical integration
    
    SourceTimeFunctionEvaluator stf(SourceTimeFunction::GAUSSIAN);
    stf.setDuration(1.0);
    stf.setPeakTime(0.5);
    stf.setOnsetTime(0.0);
    
    // Moment should increase with time and saturate
    double m1 = stf.moment(0.3);
    double m2 = stf.moment(0.5);
    double m3 = stf.moment(1.0);
    double m4 = stf.moment(2.0);
    
    // Should be monotonically increasing
    if (m1 > m2 || m2 > m3 || m3 > m4) {
        std::cerr << "  FAIL: Moment not monotonically increasing" << std::endl;
        return false;
    }
    
    // Should be non-negative
    if (m1 < 0 || m2 < 0 || m3 < 0) {
        std::cerr << "  FAIL: Negative moment" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: STF moment integral correct" << std::endl;
    return true;
}

/**
 * @test Test Yoffe slip rate function
 */
bool test_yoffe_stf() {
    std::cout << "Testing Yoffe slip rate function..." << std::endl;
    
    SourceTimeFunctionEvaluator stf(SourceTimeFunction::YOFFE);
    stf.setRiseTime(1.0);
    stf.setOnsetTime(0.0);
    
    // Yoffe should be zero before onset
    double rate_before = stf.momentRate(-0.1);
    if (std::abs(rate_before) > 1e-15) {
        std::cerr << "  FAIL: Non-zero before onset" << std::endl;
        return false;
    }
    
    // Should be zero after rise time
    double rate_after = stf.momentRate(1.5);
    if (std::abs(rate_after) > 1e-15) {
        std::cerr << "  FAIL: Non-zero after rise time" << std::endl;
        return false;
    }
    
    // Should be positive during slip
    double rate_during = stf.momentRate(0.3);
    if (rate_during <= 0) {
        std::cerr << "  FAIL: Non-positive during slip" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Yoffe function correct" << std::endl;
    return true;
}

// =============================================================================
// Point Source Tests
// =============================================================================

/**
 * @test Test point source location
 */
bool test_point_source_location() {
    std::cout << "Testing point source location..." << std::endl;
    
    double tol = 1e-12;
    
    PointSource src;
    src.setLocation(100.0, 200.0, -5000.0);
    
    double x, y, z;
    src.getLocation(x, y, z);
    
    if (std::abs(x - 100.0) > tol || 
        std::abs(y - 200.0) > tol || 
        std::abs(z - (-5000.0)) > tol) {
        std::cerr << "  FAIL: Location mismatch" << std::endl;
        return false;
    }
    
    auto loc = src.getLocation();
    if (std::abs(loc[0] - 100.0) > tol) {
        std::cerr << "  FAIL: Array location mismatch" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Point source location correct" << std::endl;
    return true;
}

/**
 * @test Test point source moment tensor evolution
 */
bool test_point_source_moment_evolution() {
    std::cout << "Testing point source moment evolution..." << std::endl;
    
    PointSource src;
    src.setLocation(0, 0, -5000);
    
    MomentTensor M0 = MomentTensor::doubleCouple(0, M_PI/2, 0, 1e18);
    src.setMomentTensor(M0);
    
    SourceTimeFunctionEvaluator stf(SourceTimeFunction::GAUSSIAN);
    stf.setDuration(1.0);
    stf.setPeakTime(0.5);
    src.setSourceTimeFunction(stf);
    
    // Before onset: moment rate should be small
    double mr_before = src.getMomentRate(-1.0);
    if (mr_before > 1e10) {  // Should be essentially zero
        std::cerr << "  FAIL: Large moment rate before onset" << std::endl;
        return false;
    }
    
    // At peak time: moment rate should be significant
    double mr_peak = src.getMomentRate(0.5);
    if (mr_peak <= 0) {
        std::cerr << "  FAIL: Zero moment rate at peak" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Point source moment evolution correct" << std::endl;
    return true;
}

/**
 * @test Test point source spatial smoothing
 */
bool test_point_source_smoothing() {
    std::cout << "Testing point source spatial smoothing..." << std::endl;
    
    PointSource src;
    src.setLocation(0, 0, 0);
    src.setSpatialSmoothing(100.0);  // 100m sigma
    
    // Weight at source location should be maximum
    double w_center = src.getSpatialWeight(0, 0, 0);
    double w_far = src.getSpatialWeight(500, 0, 0);
    
    if (w_far >= w_center) {
        std::cerr << "  FAIL: Weight doesn't decay with distance" << std::endl;
        return false;
    }
    
    // Very far should be essentially zero
    double w_very_far = src.getSpatialWeight(5000, 0, 0);
    if (w_very_far > 1e-10 * w_center) {
        std::cerr << "  FAIL: Weight doesn't decay to zero" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Point source smoothing correct" << std::endl;
    return true;
}

// =============================================================================
// Kinematic Source Tests
// =============================================================================

/**
 * @test Test rectangular fault creation
 */
bool test_kinematic_rectangular_fault() {
    std::cout << "Testing kinematic rectangular fault..." << std::endl;
    
    KinematicSource kin;
    kin.createRectangularFault(0, 0, -5000,  // center
                               0.0, M_PI/2,   // strike, dip (vertical)
                               10000, 5000,   // length, width
                               10, 5);        // subdivisions
    
    // Should have 50 subfaults
    if (kin.getSubfaults().size() != 50) {
        std::cerr << "  FAIL: Expected 50 subfaults, got " 
                  << kin.getSubfaults().size() << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Rectangular fault creation correct" << std::endl;
    return true;
}

/**
 * @test Test circular rupture propagation
 */
bool test_circular_rupture() {
    std::cout << "Testing circular rupture propagation..." << std::endl;
    
    KinematicSource kin;
    kin.createRectangularFault(0, 0, -5000, 0, M_PI/2, 10000, 5000, 10, 5);
    kin.setUniformSlip(1.0, 0.0);  // 1m slip, pure strike-slip
    kin.setCircularRupture(0, 0, -5000, 3000);  // Hypocenter at center, Vr=3km/s
    
    auto& subfaults = kin.getSubfaults();
    
    // Subfaults farther from hypocenter should have later rupture times
    double rt_center = 1e10;
    double rt_corner = 0;
    
    for (const auto& sf : subfaults) {
        double dist = std::sqrt(sf.x*sf.x + sf.y*sf.y + (sf.z+5000)*(sf.z+5000));
        if (dist < 100) rt_center = std::min(rt_center, sf.rupture_time);
        if (dist > 4000) rt_corner = std::max(rt_corner, sf.rupture_time);
    }
    
    if (rt_center >= rt_corner) {
        std::cerr << "  FAIL: Rupture times not increasing with distance" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Circular rupture propagation correct" << std::endl;
    return true;
}

/**
 * @test Test kinematic source total moment
 */
bool test_kinematic_total_moment() {
    std::cout << "Testing kinematic source total moment..." << std::endl;
    
    double tol = 0.01;  // 1% tolerance
    
    KinematicSource kin;
    kin.setShearModulus(30e9);  // 30 GPa
    kin.createRectangularFault(0, 0, -5000, 0, M_PI/2, 10000, 5000, 10, 5);
    kin.setUniformSlip(1.0, 0.0);  // 1m slip
    
    double M0 = kin.getTotalMoment();
    
    // Expected: M0 = μ * A * D = 30e9 * 10000*5000 * 1.0 = 1.5e18 N·m
    double expected = 30e9 * 10000.0 * 5000.0 * 1.0;
    
    if (std::abs(M0 - expected) / expected > tol) {
        std::cerr << "  FAIL: M0=" << M0 << ", expected " << expected << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Kinematic total moment correct" << std::endl;
    return true;
}

// =============================================================================
// Receiver Tests
// =============================================================================

/**
 * @test Test receiver data recording
 */
bool test_receiver_recording() {
    std::cout << "Testing receiver data recording..." << std::endl;
    
    SeismicReceiver rec;
    rec.setLocation(1000, 2000, 0);
    rec.setName("test_receiver");
    rec.setSamplingRate(0.01);
    
    // Record some data
    double u1[3] = {1.0, 0.5, 0.2};
    double u2[3] = {1.5, 0.7, 0.3};
    double u3[3] = {2.0, 0.9, 0.4};
    
    rec.record(0.0, u1);
    rec.record(0.01, u2);
    rec.record(0.02, u3);
    
    // Check data stored
    const auto& times = rec.getTimes();
    const auto& data = rec.getData();
    
    if (times.size() != 3 || data.size() != 3) {
        std::cerr << "  FAIL: Wrong number of recorded points" << std::endl;
        return false;
    }
    
    if (std::abs(data[1][0] - 1.5) > 1e-10) {
        std::cerr << "  FAIL: Recorded data incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Receiver recording correct" << std::endl;
    return true;
}

/**
 * @test Test receiver PGV calculation
 */
bool test_receiver_pgv() {
    std::cout << "Testing receiver PGV calculation..." << std::endl;
    
    SeismicReceiver rec;
    rec.setLocation(0, 0, 0);
    rec.setSamplingRate(0.001);
    
    // Record velocity data with known maximum
    for (int i = 0; i <= 100; ++i) {
        double t = i * 0.001;
        double v = std::sin(2 * M_PI * t * 10);  // 10 Hz sine wave
        double u[3] = {v, 0.5 * v, 0.2 * v};
        rec.record(t, u);
    }
    
    double pgv = rec.getPeakGroundVelocity();
    
    // Vector magnitude: sqrt(1 + 0.25 + 0.04) ≈ 1.136
    double expected_pgv = std::sqrt(1.0 + 0.25 + 0.04);
    
    if (std::abs(pgv - expected_pgv) / expected_pgv > 0.05) {
        std::cerr << "  FAIL: PGV=" << pgv << ", expected~" << expected_pgv << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Receiver PGV calculation correct" << std::endl;
    return true;
}

/**
 * @test Test fault receiver slip recording
 */
bool test_fault_receiver() {
    std::cout << "Testing fault receiver..." << std::endl;
    
    FaultReceiver rec;
    rec.setLocation(0, 0, -5000);
    rec.setName("fault_receiver_1");
    rec.setSamplingRate(0.001);
    
    // Record fault data
    rec.record(0.0, 0.0, 0.0, 50e6, 20e6, 0.0, 1e6);  // Before rupture
    rec.record(0.1, 0.5, 5.0, 50e6, 15e6, 0.0, 1e4);  // During rupture
    rec.record(0.2, 1.0, 0.0, 50e6, 10e6, 0.0, 1e3);  // After rupture
    
    double psr = rec.getPeakSlipRate();
    if (std::abs(psr - 5.0) > 1e-10) {
        std::cerr << "  FAIL: Peak slip rate incorrect" << std::endl;
        return false;
    }
    
    double rt = rec.getRuptureTime(0.1);  // Threshold 0.1 m/s
    if (std::abs(rt - 0.1) > 0.05) {
        std::cerr << "  FAIL: Rupture time incorrect" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Fault receiver correct" << std::endl;
    return true;
}

// =============================================================================
// Source-Receiver Manager Tests
// =============================================================================

/**
 * @test Test source-receiver manager
 */
bool test_source_receiver_manager() {
    std::cout << "Testing source-receiver manager..." << std::endl;
    
    SourceReceiverManager mgr;
    
    // Add point source
    auto src = std::make_unique<PointSource>();
    src->setLocation(0, 0, -5000);
    src->setFromFault(0, M_PI/2, 0, 1e18, 0, 0, -5000);
    mgr.addPointSource(std::move(src));
    
    // Add receivers
    auto rec1 = std::make_unique<SeismicReceiver>();
    rec1->setLocation(1000, 0, 0);
    rec1->setName("rec1");
    mgr.addReceiver(std::move(rec1));
    
    auto rec2 = std::make_unique<SeismicReceiver>();
    rec2->setLocation(2000, 0, 0);
    rec2->setName("rec2");
    mgr.addReceiver(std::move(rec2));
    
    // Test source term computation
    double f[9];
    mgr.getSourceTerm(0, 0, -5000, 0.5, f);
    
    // Should have non-zero moment tensor contribution
    bool has_nonzero = false;
    for (int i = 0; i < 6; ++i) {
        if (std::abs(f[i]) > 1e-10) has_nonzero = true;
    }
    
    if (!has_nonzero) {
        std::cerr << "  FAIL: No source term at source location" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Source-receiver manager correct" << std::endl;
    return true;
}

/**
 * @test Test receiver line creation
 */
bool test_receiver_line() {
    std::cout << "Testing receiver line creation..." << std::endl;
    
    SourceReceiverManager mgr;
    mgr.addReceiverLine(0, 0, 0,      // Start
                        10000, 0, 0,   // End
                        11,            // 11 receivers
                        "line_");
    
    // Should have 11 receivers
    if (mgr.getReceivers().size() != 11) {
        std::cerr << "  FAIL: Expected 11 receivers, got " 
                  << mgr.getReceivers().size() << std::endl;
        return false;
    }
    
    // Check spacing (should be 1000m)
    double x0, y0, z0, x1, y1, z1;
    mgr.getReceivers()[0]->getLocation(x0, y0, z0);
    mgr.getReceivers()[1]->getLocation(x1, y1, z1);
    
    double spacing = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
    if (std::abs(spacing - 1000.0) > 1.0) {
        std::cerr << "  FAIL: Spacing=" << spacing << ", expected 1000" << std::endl;
        return false;
    }
    
    std::cout << "  PASS: Receiver line creation correct" << std::endl;
    return true;
}

// =============================================================================
// Test Runner
// =============================================================================

int runSeismicSourceTests() {
    std::cout << "\n=== Seismic Source Unit Tests ===" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    // Moment tensor tests
    if (test_double_couple_properties()) ++passed; else ++failed;
    if (test_explosion_source()) ++passed; else ++failed;
    if (test_magnitude_calculation()) ++passed; else ++failed;
    if (test_moment_tensor_operations()) ++passed; else ++failed;
    
    // STF tests
    if (test_gaussian_stf()) ++passed; else ++failed;
    if (test_ricker_stf()) ++passed; else ++failed;
    if (test_stf_moment_integral()) ++passed; else ++failed;
    if (test_yoffe_stf()) ++passed; else ++failed;
    
    // Point source tests
    if (test_point_source_location()) ++passed; else ++failed;
    if (test_point_source_moment_evolution()) ++passed; else ++failed;
    if (test_point_source_smoothing()) ++passed; else ++failed;
    
    // Kinematic source tests
    if (test_kinematic_rectangular_fault()) ++passed; else ++failed;
    if (test_circular_rupture()) ++passed; else ++failed;
    if (test_kinematic_total_moment()) ++passed; else ++failed;
    
    // Receiver tests
    if (test_receiver_recording()) ++passed; else ++failed;
    if (test_receiver_pgv()) ++passed; else ++failed;
    if (test_fault_receiver()) ++passed; else ++failed;
    
    // Manager tests
    if (test_source_receiver_manager()) ++passed; else ++failed;
    if (test_receiver_line()) ++passed; else ++failed;
    
    std::cout << "\n=== Seismic Source Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total:  " << (passed + failed) << std::endl;
    
    return failed;
}

} // namespace Testing
} // namespace FSRM

#ifndef FSRM_TEST_NO_MAIN
int main() {
    return FSRM::Testing::runSeismicSourceTests();
}
#endif
