#pragma once
#include "ks.h"
#include <algorithm>
#include <limits>
#include <type_traits>

struct KSIntegrationState {
    KSVector u{};
    KSVector uprime{};
    PS::F64 h = 0;
    PS::F64 time = 0;
    PS::F64 gm = 0; // G(m_1+m_2)
    KSVector fourth{};
    KSVector fifth{};
    PS::F64 hthird = 0;
    PS::F64 hfourth = 0;
    bool isPredictionCached = false;
    PS::F64vec pert{0};
    PS::F64vec pert_jerk{0};
    bool isPertCached = false;
};

/*
 * PERTURBATION EVALUATION
 */

// Boolean callbacks for KS perturbation evaluation
enum class KSForceResult {
    success,
    retry,
    failure
};

// Result of perturbation evaluation
template<class Result>
inline KSForceResult classifyKSForceResult(const Result result)
{
    if constexpr (std::is_same<Result, KSForceResult>::value) return result;
    else {
        static_assert(std::is_same<Result, bool>::value, "KS force must return bool or KSForceResult");
        return result ? KSForceResult::success : KSForceResult::failure;
    }
}

// Evaluate a KS perturbation force
template<class Force>
class KSForceAttempt {
public:
    Force& force;
    PS::F64 accepted_time = 0;
    KSForceResult failure = KSForceResult::success;

    void reset(const PS::F64 time)
    {
        accepted_time = time;
        failure = KSForceResult::success;
    }
    bool operator()(const KSIntegrationState& state, PS::F64vec& p, PS::F64vec& j)
    {
        const auto result = classifyKSForceResult(force(state, p, j));
        if (result != KSForceResult::success) {
            failure = state.time > accepted_time ? result : KSForceResult::failure;
            return false;
        }
        return true;
    }
    bool retryable() const { return failure == KSForceResult::retry; }
};

/*
 * INTEGRATOR
 */

struct KSDerivatives {
    KSVector acc{}; // u''
    KSVector jerk{}; // u'''
    PS::F64 hprime = 0; // h'
    PS::F64 hsecond = 0; // h''
    PS::F64 r = 0;
    PS::F64 rprime = 0;
};

struct KSIntegrationStatistics {
    PS::U64 step_calls = 0;
    PS::U64 correction_iterations = 0;
    PS::U64 force_evaluations = 0;
    PS::U64 accepted_steps = 0;
    PS::U64 rejected_steps = 0;
    PS::U64 reject_iteration = 0;
    PS::U64 reject_target = 0;
    PS::U64 reject_error = 0;
};

struct KSIntegratorParameters {
    KSIntegrationStatistics* statistics = nullptr;
    PS::F64 iteration_tolerance = 2e-13;
    PS::F64 time_tolerance = 2e-14;
    PS::F64 ma98_eta = 0.2;
    PS::S32 ma98_corrections = 2;
    bool high_order_predictor = true;
    PS::S32 max_iterations = 24;
    PS::S32 max_steps = 100000;
};

// initialize function of KS integrator
inline KSIntegrationState initializeKSBinary(
    const PS::F64vec& x,
    const PS::F64vec& v,
    const PS::F64 gm,
    const PS::F64 time)
{
    KSIntegrationState state;
    posToKS(x, state.u.data());
    velocityToKS(v, state.u.data(), state.uprime.data());
    state.h = calcKSSpecificEnergy(state.u.data(), state.uprime.data(), gm);
    state.time = time;
    state.gm = gm;
    return state;
}

// check u, u', h, time
inline bool isFiniteKSState(const KSIntegrationState& state)
{
    for (int k = 0; k < 4; ++k)
        if (!std::isfinite(state.u[k]) || !std::isfinite(state.uprime[k])) return false;
    return std::isfinite(state.h) && std::isfinite(state.time);
}

// Calculate the difference between two KS states.
// The values are normalized
// use for iteration
inline PS::F64 differenceKSState(const KSIntegrationState& a, const KSIntegrationState& b)
{
    PS::F64 error = 0;
    const PS::F64 u_scale = std::max(1e-100,
        std::sqrt(std::max(calcDotKSVector(a.u, a.u), calcDotKSVector(b.u, b.u))));
    const PS::F64 w_scale = std::max(1e-100,
        std::sqrt(std::max(calcDotKSVector(a.uprime, a.uprime), calcDotKSVector(b.uprime, b.uprime))));
    PS::F64 u_difference = 0, w_difference = 0;
    for (int k = 0; k < 4; ++k) {
        u_difference = std::max(u_difference, std::abs(a.u[k]-b.u[k]));
        w_difference = std::max(w_difference, std::abs(a.uprime[k]-b.uprime[k]));
    }
    // Each positive normalization is shared by all four components.
    error = std::max(u_difference/u_scale, w_difference/w_scale);
    error = std::max(error, std::abs(a.h-b.h)/std::max(1e-100, std::max(std::abs(a.h), std::abs(b.h))));
    error = std::max(error, std::abs(a.time-b.time)/(1+std::max(std::abs(a.time), std::abs(b.time))));
    return error;
}

// Calculate perturbation and its time derivative
// then calculate u'', u''', h', h'', r', and r,''
// If the solution is not finite, return false
template<class Force>
inline bool evaluateKSDerivatives(
    const KSIntegrationState& state,
    Force& force,
    KSDerivatives& derivatives)
{
    if (!isFiniteKSState(state)) return false;
    PS::F64vec p(0), j(0);
    if (classifyKSForceResult(force(state, p, j)) != KSForceResult::success) return false;
    calcKSBinaryHigherDerivatives(state.u.data(), state.uprime.data(), state.h, p, j,
        derivatives.acc.data(), derivatives.jerk.data(), derivatives.hprime,
        derivatives.hsecond, derivatives.r, derivatives.rprime);
    for (int k = 0; k < 4; ++k)
        if (!std::isfinite(derivatives.acc[k]) || !std::isfinite(derivatives.jerk[k])) return false;
    return derivatives.r > 0 && std::isfinite(derivatives.r)
        && std::isfinite(derivatives.hprime) && std::isfinite(derivatives.hsecond);
}

// Calculate Stumpff functions c0 through c5.
// see Mikkola (2020), section 2.6
inline bool calcKSStumpffFunctions(const PS::F64 z, PS::F64 c[6])
{
    if (!std::isfinite(z) || std::abs(z) > 4) return false;
    PS::F64 h = -z;
    PS::S32 i = 0;
    while (std::abs(h) >= 0.9) {
        h *= 0.25;
        ++i;
    }
    c[4] = (201859257600.0 + h*(-3741257520.0
        + h*(40025040.0 - 147173.0*h)))
        / (240.0*(20185925760.0 + h*(298738440.0
        + h*(1945020.0 + 5801.0*h))));
    c[5] = (3750361655040.0 + h*(-40967886960.0
        + h*(358614256.0 - 1029037.0*h)))
        / (55440.0*(8117665920.0 + h*(104602680.0
        + h*(582348.0 + 1451.0*h))));
    for (PS::S32 k = 0; k < i; ++k) {
        const PS::F64 c3 = 1.0/6.0 - h*c[5];
        const PS::F64 c2 = 0.5 - h*c[4];
        c[5] = (c[5] + c[4] + c2*c3)/16.0;
        c[4] = c3*(2.0 - h*c3)/8.0;
        h *= 4.0;
    }
    c[3] = 1.0/6.0 + z*c[5];
    c[2] = 0.5 + z*c[4];
    c[1] = 1.0 + z*c[3];
    c[0] = 1.0 + z*c[2];
    return true;
}

/*
 * ANALYTICAL KEPLER SOLUTION
 */

struct KSKeplerStep {
    KSIntegrationState base;
    PS::F64 lambda = 0;
    PS::F64 c[6]{};
};

// calculate analytical solution of Kepler problem
// if no perturbation, u''= lambda*u (harmonic oscillator)
inline bool prepareKSKeplerStep(
    const KSIntegrationState& initial,
    const PS::F64 ds,
    KSKeplerStep& step)
{
    step.lambda = initial.h/2;
    const PS::F64 z = step.lambda*ds*ds;
    PS::F64 doubled[6];
    if (!calcKSStumpffFunctions(z, step.c)
        || !calcKSStumpffFunctions(4*z, doubled)) return false;
    step.base = initial;
    for (int k = 0; k < 4; ++k) {
        step.base.u[k] = step.c[0]*initial.u[k]
            +ds*step.c[1]*initial.uprime[k];
        step.base.uprime[k] = step.lambda*ds*step.c[1]*initial.u[k]
            +step.c[0]*initial.uprime[k];
    }
    // calculate real time using integral of |u_Kepler(s)|^2
    step.base.time
        += ds*(1+doubled[1])*calcDotKSVector(initial.u, initial.u)/2
        +ds*ds*step.c[1]*step.c[1]*calcDotKSVector(initial.u, initial.uprime)
        +2*ds*ds*ds*doubled[3]*calcDotKSVector(initial.uprime, initial.uprime);
    return isFiniteKSState(step.base);
}

// Predict the endpoint using the initial perturbation derivatives.
inline KSIntegrationState predictKSKepler(
    const KSIntegrationState& initial,
    const KSDerivatives& a,
    const PS::F64 ds,
    const KSKeplerStep& step)
{
    auto result = step.base;
    for (int k = 0; k < 4; ++k) {
        const PS::F64 f0 = a.acc[k]-step.lambda*initial.u[k];
        const PS::F64 d0 = ds*(a.jerk[k]-step.lambda*initial.uprime[k]);
        result.u[k] += ds*ds*(f0*step.c[2]+d0*step.c[3]);
        result.uprime[k] += ds*(f0*step.c[1]+d0*step.c[2]);
    }
    result.h += ds*(a.hprime+ds*a.hsecond/2);
    return result;
}


/*
 * PREDICTOR AND CORRECTOR WITH PERTURBATION
 */

// Predictor
inline KSIntegrationState predictKSKeplerHighOrder(
    const KSIntegrationState& initial,
    const KSDerivatives& a,
    const PS::F64 ds,
    const KSKeplerStep& step)
{
    auto result = step.base;
    const PS::F64 ds2 = ds*ds, ds3 = ds2*ds;
    for (int k = 0; k < 4; ++k) {
        result.u[k]
            = initial.u[k]+ds*initial.uprime[k]
            +ds2*a.acc[k]/2+ds3*a.jerk[k]/6
            +ds2*ds2*step.c[4]*initial.fourth[k]
            +ds3*ds2*step.c[5]*initial.fifth[k];
        result.uprime[k]
            = initial.uprime[k]+ds*a.acc[k]
            +ds2*a.jerk[k]/2+ds3*step.c[3]*initial.fourth[k]
            +ds2*ds2*step.c[4]*initial.fifth[k];
    }
    result.h
        = initial.h +ds*a.hprime +ds2*a.hsecond/2
        +ds3*initial.hthird/6 +ds2*ds2*initial.hfourth/24;
    // Predict physical time including perturbations before evaluating the force.
    KSVector base_a{}, base_j{}, base_fourth{};
    for (int k = 0; k < 4; ++k) {
        base_a[k] = step.lambda*initial.u[k];
        base_j[k] = step.lambda*initial.uprime[k];
        base_fourth[k] = step.lambda*base_a[k];
    }
    PS::F64 delta_r2 = 0, delta_r3 = 0, delta_r4 = 0;
    for (int k = 0; k < 4; ++k) {
        delta_r2 += 2*initial.u[k]*(a.acc[k]-base_a[k]);
        delta_r3 += 6*initial.uprime[k]*(a.acc[k]-base_a[k])
            +2*initial.u[k]*(a.jerk[k]-base_j[k]);
        delta_r4 += 6*(a.acc[k]*a.acc[k]-base_a[k]*base_a[k])
            +8*initial.uprime[k]*(a.jerk[k]-base_j[k])
            +2*initial.u[k]*(initial.fourth[k]-base_fourth[k]);
    }
    result.time += ds3*delta_r2/6+ds2*ds2*delta_r3/24+ds3*ds2*delta_r4/120;
    return result;
}

// Corrector
inline KSIntegrationState correctKSKepler(
    const KSIntegrationState& initial,
    const KSIntegrationState& trial,
    const KSDerivatives& a,
    const KSDerivatives& b,
    const PS::F64 ds,
    const KSKeplerStep& step)
{
    auto result = step.base;
    for (int k = 0; k < 4; ++k) {
        const PS::F64 f0 = a.acc[k]-step.lambda*initial.u[k];
        const PS::F64 d0 = ds*(a.jerk[k]-step.lambda*initial.uprime[k]);
        const PS::F64 f1 = b.acc[k]-step.lambda*trial.u[k];
        const PS::F64 d1 = ds*(b.jerk[k]-step.lambda*trial.uprime[k]);
        const PS::F64 q2 = 3*(f1-f0)-2*d0-d1;
        const PS::F64 q3 = 2*(f0-f1)+d0+d1;
        result.u[k] += ds*ds*(f0*step.c[2]+d0*step.c[3]
            +2*q2*step.c[4]+6*q3*step.c[5]);
        result.uprime[k] += ds*(f0*step.c[1]+d0*step.c[2]
            +2*q2*step.c[3]+6*q3*step.c[4]);
    }
    result.h += ds*(a.hprime+b.hprime)/2-ds*ds*(b.hsecond-a.hsecond)/12;
    const PS::F64 base_r = calcDotKSVector(step.base.u, step.base.u);
    const PS::F64 base_rprime = 2*calcDotKSVector(step.base.u, step.base.uprime);
    result.time += ds*(b.r-base_r)/2-ds*ds*(b.rprime-base_rprime)/12;
    return result;
}

/*
 * INTEGRATION
 */

// Mikkola & Aarseth (1998), section 4.1; our h=-2 Omega and C_n(z)=c_n(-z).
// Valid only for this initial state, initial derivatives, and trial ds.
// Corrector iterations change fourth, but not these coefficients.
struct KSClockCoefficients {
    PS::F64 r2 = 0; // r''
    PS::F64 r3 = 0; // r'''
    PS::F64 r4_base = 0; // first two terms of r'''' (unchange while PC loop)
    PS::F64 c4 = 0; // c_4
    PS::F64 c5 = 0; // c_5
};

inline bool prepareKSClock(
    const KSIntegrationState& s,
    const KSDerivatives& a,
    const PS::F64 ds,
    KSClockCoefficients& clock)
{
    PS::F64 c[6];
    // Preserve the original argument evaluation order; prepareKSKeplerStep
    // forms its doubled argument through a different expression.
    if (!calcKSStumpffFunctions(2*s.h*ds*ds, c)) return false;
    clock.c4 = c[4];
    clock.c5 = c[5];
    clock.r2 = 2*(calcDotKSVector(s.uprime, s.uprime)
        +calcDotKSVector(s.u, a.acc));
    clock.r3 = 6*calcDotKSVector(s.uprime, a.acc)
        +2*calcDotKSVector(s.u, a.jerk);
    clock.r4_base = 6*calcDotKSVector(a.acc, a.acc)
        +8*calcDotKSVector(s.uprime, a.jerk);
    return true;
}

inline PS::F64 calcKSClock(
    const KSIntegrationState& s,
    const KSDerivatives& a,
    const KSVector& fourth,
    const PS::F64 ds,
    const KSClockCoefficients& clock)
{
    const PS::F64 r4 = clock.r4_base+2*calcDotKSVector(s.u, fourth);
    // Eq. (23) in MA98, n=5: t^(k)=r^(k-1); the last two terms use c_4, c_5 at 4z.
    return s.time+ds*(a.r+ds*(a.rprime/2+ds*(clock.r2/6
        +ds*(clock.c4*clock.r3+ds*clock.c5*r4))));
}

// propose a final step from the MA98 predictor clock.
inline PS::F64 proposeKSTargetSeed(const KSIntegrationState& s,
    const KSDerivatives& a, const PS::F64 target, const PS::F64 phase_limit,
    const PS::F64 original, const PS::F64 tolerance)
{
    KSVector fourth;
    for (int k = 0; k < 4; ++k)
        fourth[k] = s.isPredictionCached ? s.fourth[k] : .5*s.h*a.acc[k];
    KSClockCoefficients coefficients;
    if (!prepareKSClock(s,a,phase_limit,coefficients)) return original;
    const PS::F64 r4 = coefficients.r4_base+2*calcDotKSVector(s.u,fourth);
    const auto evaluate = [&](PS::F64 ds, PS::F64& residual, PS::F64& derivative) {
        PS::F64 c[6];
        if (!calcKSStumpffFunctions(2*s.h*ds*ds,c)) return false;
        residual = s.time+ds*(a.r+ds*(a.rprime/2+ds*(coefficients.r2/6
            +ds*(c[4]*coefficients.r3+ds*c[5]*r4))))-target;
        derivative = a.r+ds*(a.rprime+ds*(coefficients.r2/2
            +ds*(c[3]*coefficients.r3+ds*c[4]*r4)));
        return std::isfinite(residual) && std::isfinite(derivative) && derivative > 0;
    };
    PS::F64 residual, derivative;
    if (!evaluate(phase_limit,residual,derivative) || residual < 0) return original;
    PS::F64 lower = 0, upper = phase_limit, ds = original;
    for (int iteration = 0; iteration < 8; ++iteration) {
        if (!evaluate(ds,residual,derivative)) return original;
        if (std::abs(residual) <= tolerance) return ds;
        if (residual > 0) upper = ds; else lower = ds;
        PS::F64 next = ds-residual/derivative;
        if (!(next > lower && next < upper)) next = .5*(lower+upper);
        ds = next;
    }
    return ds;
}

// Integrate one step in KS time.
template<class Force>
inline bool integrateKSStep(const KSIntegrationState& initial, const PS::F64 ds,
    Force& force, KSIntegrationState& output, const KSIntegratorParameters& parameters = {},
    bool* evaluation_failed = nullptr)
{
    if (evaluation_failed) *evaluation_failed = false;
    if (parameters.statistics) ++parameters.statistics->step_calls;
    if (!(ds > 0) || !std::isfinite(ds) || parameters.ma98_corrections < 1) return false;
    auto physical = [&](const KSIntegrationState& s, PS::F64vec& p, PS::F64vec& j) {
        if (parameters.statistics) ++parameters.statistics->force_evaluations;
        if (classifyKSForceResult(force(s, p, j)) != KSForceResult::success) {
            if (evaluation_failed) *evaluation_failed = true;
            return false;
        }
        for (int k = 0; k < 3; ++k) {
            if (!std::isfinite(p[k]) || !std::isfinite(j[k])) {
                if (evaluation_failed) *evaluation_failed = true;
                return false;
            }
        }
        return true;
    };
    PS::F64vec p0(0), j0(0);
    if (initial.isPertCached) {
        p0 = initial.pert;
        j0 = initial.pert_jerk;
    } else if (!physical(initial, p0, j0)) return false;
    auto start_force = [&](const KSIntegrationState&, PS::F64vec& p, PS::F64vec& j) {
        p = p0; j = j0; return true;
    };
    KSDerivatives a, b;
    if (!evaluateKSDerivatives(initial, start_force, a)) return false;
    KSKeplerStep step;
    if (!prepareKSKeplerStep(initial, ds, step)) return false;
    KSClockCoefficients clock;
    if (!prepareKSClock(initial, a, ds, clock)) return false;
    auto trial = initial.isPredictionCached
        ? predictKSKeplerHighOrder(initial, a, ds, step)
        : predictKSKepler(initial, a, ds, step);
    KSVector u4{}, u5{};
    for (int k = 0; k < 4; ++k) {
        u4[k] = initial.isPredictionCached ? initial.fourth[k]
            : step.lambda*a.acc[k];
        u5[k] = initial.isPredictionCached ? initial.fifth[k]
            : step.lambda*a.jerk[k];
    }
    trial.time = calcKSClock(initial, a, u4, ds, clock);
    PS::F64vec frozen_p(0), frozen_j(0);
    if (!physical(trial, frozen_p, frozen_j)) return false;
    auto endpoint_force = [&](const KSIntegrationState&, PS::F64vec& p, PS::F64vec& j) {
        p = frozen_p; j = frozen_j; return true;
    };
    const int count = initial.isPredictionCached
        ? parameters.ma98_corrections : parameters.max_iterations;
    for (int iteration = 0; iteration < count; ++iteration) {
        if (parameters.statistics) ++parameters.statistics->correction_iterations;
        if (iteration && !initial.isPredictionCached
            && !physical(trial, frozen_p, frozen_j)) return false;
        if (!evaluateKSDerivatives(trial, endpoint_force, b)) return false;
        // Eqs. (33)-(36): Hermite reconstruction of the cubic residual.
        for (int k = 0; k < 4; ++k) {
            const PS::F64 f0 = a.acc[k]-step.lambda*initial.u[k];
            const PS::F64 f1 = b.acc[k]-step.lambda*trial.u[k];
            const PS::F64 d0 = ds*(a.jerk[k]-step.lambda*initial.uprime[k]);
            const PS::F64 d1 = ds*(b.jerk[k]-step.lambda*trial.uprime[k]);
            u4[k] = step.lambda*a.acc[k]
                +2*(3*(f1-f0)-2*d0-d1)/(ds*ds);
            u5[k] = step.lambda*a.jerk[k]
                +6*(2*(f0-f1)+d0+d1)/(ds*ds*ds);
        }
        auto next = correctKSKepler(initial, trial, a, b, ds, step);
        next.time = calcKSClock(initial, a, u4, ds, clock);
        if (!isFiniteKSState(next) || !(next.time > initial.time)) return false;
        const bool done = initial.isPredictionCached ? iteration+1 == count
            : differenceKSState(next, trial) <= parameters.iteration_tolerance;
        if (done) {
            // Eqs. (38)-(39), rather than differentiating the endpoint residual.
            for (int k = 0; k < 4; ++k) {
                next.fourth[k] = u4[k]+ds*u5[k];
                next.fifth[k] = u5[k];
            }
            next.hfourth = (12*(a.hprime-b.hprime)
                +6*ds*(a.hsecond+b.hsecond))/(ds*ds*ds);
            const PS::F64 h3 = (6*(b.hprime-a.hprime)
                -ds*(4*a.hsecond+2*b.hsecond))/(ds*ds);
            next.hthird = h3+ds*next.hfourth;
            next.isPredictionCached = true;
            next.pert = frozen_p;
            next.pert_jerk = frozen_j;
            next.isPertCached = true;
            // GPLUM-specific admission/contact guard; does not change the
            // frozen perturbation saved for the next predictor/corrector.
            PS::F64vec checked_p(0), checked_j(0);
            if (!physical(next, checked_p, checked_j)) return false;
            output = next;
            return true;
        }
        trial = next;
    }
    return false;
}

// Integrate the binary to the target physical time.
template<class Force>
inline bool integrateKSBinaryToTime(const KSIntegrationState& initial, const PS::F64 target,
    Force& force, KSIntegrationState& output, const KSIntegratorParameters& parameters = {})
{
    if (!std::isfinite(target) || target < initial.time || !(parameters.ma98_eta > 0)
        || !std::isfinite(parameters.ma98_eta) || !(parameters.time_tolerance > 0)) return false;
    KSForceAttempt<Force> attempted_force{force};
    auto state = initial;
    // A new hard coupling interval may change the external force interpolant.
    state.isPredictionCached = false;
    state.isPertCached = false;
    const PS::F64 tolerance = std::max(parameters.time_tolerance*(target-initial.time),
        16*std::numeric_limits<PS::F64>::epsilon()*std::max(std::abs(target), std::abs(initial.time)));
    for (int n = 0; n < parameters.max_steps; ++n) {
        if (std::abs(state.time-target) <= tolerance) { output = state; return true; }
        const PS::F64 r = calcDotKSVector(state.u, state.u);
        const PS::F64 gm = state.gm;
        if (!(r > 0) || !(gm > 0)) return false;
        if (!state.isPertCached) {
            if (parameters.statistics) ++parameters.statistics->force_evaluations;
            if (classifyKSForceResult(force(state, state.pert, state.pert_jerk)) != KSForceResult::success) return false;
            for (int k = 0; k < 3; ++k)
                if (!std::isfinite(state.pert[k]) || !std::isfinite(state.pert_jerk[k])) return false;
            state.isPertCached = true;
        }
        const PS::F64 gamma = std::sqrt(state.pert*state.pert)*r*r/gm;
        // Eq. (30). The parabolic limit needs a local crossing-time cap,
        // since the paper's energy-based formula diverges at h=0.
        PS::F64 ds = parameters.ma98_eta/std::sqrt(2*std::max(std::abs(state.h), 1e-100))
            /std::cbrt(1+1000*gamma);
        if (state.h >= 0) ds = std::min(ds, parameters.ma98_eta*std::sqrt(r/gm));
        const PS::F64 phase_limit = ds;
        ds = std::min(ds, (target-state.time)/r);
        auto seed_force = [&](const KSIntegrationState&, PS::F64vec& p, PS::F64vec& j) {
            p=state.pert; j=state.pert_jerk; return true;
        };
        KSDerivatives seed_derivatives;
        if (evaluateKSDerivatives(state,seed_force,seed_derivatives))
            ds=proposeKSTargetSeed(state,seed_derivatives,target,phase_limit,ds,tolerance);
        PS::F64 lower = 0, upper = 0;
        bool accepted = false;
        for (int retry = 0; retry < 60; ++retry) {
            KSIntegrationState next;
            bool force_failed = false;
            attempted_force.reset(state.time);
            if (!integrateKSStep(state, ds, attempted_force, next, parameters, &force_failed)) {
                if (force_failed && !attempted_force.retryable()) return false;
                if (parameters.statistics) {
                    ++parameters.statistics->rejected_steps;
                    ++parameters.statistics->reject_iteration;
                }
                ds *= 0.5; lower = upper = 0; continue;
            }
            const PS::F64 residual = next.time-target;
            if (residual > tolerance || (upper > 0 && residual < -tolerance)) {
                if (residual > 0) upper = ds; else lower = ds;
                PS::F64 proposed = ds-residual/calcDotKSVector(next.u, next.u);
                if (!(proposed > lower && proposed < upper)) proposed = (lower+upper)/2;
                ds = proposed;
                if (parameters.statistics) {
                    ++parameters.statistics->rejected_steps;
                    ++parameters.statistics->reject_target;
                }
                continue;
            }
            state = next;
            if (parameters.statistics) ++parameters.statistics->accepted_steps;
            accepted = true;
            break;
        }
        if (!accepted) return false;
    }
    return false;
}
