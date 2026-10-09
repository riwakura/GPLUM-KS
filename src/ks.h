#pragma once
#include <particle_simulator.hpp>
#include <array>
#include <cmath>
#include <stdexcept>

// KS transformations, binary equations of motion, and equations of motion.
// KS vector is 4th dementional vector u = (u1, u2, u3, u4).
using KSVector = std::array<PS::F64, 4>;

// x=L(u)u; dt/ds=r=|u|^2; uprime denotes du/ds.
inline PS::F64 calcDotKSVector(const KSVector& a, const KSVector& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

inline PS::F64vec calcCrossKSVector(const PS::F64vec& a, const PS::F64vec& b)
{
    return PS::F64vec(a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x);
}

inline PS::F64 calcNormOfKSVector(const PS::F64 u[4])
{
    return std::sqrt(u[0] * u[0] + u[1] * u[1]
        + u[2] * u[2] + u[3] * u[3]);
}

inline PS::F64vec ksToPos(const PS::F64 u[4])
{
    return PS::F64vec(
        u[0] * u[0] - u[1] * u[1] - u[2] * u[2] + u[3] * u[3],
        2 * (u[0] * u[1] - u[2] * u[3]),
        2 * (u[0] * u[2] + u[1] * u[3]));
}

inline PS::F64vec ksToVel(const PS::F64 u[4], const PS::F64 uprime[4])
{
    const PS::F64 r = u[0] * u[0] + u[1] * u[1]
        + u[2] * u[2] + u[3] * u[3];
    if (!(r > 0) || !std::isfinite(r))
        throw std::domain_error("KS: invalid radius");
    const PS::F64 r_inv = 1 / r;
    return PS::F64vec(
        2 * (u[0] * uprime[0] - u[1] * uprime[1] - u[2] * uprime[2] + u[3] * uprime[3]) * r_inv,
        2 * (u[1] * uprime[0] + u[0] * uprime[1] - u[3] * uprime[2] - u[2] * uprime[3]) * r_inv,
        2 * (u[2] * uprime[0] + u[3] * uprime[1] + u[0] * uprime[2] + u[1] * uprime[3]) * r_inv);
}

// impose KS constraint, this value should be zero.
// see for eq 2.144 Tremaine (2023), p125
inline PS::F64 calcKSGaugeResidual(const PS::F64 u[4], const PS::F64 uprime[4])
{
    return u[3] * uprime[0] - u[0] * uprime[3] + u[1] * uprime[2] - u[2] * uprime[1];
}

// Convert Cartesian position (3D) to KS position u (4D).
inline void posToKS(const PS::F64vec& x, PS::F64 u[4])
{
    const PS::F64 r = std::sqrt(x.x * x.x + x.y * x.y + x.z * x.z);
    if (!(r > 0) || !std::isfinite(r))
        throw std::domain_error("KS: invalid Cartesian position");
    if (x.x >= 0) {
        u[0] = std::sqrt(0.5 * r + 0.5 * x.x);
        u[1] = x.y / (2 * u[0]);
        u[2] = x.z / (2 * u[0]);
        u[3] = 0;
    } else {
        u[1] = std::sqrt(0.5 * r - 0.5 * x.x);
        u[0] = x.y / (2 * u[1]);
        u[2] = 0;
        u[3] = x.z / (2 * u[1]);
    }
}

// Convert KS position vector uprime (4D) to cartesian velocity v (3D)
// uprime = du/ds = (1/2) L(u)^T v, with dt/ds = |u|^2.
inline void velocityToKS(const PS::F64vec& v, const PS::F64 u[4], PS::F64 uprime[4])
{
    const PS::F64 norm = calcNormOfKSVector(u);
    if (!(norm > 0) || !std::isfinite(norm))
        throw std::domain_error("KS: invalid KS position (u < 0 or u is infty).");
    for (PS::S32 k = 0; k < 3; ++k)
        if (!std::isfinite(v[k]))
            throw std::domain_error("KS: invalid velocity (uprime is infty)");
    // Cache the input so even aliased input/output buffers remain well-defined.
    const PS::F64 u0 = u[0], u1 = u[1], u2 = u[2], u3 = u[3];
    uprime[0] = 0.5 * (u0 * v.x + u1 * v.y + u2 * v.z);
    uprime[1] = 0.5 * (-u1 * v.x + u0 * v.y + u3 * v.z);
    uprime[2] = 0.5 * (-u2 * v.x - u3 * v.y + u0 * v.z);
    uprime[3] = 0.5 * (u3 * v.x - u2 * v.y + u1 * v.z);
}

// Specific Kepler energy per unit reduced mass (not total hard energy).
// Initialize h with this function; subsequently evolve h as an independent variable.
inline PS::F64 calcKSSpecificEnergy(
    const PS::F64 u[4],
    const PS::F64 uprime[4],
    const PS::F64 gm)
{
    const PS::F64 r = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
    if (!(r > 0) || !std::isfinite(r) || !(gm > 0) || !std::isfinite(gm))
        throw std::domain_error("KS: invalid radius or G(m1+m2)");
    const PS::F64 w2 = uprime[0]*uprime[0] + uprime[1]*uprime[1]
        + uprime[2]*uprime[2] + uprime[3]*uprime[3];
    return (2 * w2 - gm) / r;
}

// calculate the KS equation of motion
inline void calcKSBinaryDerivatives(
    const PS::F64 u[4],
    const PS::F64 uprime[4],
    const PS::F64 h,
    const PS::F64vec& pert,
    PS::F64 usecond[4],
    PS::F64& hprime,
    PS::F64& tprime)
{
    const PS::F64 r = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
    const PS::F64 lp[4] = {
        u[0]*pert.x + u[1]*pert.y + u[2]*pert.z,
        -u[1]*pert.x + u[0]*pert.y + u[3]*pert.z,
        -u[2]*pert.x - u[3]*pert.y + u[0]*pert.z,
        u[3]*pert.x - u[2]*pert.y + u[1]*pert.z
    };
    const PS::F64 dh = 2 * (uprime[0]*lp[0] + uprime[1]*lp[1]
        + uprime[2]*lp[2] + uprime[3]*lp[3]);
    for (PS::S32 k = 0; k < 4; ++k)
        usecond[k] = 0.5 * h * u[k] + 0.5 * r * lp[k];
    hprime = dh;
    tprime = r;
}

// Calculate higher s derivatives for fourth-order Hermite Predictor-Corrector iteration.
inline void calcKSBinaryHigherDerivatives(
    const PS::F64 u[4],
    const PS::F64 uprime[4],
    const PS::F64 h,
    const PS::F64vec& pert,
    const PS::F64vec& pert_jerk,
    PS::F64 usecond[4],
    PS::F64 uthird[4],
    PS::F64& hprime,
    PS::F64& hsecond,
    PS::F64& r,
    PS::F64& rprime)
{
    calcKSBinaryDerivatives(u, uprime, h, pert, usecond, hprime, r);
    rprime = 0;
    for (PS::S32 k = 0; k < 4; ++k) {
        rprime += 2 * u[k] * uprime[k];
    }
    const auto transpose = [](const PS::F64 q[4],
                              const PS::F64vec& p,
                              PS::F64 out[4]) {
        out[0] = q[0]*p.x + q[1]*p.y + q[2]*p.z;
        out[1] = -q[1]*p.x + q[0]*p.y + q[3]*p.z;
        out[2] = -q[2]*p.x - q[3]*p.y + q[0]*p.z;
        out[3] = q[3]*p.x - q[2]*p.y + q[1]*p.z;
    };
    PS::F64 q[4], dq[4], lj[4];
    transpose(u, pert, q);
    transpose(uprime, pert, dq);
    transpose(u, pert_jerk, lj);
    hsecond = 0;
    for (PS::S32 k = 0; k < 4; ++k) {
        dq[k] += r * lj[k];
        uthird[k] = 0.5 * (hprime*u[k] + h*uprime[k] + rprime*q[k] + r*dq[k]);
        hsecond += 2 * (usecond[k]*q[k] + uprime[k]*dq[k]);
    }
}
