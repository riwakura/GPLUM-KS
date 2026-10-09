#pragma once
#include <particle_simulator.hpp>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>

struct BinarySearchParameters {
    PS::F64 central_mass = FP_t::m_sun;
    PS::F64 hill_factor = 1;
    PS::F64 gravitational_constant = 1;
    PS::F64 inner_cutoff_ratio = FP_t::gamma;
    bool include_unbound = false;
};

struct BinaryInfo {
    std::size_t first = 0;
    std::size_t second = 0;
    PS::S64 id1 = 0;
    PS::S64 id2 = 0;
    PS::F64 separation = 0;
    PS::F64 hill_radius = 0;
    PS::F64 specific_energy = 0;
    PS::F64 semi = 0;
    PS::F64 ecc = 0;
    PS::F64 peri = 0;
    PS::F64 apo = 0;
    PS::F64 period = 0;
};

inline PS::F64 calcMutualHillRadius(
    const PS::F64 mass,
    const PS::F64vec& center,
    const PS::F64 central_mass)
{
    if (!(mass > 0) || !(central_mass > 0)|| !std::isfinite(mass) || !std::isfinite(central_mass))
        throw std::domain_error("binary: invalid mass for Hill radius");

    return std::sqrt(center * center) * std::cbrt(mass / (3 * central_mass));
}

// Determine whether two particles should be considered as a binary candidate.
template<class Particle1, class Particle2>
inline bool checkBinaryCandidate(
    const Particle1& a,
    const Particle2& b,
    const BinarySearchParameters& param,
    BinaryInfo& binary)
{
    if (a.isDead || b.isDead || a.id == b.id || !(a.mass > 0) || !(b.mass > 0))
        return false;
    const auto r_rel = a.pos-b.pos;
    const PS::F64 inner = param.inner_cutoff_ratio*std::max(a.getROut(), b.getROut());
    // Reject the changeover shell before computing Hill/orbital quantities.
    if (!(inner > 0) || !std::isfinite(inner)
        || !(r_rel*r_rel < inner*inner)) return false;
    const PS::F64 mass = a.mass + b.mass;
    const auto com = (a.mass / mass) * a.pos + (b.mass / mass) * b.pos;
    const auto x = a.pos - b.pos, v = a.vel - b.vel;
    const PS::F64 r = std::sqrt(x * x);
    if (!(r > 0) || !std::isfinite(r) || !std::isfinite(mass)) return false;
    const PS::F64 hill = calcMutualHillRadius(mass, com, param.central_mass);
    if (!std::isfinite(hill) || !(r / hill < param.hill_factor)) return false;
    if (r <= a.f * a.r_planet + b.f * b.r_planet) return false;
    const PS::F64 gm = param.gravitational_constant * mass;
    const PS::F64 energy = 0.5 * (v * v) - gm / r;
    if (!std::isfinite(energy)) return false;
    if (energy >= 0 && (!param.include_unbound || !(x*v < 0))) return false;
    binary = BinaryInfo{};
    binary.id1 = a.id;
    binary.id2 = b.id;
    binary.separation = r;
    binary.hill_radius = hill;
    binary.specific_energy = energy;
    return true;
}

// Calculate orbital parameters of a binary candidate
template<class Particle1, class Particle2>
inline bool calcBinaryOrbit(
    const Particle1& a,
    const Particle2& b,
    const BinarySearchParameters& param,
    BinaryInfo& binary)
{
    const auto x = a.pos-b.pos;
    const auto v = a.vel-b.vel;
    const PS::F64 r = binary.separation;
    const PS::F64 energy = binary.specific_energy;
    const PS::F64 gm = param.gravitational_constant*(a.mass+b.mass);
    const PS::F64vec angular(x.y*v.z-x.z*v.y, x.z*v.x-x.x*v.z, x.x*v.y-x.y*v.x);
    const PS::F64vec vxh(v.y*angular.z-v.z*angular.y,
        v.z*angular.x-v.x*angular.z, v.x*angular.y-v.y*angular.x);
    const auto evec = vxh / gm - x / r;
    BinaryInfo result = binary;
    result.semi = energy == 0 ? 0 : -gm / (2 * energy);
    result.ecc = std::sqrt(evec * evec);
    result.apo = energy < 0 ? result.semi * (1 + result.ecc) : r;
    result.peri = (angular * angular) / (gm * (1 + result.ecc));
    result.period = energy < 0
        ? 2 * std::acos(-1.) * result.semi * std::sqrt(result.semi / gm) : 0;
    if (!std::isfinite(result.apo) || !std::isfinite(result.period)) return false;
    const PS::F64 inner = param.inner_cutoff_ratio*std::max(a.getROut(), b.getROut());
    if (!(result.apo < inner)) return false;
    binary = result;
    return true;
}

// 1. calculate the orbital parameters
// 2. if the apocenter is outside the inner cutoff, reject the candidate
template<class Particle1, class Particle2>
inline bool isBinaryCandidate(
    const Particle1& a,
    const Particle2& b,
    const BinarySearchParameters& param,
    BinaryInfo& binary)
{
    BinaryInfo result;
    if (!checkBinaryCandidate(a, b, param, result)
        || !calcBinaryOrbit(a, b, param, result)) return false;
    binary = result;
    return true;
}

// search for binary candidates in the given hard cluster
template<class Particles>
inline std::vector<BinaryInfo> detectBinaries(
    const Particles& particles,
    const BinarySearchParameters& param)
{
    std::vector<BinaryInfo> result;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        for (std::size_t j = i + 1; j < particles.size(); ++j) {
            BinaryInfo candidate;
            if (!isBinaryCandidate(particles[i], particles[j], param, candidate)) continue;
            candidate.first = i;
            candidate.second = j;
            result.push_back(candidate);
        }
    }
    return result;
}

