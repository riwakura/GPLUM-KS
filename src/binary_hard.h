#pragma once
#include "binary.h"
#include "ks_integrator.h"
#include "binary_profile.h"
#include <algorithm>

class BinaryPairForceCoefficients {
public:
    PS::F64 inv, alpha, weight, derivative;
    BinaryPairForceCoefficients(
        const PS::F64vec& dx,
        const PS::F64vec& dv,
        PS::F64 radius,
        PS::F64 rout_inv)
    {
        inv = 1/radius;
        alpha = (dx*dv)*inv*inv;
        weight = 1-cutoff_K(radius, rout_inv);
        derivative = cutoff_dKdt(radius, rout_inv, alpha);
    }
    PS::F64vec acceleration(const PS::F64vec& dx, PS::F64 factor) const
    {
        return factor*weight*dx;
    }
    PS::F64vec jerk(const PS::F64vec& dx, const PS::F64vec& dv, PS::F64 factor) const
    {
        return factor*(weight*dv-(3*alpha*weight+derivative)*dx);
    }
};

struct PairForceStatistics {
    PS::U64 calls = 0;
    PS::U64 geometry_pairs = 0;
    PS::U64 force_pairs = 0;
    PS::U64 target_rows = 0;
};

inline PairForceStatistics& pairForceStatistics()
{
    static thread_local PairForceStatistics statistics;
    return statistics;
}

// Calculate force from particles within the same hard cluster and central body
// If the particle is a binary member, except the force from the other binary member
template<class Particles, class Selected>
inline bool calcExternalHardForce(Particles& particles, const std::vector<int>& owner,
    Selected selected, std::vector<PS::F64vec>& acceleration,
    std::vector<PS::F64vec>& jerk, const bool record_statistics = false)
{
    const std::size_t count = particles.size();
    if (owner.size() != count) return false;
    for (std::size_t i = 0; i < count; ++i) {
        const auto& p = particles[i];
        if (!std::isfinite(p.mass) || !std::isfinite(p.getROut()) || !(p.getROut() > 0)
            || !std::isfinite(p.f*p.r_planet)) return false;
        for (int k = 0; k < 3; ++k) {
            if (!std::isfinite(p.pos[k]) || !std::isfinite(p.vel[k])) return false;
        }
        if (!(p.pos*p.pos+FP_t::eps2_sun > 0)) return false;
    }
    auto& statistics = pairForceStatistics();
    if (record_statistics) ++statistics.calls;
    acceleration.assign(count, PS::F64vec(0));
    jerk.assign(count, PS::F64vec(0));
    for (std::size_t i = 0; i < count; ++i) {
        if (selected(i)) {
            if (record_statistics) ++statistics.target_rows;
            calcStarGravity(particles[i]);
            particles[i].acc_d = particles[i].jerk_d = PS::F64vec(0);
            acceleration[i] = particles[i].acc_s;
            jerk[i] = particles[i].jerk_s;
        }
        for (std::size_t k = 0; k < i; ++k) {
            if (record_statistics) ++statistics.geometry_pairs;
            const auto dx = particles[k].pos-particles[i].pos;
            const auto dv = particles[k].vel-particles[i].vel;
            const PS::F64 radius = std::sqrt(dx*dx);
            const PS::F64 contact = particles[i].f*particles[i].r_planet
                +particles[k].f*particles[k].r_planet;
            if (!(radius > contact) || !std::isfinite(radius)) return false;
            if (owner[i] == owner[k]) {
                const PS::F64 inner = FP_t::gamma*std::max(particles[i].getROut(), particles[k].getROut());
                if (!(radius < inner)) return false;
                continue;
            }
            if (!selected(i) && !selected(k)) continue;
            if (record_statistics) ++statistics.force_pairs;
            const BinaryPairForceCoefficients force(dx, dv, radius,
                1/std::max(particles[i].getROut(), particles[k].getROut()));
            const PS::F64 factor = force.inv*force.inv*force.inv;
            if (force.weight == 0 && force.derivative == 0 && std::isfinite(3*force.alpha)
                && std::isfinite(factor) && std::isfinite(dv.x)
                && std::isfinite(dv.y) && std::isfinite(dv.z)
                && std::isfinite(particles[k].mass) && std::isfinite(particles[i].mass)) continue;
            const auto a = force.acceleration(dx, factor);
            const auto j = force.jerk(dx, dv, factor);
            if (selected(i)) {
                const auto da = particles[k].mass*a, dj = particles[k].mass*j;
                acceleration[i] += da;
                jerk[i] += dj;
                particles[i].acc_d += da;
                particles[i].jerk_d += dj;
            }
            if (selected(k)) {
                const auto da = particles[i].mass*a, dj = particles[i].mass*j;
                acceleration[k] -= da;
                jerk[k] -= dj;
                particles[k].acc_d -= da;
                particles[k].jerk_d -= dj;
            }
        }
    }
    for (std::size_t i = 0; i < count; ++i) {
        if (!selected(i)) continue;
        for (int k = 0; k < 3; ++k) {
            if (!std::isfinite(acceleration[i][k]) || !std::isfinite(jerk[i][k])) return false;
        }
    }
    return true;
}

// Calculate acceleration and jerk of binary particles
template<class Particles>
inline bool calcBinaryPerturbation(const Particles& particles,
    const std::vector<int>& owner, std::size_t first, std::size_t second,
    const std::vector<PS::F64vec>& position, const std::vector<PS::F64vec>& velocity,
    typename Particles::value_type& left, typename Particles::value_type& right,
    PS::F64vec& perturbation, PS::F64vec& perturbation_jerk)
{
    PS::F64vec acceleration[2], jerk[2];
    int component = 0;
    for (const auto i : {first, second}) {
        auto& target = component == 0 ? left : right;
        calcStarGravity(target);
        acceleration[component] = target.acc_s;
        jerk[component] = target.jerk_s;
        for (std::size_t k = 0; k < particles.size(); ++k) {
            if (owner[i] == owner[k]) continue;
            const auto dx = position[owner[k]]-target.pos;
            const auto dv = velocity[owner[k]]-target.vel;
            const PS::F64 radius = std::sqrt(dx*dx);
            const PS::F64 contact = particles[i].f*particles[i].r_planet
                +particles[k].f*particles[k].r_planet;
            if (!(radius > contact) || !std::isfinite(radius)) return false;
            const BinaryPairForceCoefficients force(dx, dv, radius,
                1/std::max(particles[i].getROut(), particles[k].getROut()));
            const PS::F64 factor = particles[k].mass*force.inv*force.inv*force.inv;
            // Omit exact zero hard forces, retaining rejection of non-finite results.
            if (force.weight == 0 && force.derivative == 0 && std::isfinite(factor)
                && std::isfinite(3 * force.alpha) && std::isfinite(dv.x)
                && std::isfinite(dv.y) && std::isfinite(dv.z)) continue;
            acceleration[component] += force.acceleration(dx, factor);
            jerk[component] += force.jerk(dx, dv, factor);
        }
        for (int k = 0; k < 3; ++k)
            if (!std::isfinite(acceleration[component][k]) || !std::isfinite(jerk[component][k])) return false;
        ++component;
    }
    perturbation = acceleration[0]-acceleration[1];
    perturbation_jerk = jerk[0]-jerk[1];
    return true;
}

template<class Particles>
inline bool calcAllExternalHardForces(Particles& particles, const std::vector<int>& owner,
    std::vector<PS::F64vec>& acceleration, std::vector<PS::F64vec>& jerk)
{
    return calcExternalHardForce(particles, owner,
        [](std::size_t) { return true; }, acceleration, jerk);
}

/*
 * For Diagnostic
 */

struct BinaryEncounterDiagnostic {
    PS::S32 stage = 0;
    PS::S64 ids[4] = {-1, -1, -1, -1};
    bool binary[2] = {false, false};
    PS::F64 time = 0, distance = 0, buffer = 0, threshold = 0;
    PS::F64 extent[2] = {0, 0};
};

struct BinaryHardReport {
    std::size_t groups = 0, blocks = 0, partial_blocks = 0;
    std::size_t screened_candidates = 0;
    PS::U64 coupling_iterations = 0, full_force_evaluations = 0;
    PS::U64 omitted_tide_samples = 0;
    PS::F64 max_omitted_tide_ratio = 0;
    KSIntegrationStatistics integration;
    BinaryEncounterDiagnostic encounter;
    bool accepted = false;
    int failure_line = 0;
    BinaryFailure reason = BinaryFailure::none;
    PS::F64 resume_time = 0, attempt_seconds = 0;
};

class BinaryHardParameters {
public:
    // Experimental: admit isolated incoming non-bound pairs.
    #ifdef KS_HYPERBOLIC
    bool enable_unbound_encounters = true;
#else
    bool enable_unbound_encounters = false;
#endif
    // Experimental: omit non-stellar differential hard forces inside KS only.
#ifdef KS_SOLAR_ONLY
    bool solar_only = true;
#else
    bool solar_only = false;
#endif
    PS::F64 hill_factor = FP_t::binary_hill_factor;
    PS::F64 max_perturbation = FP_t::binary_max_perturbation;
    // Geometric opening criterion
    PS::F64 max_group_opening = 0.05;
    bool stellar_tide_precheck = true;
    bool cheap_pair_screen = true;
    // Experimental contact margin; collision events still belong to Hermite.
#ifdef KS_CONTACT_FACTOR
    PS::F64 contact_factor = KS_CONTACT_FACTOR;
#else
    PS::F64 contact_factor = 2;
#endif
    PS::F64 coupling_tolerance = 1e-10;
    // Limit only binary centre-of-mass steps; 1 restores the original rule.
    // Block rounding means the resulting step is not always exactly scaled.
    PS::F64 centre_of_mass_step_factor = 1;
    bool require_step_advantage = false;
    PS::S32 max_iterations = 12;
    KSIntegratorParameters internal;
    BinaryHardParameters() {
        // The trajectory error must be below the coupling convergence tolerance.
        internal.iteration_tolerance = 2e-14;
    }
};

// Bound binaries claim their members before non-bound encounters.
inline bool binaryPriorityLess(const BinaryInfo& a, const BinaryInfo& b)
{
    const bool a_bound = a.specific_energy < 0;
    const bool b_bound = b.specific_energy < 0;
    if (a_bound != b_bound) return a_bound;
    if (a.period != b.period) return a.period < b.period;
    if (a.id1 != b.id1) return a.id1 < b.id1;
    return a.id2 < b.id2;
}

// A single particle or a binary centre of mass advanced with Hermite.
struct HermiteParticle {
    std::size_t first = 0, second = 0;
    bool binary = false;
    PS::F64 mass = 0, time = 0, dt = 0;
    PS::F64 first_weight = 0, second_weight = 0;
    PS::F64vec pos, vel, acc, jerk;
    PS::F64vec snap{0}, crackle{0};
    PS::F64 acc0 = 0;
    PS::F64vec acc_d{0}, jerk_d{0}, snap_d{0}, crackle_d{0};
    PS::F64vec acc_s{0}, jerk_s{0}, snap_s{0}, crackle_s{0};
    bool has_high_derivatives = false;
};

// Internal KS state and the index of its Hermite centre of mass.
struct KSBinary {
    std::size_t com_index = 0;
    KSIntegrationState state;
    PS::F64 encounter_radius = 0;
};

inline void predictHermiteParticle(
    const HermiteParticle& hermite_particle, 
    PS::F64 time,
    PS::F64vec& pos, PS::F64vec& vel)
{
    const PS::F64 dt = time-hermite_particle.time;
    pos = hermite_particle.pos + dt*(hermite_particle.vel + dt*(hermite_particle.acc/2 + dt*hermite_particle.jerk/6));
    vel = hermite_particle.vel + dt*(hermite_particle.acc + dt*hermite_particle.jerk/2);
}

inline PS::F64 chooseBinaryHardStep(
    const HermiteParticle& hermite_particle, 
    const PS::F64 origin,
    const PS::F64 centre_of_mass_step_factor = 1)
{
    PS::F64 speed_scale;
    if (hermite_particle.binary) {
        // These particles combine stellar and particle forces. Respect the stricter
        // accuracy setting; use the initial criterion until endpoint derivatives exist.
        const PS::F64 eta = std::min(FP_t::eta, FP_t::eta_sun);
        const PS::F64 eta_initial = std::min(FP_t::eta_0, FP_t::eta_sun0);
        speed_scale = calcDt2nd(hermite_particle.has_high_derivatives ? eta : eta_initial,
            0., 0., hermite_particle.acc, hermite_particle.jerk);
        if (hermite_particle.has_high_derivatives) {
            const PS::F64 fourth_order = calcDt4th(eta, 0., 0., hermite_particle.acc, hermite_particle.jerk,
                hermite_particle.snap, hermite_particle.crackle);
            if (std::isfinite(fourth_order) && fourth_order >= 0)
                speed_scale = std::min(speed_scale, fourth_order);
        }
    } else if (hermite_particle.has_high_derivatives) {
        speed_scale = std::min(
            calcDt4th(FP_t::eta, FP_t::alpha2, hermite_particle.acc0,
                hermite_particle.acc_d, hermite_particle.jerk_d, hermite_particle.snap_d, hermite_particle.crackle_d),
            calcDt4th(FP_t::eta_sun, FP_t::alpha2, 0.,
                hermite_particle.acc_s, hermite_particle.jerk_s, hermite_particle.snap_s, hermite_particle.crackle_s));
    } else {
        speed_scale = std::min(
            calcDt2nd(FP_t::eta_0, FP_t::alpha2, hermite_particle.acc0,
                hermite_particle.acc_d, hermite_particle.jerk_d),
            calcDt2nd(FP_t::eta_sun0, FP_t::alpha2, 0.,
                hermite_particle.acc_s, hermite_particle.jerk_s));
    }
    if (hermite_particle.binary) speed_scale *= centre_of_mass_step_factor;
    PS::F64 dt = 0.5*FP_t::dt_tree;
    if (hermite_particle.dt > 0) dt = std::min(dt, 2*hermite_particle.dt);
    while (dt > FP_t::dt_min && (dt > speed_scale || std::fmod(hermite_particle.time-origin, dt) != 0)) dt *= 0.5;
    return std::max(dt, FP_t::dt_min);
}

template<class Particles>
inline void refreshBinaryHardParticles(Particles& particles, PS::F64 time)
{
    for (std::size_t i = 0; i < particles.size(); ++i) {
        particles[i].time = time;
        particles[i].xp = particles[i].pos;
        particles[i].vp = particles[i].vel;
    }
    for (std::size_t i = 0; i < particles.size(); ++i) {
        calcGravity(particles[i], particles);
#ifdef INTEGRATE_6TH_SUN
        particles[i].setAcc_();
#endif
    }
    for (std::size_t i = 0; i < particles.size(); ++i) {
#ifdef INTEGRATE_6TH_SUN
        // Solar snap depends on the updated total acceleration, not its old value.
        // Reevaluate only the stellar term; do not repeat the particle pair loop.
        calcStarGravity(particles[i]);
#endif
        particles[i].calcDeltatInitial();
    }
}

// Exact lower bound for an unsoftened central point mass:
// |a(x)-a(y)| >= M_star |x-y| / max(|x|,|y|)^3.
// Only usable without other hard perturbers (their vector force may cancel it).
inline bool exceedsMinimumStellarTide(const PS::F64vec& x, const PS::F64vec& y,
    const PS::F64 pair_mass, const PS::F64 star_mass, const PS::F64 limit)
{
    const auto dr = x-y;
    const PS::F64 outer2 = std::max(x*x, y*y);
    if (!(outer2 > 0) || !(pair_mass > 0) || !(star_mass > 0) || !(limit > 0)) return false;
    const PS::F64 q = (dr*dr)/outer2;
    const PS::F64 threshold = limit*(pair_mass/star_mass);
    // Leave a roundoff margin at the boundary; uncertain cases use full forces.
    return q*q*q > threshold*threshold*(1+1e-10);
}

// Only used when there are two particles in a hard cluster
// if the peturbation from sun exceeds the threshold, reject the pair
template<class Particle1, class Particle2>
inline BinaryFailure screenBinaryPairGeometry(const Particle1& a, const Particle2& b,
    const BinarySearchParameters& search, const BinaryHardParameters& parameters,
    BinaryInfo& candidate)
{
    if (!checkBinaryCandidate(a, b, search, candidate))
        return BinaryFailure::no_candidate;
    if (parameters.stellar_tide_precheck && FP_t::eps2_sun == 0
        && exceedsMinimumStellarTide(a.pos, b.pos, a.mass+b.mass,
            FP_t::m_sun, parameters.max_perturbation))
        return BinaryFailure::perturbation;
    const PS::F64 contact = a.f*a.r_planet+b.f*b.r_planet;
    const PS::F64 full_hard = FP_t::gamma*std::max(a.getROut(), b.getROut());
    if (candidate.separation < parameters.contact_factor*contact*(1-1e-10)
        || candidate.separation > full_hard*(1+1e-10))
        return BinaryFailure::unsafe_candidate;
    return BinaryFailure::none;
}

// Check whether short Hermite steps make a KS attemps worthwhile
#ifndef BINARY_MIN_HERMITE_STEPS
#define BINARY_MIN_HERMITE_STEPS 32
#endif
template<class Particles>
inline bool binaryHardMaySaveWork(Particles& input, const PS::F64 remaining)
{
    if (BINARY_MIN_HERMITE_STEPS <= 0) return true;
    if (!(remaining > 0)) return false;
    const PS::F64 limit = remaining/BINARY_MIN_HERMITE_STEPS;
    for (std::size_t i = 0; i < input.size(); ++i) {
        const auto& p = input[i];
        if (p.isDead || p.isMerged || !(p.dt > 0) || p.dt > limit) continue;
#ifndef WITHOUT_SUN
        if (FP_t::m_sun > 0) {
            const PS::F64 r2 = p.pos*p.pos;
            const PS::F64 softened = r2+FP_t::eps2_sun;
            if (!(r2 > 0) || !(softened > 0)) continue;
            const PS::F64vec rate = p.vel-(3*(p.pos*p.vel)/softened)*p.pos;
            const PS::F64 eta2 = FP_t::eta_sun0*FP_t::eta_sun0;
            if (!(4*p.dt*p.dt*(rate*rate) < eta2*r2)) continue;
        }
#endif
        return true;
    }
    return false;
}

// Experimental: Set the KS encounter radius for an unbound pair (hyperbolic encounter)
template<class Particle1, class Particle2>
inline PS::F64 binaryPairEncounterRadius(
    const Particle1& a,
    const Particle2& b,
    const BinaryInfo& pair,
    const BinaryHardParameters& parameters)
{
    return pair.specific_energy >= 0
        ? std::min(2*pair.separation, std::min(pair.hill_radius*parameters.hill_factor,
            FP_t::gamma*std::max(a.getROut(), b.getROut()))) : 0;
}

// Check the binary pair stays inside the hard radius and do not collide
template<class Particle1, class Particle2>
inline bool binaryPairOrbitAdmissible(const Particle1& a, const Particle2& b,
    const BinaryInfo& pair, const BinaryHardParameters& parameters)
{
    const PS::F64 radius = binaryPairEncounterRadius(a, b, pair, parameters);
    const PS::F64 contact = a.f*a.r_planet+b.f*b.r_planet;
    return !(radius > 0 && !(pair.separation < radius))
        && pair.apo < FP_t::gamma*std::max(a.getROut(), b.getROut())
        && pair.peri > parameters.contact_factor*contact;
}

// Check whether current state permits another KS attempt.
// Only used for two-particle hard cluster
template<class Particles>
inline bool binaryHardRescanMayAdmit(Particles& input,
    const BinaryHardParameters& parameters = {})
{
    if (!(parameters.hill_factor > 0) || !std::isfinite(parameters.hill_factor)) return false;
    if (!parameters.cheap_pair_screen || input.size() != 2) return true;
    BinarySearchParameters search;
    search.include_unbound = parameters.enable_unbound_encounters;
    search.hill_factor = parameters.hill_factor;
    BinaryInfo candidate;
    if (screenBinaryPairGeometry(input[0], input[1], search, parameters, candidate)
        != BinaryFailure::none) return false;
    if (!calcBinaryOrbit(input[0], input[1], search, candidate)) return false;
    return binaryPairOrbitAdmissible(input[0], input[1], candidate, parameters);
}

// integrate hard cluster using KS inner motion and whole Hermite scheme including com particle of binary
template<class Particles>
inline bool integrateBinaryHardAttempt(Particles& input, const PS::F64 begin,
    const PS::F64 end, const BinaryHardParameters& parameters,
    BinaryHardReport* report, const bool keep_prefix)
{
    *report = BinaryHardReport{};
    report->resume_time = begin;
    std::vector<FPHard> checkpoint;
    PS::F64 checkpoint_time = begin;
    KSIntegratorParameters internal_parameters = parameters.internal;
    internal_parameters.statistics = FP_t::ks_profile ? &report->integration : nullptr;
    const auto fail = [&](int line, BinaryFailure reason) {
        if (report->failure_line == 0) {
            report->failure_line = line;
            report->reason = reason;
        }
        if (keep_prefix && checkpoint_time > begin) {
            for (std::size_t i = 0; i < input.size(); ++i) input[i] = checkpoint[i];
            refreshBinaryHardParticles(input, checkpoint_time);
            report->resume_time = checkpoint_time;
        }
        return false;
    };
#if defined(WITHOUT_SUN) || defined(MERGE_BINARY) || defined(KOMINAMI) || defined(CHAMBERS) || defined(SHIBATA) || defined(TEST_PTCL)
    return fail(__LINE__, BinaryFailure::unsupported);
#else
    if (!(end > begin) || FP_t::eps2 != 0 || !(FP_t::m_sun > 0)
        || !(FP_t::dt_min > 0) || !(parameters.max_perturbation > 0)
        || !(parameters.max_group_opening > 0 && parameters.max_group_opening < 1)
        || !(parameters.contact_factor >= 1) || !std::isfinite(parameters.contact_factor)
        || !(parameters.centre_of_mass_step_factor > 0)
        || !std::isfinite(parameters.centre_of_mass_step_factor)
        || !(parameters.hill_factor > 0) || !std::isfinite(parameters.hill_factor)) return fail(__LINE__, BinaryFailure::invalid_input);
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (input[i].isDead || input[i].isMerged || input[i].time != begin) return fail(__LINE__, BinaryFailure::invalid_input);
    }
    BinarySearchParameters search;
    search.include_unbound = parameters.enable_unbound_encounters;
    search.hill_factor = parameters.hill_factor;
    std::vector<BinaryInfo> pairs;
    if (parameters.cheap_pair_screen && input.size() == 2) {
        BinaryInfo candidate;
        const auto reason = screenBinaryPairGeometry(input[0], input[1], search,
            parameters, candidate);
        if (reason != BinaryFailure::none) return fail(__LINE__, reason);
        if (!calcBinaryOrbit(input[0], input[1], search, candidate))
            return fail(__LINE__, BinaryFailure::no_candidate);
        candidate.first = 0;
        candidate.second = 1;
        pairs.push_back(candidate);
    } else {
        pairs = detectBinaries(input, search);
        std::sort(pairs.begin(), pairs.end(), binaryPriorityLess);
    }
    if (!parameters.cheap_pair_screen && parameters.stellar_tide_precheck
        && input.size() == 2 && !pairs.empty()
        && FP_t::eps2_sun == 0
        && exceedsMinimumStellarTide(input[0].pos, input[1].pos,
            input[0].mass+input[1].mass, FP_t::m_sun, parameters.max_perturbation))
        return fail(__LINE__, BinaryFailure::perturbation);
    std::vector<HermiteParticle> hermite_particles;
    std::vector<KSBinary> groups;
    std::vector<bool> assigned(input.size(), false);
    for (const auto& pair : pairs) {
        if (assigned[pair.first] || assigned[pair.second]) continue;
        const auto& a = input[pair.first];
        const auto& b = input[pair.second];
        if (!binaryPairOrbitAdmissible(a, b, pair, parameters)) continue;
        const PS::F64 encounter_radius = binaryPairEncounterRadius(a, b, pair, parameters);
        HermiteParticle hermite_particle;
        hermite_particle.first = pair.first;
        hermite_particle.second = pair.second;
        hermite_particle.binary = true;
        hermite_particle.mass = a.mass+b.mass;
        hermite_particle.first_weight = a.mass/hermite_particle.mass;
        hermite_particle.second_weight = b.mass/hermite_particle.mass;
        hermite_particle.pos = hermite_particle.first_weight*a.pos+hermite_particle.second_weight*b.pos;
        hermite_particle.vel = hermite_particle.first_weight*a.vel+hermite_particle.second_weight*b.vel;
        hermite_particle.time = begin;
        bool isolated = true;
        const PS::F64 exclusion = std::max(pair.apo, encounter_radius)/parameters.max_group_opening;
        for (std::size_t k = 0; k < input.size(); ++k) {
            if (k == pair.first || k == pair.second) continue;
            const auto dx = input[k].pos-hermite_particle.pos;
            if (!(dx*dx > exclusion*exclusion)) { isolated = false; break; }
        }
        if (!isolated) {
            if (FP_t::ks_profile) ++report->screened_candidates;
            continue;
        }
        KSBinary group;
        group.com_index = hermite_particles.size();
        group.encounter_radius = encounter_radius;
        group.state = initializeKSBinary(a.pos-b.pos, a.vel-b.vel, hermite_particle.mass, begin);
        hermite_particles.push_back(hermite_particle);
        groups.push_back(group);
        assigned[pair.first] = assigned[pair.second] = true;
    }
    if (groups.empty()) return fail(__LINE__, pairs.empty()
        ? BinaryFailure::no_candidate : BinaryFailure::unsafe_candidate);
    if (parameters.require_step_advantage) {
        PS::F64 paired_step = std::numeric_limits<PS::F64>::infinity();
        PS::F64 single_step = std::numeric_limits<PS::F64>::infinity();
        for (std::size_t i = 0; i < input.size(); ++i) {
            if (!(input[i].dt > 0) || !std::isfinite(input[i].dt)) continue;
            auto& step = assigned[i] ? paired_step : single_step;
            step = std::min(step, input[i].dt);
        }
        // Require two larger block levels to offset KS synchronization work.
        if (std::isfinite(single_step) && single_step < 4*paired_step)
            return fail(__LINE__, BinaryFailure::cost);
    }
    std::vector<FPHard> original;
    original.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) original.push_back(input[i]);
    std::vector<int> owner(input.size(), -1);
    for (std::size_t n = 0; n < hermite_particles.size(); ++n)
        owner[hermite_particles[n].first] = owner[hermite_particles[n].second] = n;
    if (report) report->groups = groups.size();
    for (std::size_t i = 0; i < original.size(); ++i) {
        if (owner[i] >= 0) continue;
        HermiteParticle hermite_particle;
        hermite_particle.first = hermite_particle.second = i;
        hermite_particle.mass = original[i].mass;
        hermite_particle.dt = original[i].dt;
        hermite_particle.acc0 = original[i].acc0;
        hermite_particle.pos = original[i].pos;
        hermite_particle.vel = original[i].vel;
        hermite_particle.time = begin;
        owner[i] = hermite_particles.size();
        hermite_particles.push_back(hermite_particle);
    }
    std::vector<PS::F64vec> a, j;
    auto physical = original;
    if (FP_t::ks_profile) ++report->full_force_evaluations;
    if (!calcAllExternalHardForces(physical, owner, a, j)) return fail(__LINE__, BinaryFailure::force);
    const auto particleForce = [&](const HermiteParticle& hermite_particle, const std::vector<PS::F64vec>& force) {
        if (!hermite_particle.binary) return force[hermite_particle.first];
        return hermite_particle.first_weight*force[hermite_particle.first]
            + hermite_particle.second_weight*force[hermite_particle.second];
    };
    const auto setSingleForces = [&](HermiteParticle& particle) {
        if (particle.binary) return;
        const auto& physical_particle = physical[particle.first];
        particle.acc_d = physical_particle.acc_d;
        particle.jerk_d = physical_particle.jerk_d;
        particle.acc_s = physical_particle.acc_s;
        particle.jerk_s = physical_particle.jerk_s;
    };
    for (auto& hermite_particle : hermite_particles) {
        setSingleForces(hermite_particle);
        hermite_particle.acc = particleForce(hermite_particle, a);
        hermite_particle.jerk = particleForce(hermite_particle, j);
        hermite_particle.dt = chooseBinaryHardStep(hermite_particle, begin, parameters.centre_of_mass_step_factor);
    }
    std::vector<PS::F64vec> source_position(hermite_particles.size()), source_velocity(hermite_particles.size());
    struct StepCoefficients {
        PS::F64 dt = 0, inverse = 0, half = 0, square_over_twelve = 0;
    };
    std::vector<StepCoefficients> step_coefficients(hermite_particles.size());
    PS::F64 current = begin;
    for (int block = 0; current < end && block < 1000000; ++block) {
        if (FP_t::ks_profile) ++report->blocks;
        PS::F64 target = end;
        for (const auto& hermite_particle : hermite_particles) target = std::min(target, hermite_particle.time+hermite_particle.dt);
        if (!(target > current)) return fail(__LINE__, BinaryFailure::schedule);
        if (FP_t::ks_profile) {
            for (const auto& hermite_particle : hermite_particles)
                if (hermite_particle.time+hermite_particle.dt > target) { ++report->partial_blocks; break; }
        }
        const PS::F64 duration = target-current;
        std::vector<PS::F64> extent(hermite_particles.size(), 0);
        for (const auto& group : groups) {
            const auto& hermite_particle = hermite_particles[group.com_index];
            const auto x = ksToPos(group.state.u.data());
            const auto v = ksToVel(group.state.u.data(), group.state.uprime.data());
            const PS::F64 r = std::sqrt(x*x), energy = 0.5*(v*v)-hermite_particle.mass/r;
            if (group.encounter_radius > 0) {
                extent[group.com_index] = group.encounter_radius;
                continue;
            }
            if (!(energy < 0)) return fail(__LINE__, BinaryFailure::unbound);
            const auto angular = calcCrossKSVector(x, v);
            const PS::F64 eccentricity = std::sqrt(std::max(0.,
                1+2*energy*(angular*angular)/(hermite_particle.mass*hermite_particle.mass)));
            extent[group.com_index] = -hermite_particle.mass/(2*energy)*(1+eccentricity);
        }
        if (hermite_particles.size() > 1) {
            for (std::size_t n = 0; n < hermite_particles.size(); ++n)
                predictHermiteParticle(hermite_particles[n], current, source_position[n], source_velocity[n]);
        }
        for (std::size_t n = 0; n < hermite_particles.size(); ++n) {
            for (std::size_t m = 0; m < n; ++m) {
                if (!hermite_particles[n].binary && !hermite_particles[m].binary) continue;
                const auto dx = source_position[n]-source_position[m];
                const auto dv = source_velocity[n]-source_velocity[m];
                const auto da = hermite_particles[n].acc-hermite_particles[m].acc;
                const PS::F64 buffer = 2*duration*std::sqrt(dv*dv)
                    +4*duration*duration*std::sqrt(da*da);
                if (!(std::sqrt(dx*dx)-buffer >
                    (extent[n]+extent[m])/parameters.max_group_opening)) {
                    auto& diagnostic = report->encounter;
                    diagnostic.stage = 1;
                    diagnostic.time = current;
                    diagnostic.ids[0] = original[hermite_particles[n].first].id;
                    diagnostic.ids[1] = original[hermite_particles[n].second].id;
                    diagnostic.ids[2] = original[hermite_particles[m].first].id;
                    diagnostic.ids[3] = original[hermite_particles[m].second].id;
                    diagnostic.binary[0] = hermite_particles[n].binary;
                    diagnostic.binary[1] = hermite_particles[m].binary;
                    diagnostic.distance = std::sqrt(dx*dx);
                    diagnostic.buffer = buffer;
                    diagnostic.threshold = (extent[n]+extent[m])/parameters.max_group_opening;
                    diagnostic.extent[0] = extent[n];
                    diagnostic.extent[1] = extent[m];
                    return fail(__LINE__, BinaryFailure::encounter);
                }
            }
        }
        for (std::size_t n = 0; n < hermite_particles.size(); ++n) {
            if (target != end && hermite_particles[n].time+hermite_particles[n].dt > target) continue;
            auto& step = step_coefficients[n];
            step.dt = target-hermite_particles[n].time;
            step.inverse = 1/step.dt;
            step.half = 0.5*step.dt;
            step.square_over_twelve = (step.dt*step.dt)*(1.0/12.0);
        }
        auto trial_particles = hermite_particles;
        for (auto& hermite_particle : trial_particles) predictHermiteParticle(hermite_particle, target, hermite_particle.pos, hermite_particle.vel);
        auto trial_groups = groups;
        bool converged = false;
        const auto predict_sources = [&](PS::F64 time, std::size_t target_node) {
            for (std::size_t n = 0; n < hermite_particles.size(); ++n) {
                if (parameters.solar_only && n != target_node) continue;
                PS::F64vec pos, vel;
                if (target == end || hermite_particles[n].time+hermite_particles[n].dt <= target) {
                    const auto& step = step_coefficients[n];
                    const PS::F64 dt = step.dt;
                    const PS::F64 x = (time-hermite_particles[n].time)*step.inverse;
                    const auto& left = hermite_particles[n];
                    const auto& right = trial_particles[n];
                    pos = (2*x*x*x-3*x*x+1)*left.pos+(x*x*x-2*x*x+x)*dt*left.vel
                        +(-2*x*x*x+3*x*x)*right.pos+(x*x*x-x*x)*dt*right.vel;
                    vel = (6*x*x-6*x)*step.inverse*left.pos+(3*x*x-4*x+1)*left.vel
                        +(-6*x*x+6*x)*step.inverse*right.pos+(3*x*x-2*x)*right.vel;
                } else {
                    predictHermiteParticle(hermite_particles[n], time, pos, vel);
                }
                source_position[n] = pos;
                source_velocity[n] = vel;
            }
        };
        const auto materialize_endpoint = [&]() {
            physical = original;
            for (std::size_t n = 0; n < hermite_particles.size(); ++n) {
                physical[hermite_particles[n].first].pos = trial_particles[n].pos;
                physical[hermite_particles[n].first].vel = trial_particles[n].vel;
                if (hermite_particles[n].binary) {
                    physical[hermite_particles[n].second].pos = trial_particles[n].pos;
                    physical[hermite_particles[n].second].vel = trial_particles[n].vel;
                }
            }
            for (std::size_t g = 0; g < groups.size(); ++g) {
                const auto& hermite_particle = hermite_particles[groups[g].com_index];
                const auto& state = trial_groups[g].state;
                const auto x = ksToPos(state.u.data()), v = ksToVel(state.u.data(), state.uprime.data());
                const PS::F64 fraction = hermite_particle.second_weight;
                physical[hermite_particle.first].pos += fraction*x;
                physical[hermite_particle.first].vel += fraction*v;
                physical[hermite_particle.second].pos -= (1-fraction)*x;
                physical[hermite_particle.second].vel -= (1-fraction)*v;
            }
        };
        for (int iteration = 0; iteration < parameters.max_iterations; ++iteration) {
            if (FP_t::ks_profile) ++report->coupling_iterations;
            const auto previous_groups = trial_groups;
            const auto previous_nodes = trial_particles;
            for (std::size_t g = 0; g < groups.size(); ++g) {
                const auto& hermite_particle = hermite_particles[groups[g].com_index];
                auto left = original[hermite_particle.first], right = original[hermite_particle.second];
                BinaryFailure force_reason = BinaryFailure::none;
                int force_failure_line = 0;
                const auto reject_force = [&](int line, BinaryFailure reason, bool retryable) {
                    force_reason = reason;
                    force_failure_line = line;
                    return retryable ? KSForceResult::retry : KSForceResult::failure;
                };
                auto force = [&](const KSIntegrationState& state, PS::F64vec& p, PS::F64vec& jp) {
                    force_reason = BinaryFailure::none;
                    force_failure_line = 0;
                    predict_sources(state.time, groups[g].com_index);
                    const auto& hermite_particle = hermite_particles[groups[g].com_index];
                    const auto x = ksToPos(state.u.data());
                    const auto v = ksToVel(state.u.data(), state.uprime.data());
                    const PS::F64 fraction = hermite_particle.second_weight;
                    left.pos = right.pos = source_position[groups[g].com_index];
                    left.vel = right.vel = source_velocity[groups[g].com_index];
                    left.pos += fraction*x;
                    left.vel += fraction*v;
                    right.pos -= (1-fraction)*x;
                    right.vel -= (1-fraction)*v;
                    if (parameters.solar_only) {
                        calcStarGravity(left);
                        calcStarGravity(right);
                        p = left.acc_s-right.acc_s;
                        jp = left.jerk_s-right.jerk_s;
                    } else if (!calcBinaryPerturbation(original, owner, hermite_particle.first, hermite_particle.second,
                        source_position, source_velocity, left, right, p, jp)) {
                        return reject_force(__LINE__, BinaryFailure::force, false);
                    }
                    const PS::F64 r = std::sqrt(x*x);
                    if (std::sqrt(p*p)*r*r/hermite_particle.mass > parameters.max_perturbation) return reject_force(__LINE__, BinaryFailure::perturbation, state.time > target);
                    const PS::F64 energy = 0.5*(v*v)-hermite_particle.mass/r;
                    if (groups[g].encounter_radius > 0) {
                        if (!(r < groups[g].encounter_radius))
                            return reject_force(__LINE__, BinaryFailure::encounter, state.time > target);
                    } else if (!(energy < 0)) return reject_force(__LINE__, BinaryFailure::unbound, state.time > target);
                    const auto angular = calcCrossKSVector(x, v);
                    const PS::F64 ecc = std::sqrt(std::max(0., 1+2*energy*(angular*angular)/(hermite_particle.mass*hermite_particle.mass)));
                    const PS::F64 peri = (angular*angular)/(hermite_particle.mass*(1+ecc));
                    const PS::F64 contact = original[hermite_particle.first].f*original[hermite_particle.first].r_planet
                        + original[hermite_particle.second].f*original[hermite_particle.second].r_planet;
                    if (!(peri > parameters.contact_factor*contact)) {
                        return reject_force(__LINE__, BinaryFailure::contact, false);
                    }
                    return KSForceResult::success;
                };
                if (!integrateKSBinaryToTime(groups[g].state, target, force,
                    trial_groups[g].state, internal_parameters)) {
                    return fail(force_failure_line ? force_failure_line : __LINE__,
                        force_reason == BinaryFailure::none ? BinaryFailure::integrator : force_reason);
                }
            }
            materialize_endpoint();
            if (FP_t::ks_profile) ++report->full_force_evaluations;
            if (!calcAllExternalHardForces(physical, owner, a, j)) return fail(__LINE__, BinaryFailure::force);
            for (std::size_t n = 0; n < hermite_particles.size(); ++n) {
                if (target != end && hermite_particles[n].time+hermite_particles[n].dt > target) continue;
                const auto& old = hermite_particles[n];
                auto& next = trial_particles[n];
                const auto& step = step_coefficients[n];
                next.acc = particleForce(old, a);
                next.jerk = particleForce(old, j);
                next.vel = old.vel+step.half*(old.acc+next.acc)-step.square_over_twelve*(next.jerk-old.jerk);
                next.pos = old.pos+step.half*(old.vel+next.vel)-step.square_over_twelve*(next.acc-old.acc);
            }
            PS::F64 error = 0;
            for (std::size_t g = 0; g < groups.size(); ++g)
                error = std::max(error, differenceKSState(trial_groups[g].state, previous_groups[g].state));
            for (std::size_t n = 0; n < hermite_particles.size(); ++n)
                for (int k = 0; k < 3; ++k) {
                    error = std::max(error, std::abs(trial_particles[n].pos[k]-previous_nodes[n].pos[k])
                        /(1+std::abs(trial_particles[n].pos[k])));
                    error = std::max(error, std::abs(trial_particles[n].vel[k]-previous_nodes[n].vel[k])
                        /(1+std::abs(trial_particles[n].vel[k])));
                }
            if (error < parameters.coupling_tolerance) { converged = true; break; }
        }
        if (!converged) return fail(__LINE__, BinaryFailure::coupling);
        materialize_endpoint();
        if (FP_t::ks_profile) ++report->full_force_evaluations;
        if (!calcAllExternalHardForces(physical, owner, a, j)) return fail(__LINE__, BinaryFailure::force);
        for (std::size_t i = 0; i < physical.size(); ++i) {
            for (std::size_t k = 0; k < i; ++k) {
                if (owner[i] == owner[k]) continue;
                const auto dx = physical[i].pos-physical[k].pos;
                const auto dv = physical[i].vel-physical[k].vel;
                const auto da = a[i]-a[k];
                const PS::F64 guard = 2*duration*std::sqrt(dv*dv)+4*duration*duration*std::sqrt(da*da)
                    + physical[i].f*physical[i].r_planet+physical[k].f*physical[k].r_planet;
                if (std::sqrt(dx*dx) < guard) {
                    auto& diagnostic = report->encounter;
                    diagnostic.stage = 2;
                    diagnostic.time = target;
                    diagnostic.ids[0] = diagnostic.ids[1] = physical[i].id;
                    diagnostic.ids[2] = diagnostic.ids[3] = physical[k].id;
                    diagnostic.binary[0] = hermite_particles[owner[i]].binary;
                    diagnostic.binary[1] = hermite_particles[owner[k]].binary;
                    diagnostic.distance = std::sqrt(dx*dx);
                    diagnostic.threshold = guard;
                    return fail(__LINE__, BinaryFailure::encounter);
                }
            }
        }
        for (std::size_t g = 0; g < groups.size(); ++g) {
            const auto& hermite_particle = hermite_particles[groups[g].com_index];
            const auto& state = trial_groups[g].state;
            const auto x = ksToPos(state.u.data()), v = ksToVel(state.u.data(), state.uprime.data());
            const PS::F64 r = std::sqrt(x*x), energy = 0.5*(v*v)-hermite_particle.mass/r;
            if (FP_t::ks_profile && parameters.solar_only) {
                auto left = physical[hermite_particle.first], right = physical[hermite_particle.second];
                calcStarGravity(left);
                calcStarGravity(right);
                const auto omitted = a[hermite_particle.first]-a[hermite_particle.second]-(left.acc_s-right.acc_s);
                const PS::F64 ratio = std::sqrt(omitted*omitted)*r*r/hermite_particle.mass;
                if (!std::isfinite(ratio)) return fail(__LINE__, BinaryFailure::force);
                ++report->omitted_tide_samples;
                report->max_omitted_tide_ratio = std::max(report->max_omitted_tide_ratio, ratio);
            }
            if (groups[g].encounter_radius > 0) {
                if (!(r < groups[g].encounter_radius)) return fail(__LINE__, BinaryFailure::encounter);
            } else if (!(energy < 0)) return fail(__LINE__, BinaryFailure::unbound);
            const auto angular = calcCrossKSVector(x, v);
            const PS::F64 e2 = std::max(0., 1+2*energy*(angular*angular)/(hermite_particle.mass*hermite_particle.mass));
            const PS::F64 peri = (angular*angular)/(hermite_particle.mass*(1+std::sqrt(e2)));
            const PS::F64 contact = physical[hermite_particle.first].f*physical[hermite_particle.first].r_planet
                + physical[hermite_particle.second].f*physical[hermite_particle.second].r_planet;
            if (!(peri > parameters.contact_factor*contact)) return fail(__LINE__, BinaryFailure::contact);
        }
        groups = trial_groups;
        for (std::size_t n = 0; n < hermite_particles.size(); ++n) {
            if (target != end && hermite_particles[n].time+hermite_particles[n].dt > target) continue;
            auto& next = trial_particles[n];
            const auto& previous = hermite_particles[n];
            const auto& step = step_coefficients[n];
            const PS::F64 dt = step.dt;
            const PS::F64 inverse2 = step.inverse*step.inverse;
            const PS::F64 inverse3 = inverse2*step.inverse;
            const auto delta_acc = next.acc-previous.acc;
            next.snap = (-6*delta_acc+dt*(2*previous.jerk+4*next.jerk))*inverse2;
            next.crackle = (-12*delta_acc+6*dt*(previous.jerk+next.jerk))*inverse3;
            if (!next.binary) {
                setSingleForces(next);
                const auto delta_d = next.acc_d-previous.acc_d;
                const auto delta_s = next.acc_s-previous.acc_s;
                next.snap_d = (-6*delta_d+dt*(2*previous.jerk_d+4*next.jerk_d))*inverse2;
                next.crackle_d = (-12*delta_d+6*dt*(previous.jerk_d+next.jerk_d))*inverse3;
                next.snap_s = (-6*delta_s+dt*(2*previous.jerk_s+4*next.jerk_s))*inverse2;
                next.crackle_s = (-12*delta_s+6*dt*(previous.jerk_s+next.jerk_s))*inverse3;
            }
            next.has_high_derivatives = true;
            hermite_particles[n] = next;
            hermite_particles[n].time = target;
            hermite_particles[n].dt = chooseBinaryHardStep(hermite_particles[n], begin, parameters.centre_of_mass_step_factor);
        }
        current = target;
        bool synchronized = true;
        for (const auto& hermite_particle : hermite_particles) synchronized = synchronized && hermite_particle.time == current;
        if (keep_prefix && synchronized) {
            checkpoint = physical;
            checkpoint_time = current;
        }
    }
    if (current != end) return fail(__LINE__, BinaryFailure::step_limit);
    refreshBinaryHardParticles(physical, end);
    for (std::size_t i = 0; i < physical.size(); ++i) input[i] = physical[i];
    if (report) { report->accepted = true; report->failure_line = 0;
        report->reason = BinaryFailure::none; report->resume_time = end; }
    return true;
#endif
}

// attempt to integrate hard cluster and record the result
template<class Particles>
inline bool tryIntegrateBinaryHard(Particles& input, const PS::F64 begin,
    const PS::F64 end, const BinaryHardParameters& parameters = {},
    BinaryHardReport* output_report = nullptr, const bool keep_prefix = false)
{
    BinaryHardReport report;
    BinaryAttemptTimer timer;
    const bool success = integrateBinaryHardAttempt(input, begin, end, parameters, &report, keep_prefix);
    report.attempt_seconds = timer.elapsed();
    if (output_report) *output_report = report;
    recordBinaryHardProfile(report, begin, end, success);
    return success;
}

// Decide when to try KS integration
class BinaryWorkAdmission {
public:
    PS::S64 next_probe = 32;
    PS::U64 attempts = 0, skipped = 0, accepted = 0;
    bool due(const PS::S64 blocks, const bool synchronized) const {
        return blocks >= next_probe && synchronized;
    }
    void defer(const PS::S64 blocks) {
        next_probe = blocks > std::numeric_limits<PS::S64>::max()/4
            ? std::numeric_limits<PS::S64>::max() : 4*blocks;
    }
    ~BinaryWorkAdmission() {
        recordBinaryRescanProfile(attempts, skipped, accepted);
    }
};


// Try switching KS integration from Hermite scheme during hard integration.
template<class Particles>
inline bool tryBinaryHardRescan(Particles& pp, const bool allow_binary,
    const PS::S32 collisions, const PS::S64 blocks, const bool synchronized,
    PS::F64& time, const PS::F64 end, BinaryWorkAdmission& admission,
    BinaryFallbackTimer& fallback)
{
    if (!allow_binary || collisions != 0 || !admission.due(blocks, synchronized)) return false;
    admission.defer(blocks);
    if (!binaryHardMaySaveWork(pp, end-time) || !binaryHardRescanMayAdmit(pp)) {
        if (FP_t::ks_profile) ++admission.skipped;
        return false;
    }
    if (FP_t::ks_profile) ++admission.attempts;
    BinaryHardReport report;
    BinaryHardParameters parameters;
    parameters.require_step_advantage = true;
    const bool resumed = tryIntegrateBinaryHard(pp, time, end, parameters, &report, true);
    if (fallback.enabled) fallback.excluded_seconds += report.attempt_seconds;
    if (resumed) {
        if (FP_t::ks_profile) ++admission.accepted;
        return true;
    }
    if (report.resume_time > time) time = report.resume_time;
    fallback.activate();
    return false;
}
