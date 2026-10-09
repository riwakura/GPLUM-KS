#pragma once
#include "binary.h"
#include "ks_integrator.h"
#include <chrono>
#include <mutex>
#include <ostream>
#ifdef BINARY_ENCOUNTER_TRACE
#include <iostream>
#endif

enum class BinaryFailure {
    none, unsupported, invalid_input, no_candidate, unsafe_candidate,
    force, unbound, encounter, perturbation, contact, integrator,
    coupling, schedule, step_limit, cost, count
};

inline const char* binaryFailureName(BinaryFailure reason)
{
    static const char* names[] = {"none", "unsupported", "invalid_input", "no_candidate",
        "unsafe_candidate", "force", "unbound", "encounter", "perturbation", "contact",
        "integrator", "coupling", "schedule", "step_limit", "cost"};
    return names[static_cast<int>(reason)];
}

struct BinaryProfile {
    PS::U64 attempts = 0;
    PS::U64 accepted = 0;
    PS::U64 groups = 0;
    PS::U64 prefixes = 0;
    PS::U64 screened_candidates = 0;
    PS::U64 blocks = 0;
    PS::U64 coupling_iterations = 0;
    PS::U64 full_force_evaluations = 0;
    PS::U64 fallback_intervals = 0;
    PS::U64 rescan_attempts = 0;
    PS::U64 reactivations = 0;
    PS::U64 rescan_precheck_skips = 0;
    PS::U64 reasons[static_cast<int>(BinaryFailure::count)] = {};
    PS::U64 omitted_tide_samples = 0;
    PS::F64 max_omitted_tide_ratio = 0;
    KSIntegrationStatistics integration;
    PS::F64 attempt_seconds = 0, fallback_seconds = 0;
};

inline BinaryProfile& binaryProfileStorage()
{
    static BinaryProfile profile;
    return profile;
}

inline std::mutex& binaryProfileMutex()
{
    static std::mutex mutex;
    return mutex;
}

inline BinaryProfile getBinaryProfile()
{
    std::lock_guard<std::mutex> lock(binaryProfileMutex());
    return binaryProfileStorage();
}

template<class Report>
inline void recordBinaryHardProfile(
    const Report& report,
    const PS::F64 begin,
    const PS::F64 end,
    const bool success)
{
#ifndef BINARY_ENCOUNTER_TRACE
    if (!FP_t::ks_profile) return;
#endif
    std::lock_guard<std::mutex> lock(binaryProfileMutex());
    auto& p = binaryProfileStorage();
#ifdef BINARY_ENCOUNTER_TRACE
    if (report.reason == BinaryFailure::encounter) {
        const auto& d = report.encounter;
        const auto precision = std::cerr.precision();
        std::cerr.precision(17);
        std::cerr << "BINARY_ENCOUNTER," << begin << ',' << end << ',' << d.time
            << ',' << d.stage << ',' << report.groups;
        for (const auto id : d.ids) std::cerr << ',' << id;
        std::cerr << ',' << d.binary[0] << ',' << d.binary[1] << ',' << d.distance
            << ',' << d.buffer << ',' << d.threshold << ',' << d.extent[0]
            << ',' << d.extent[1] << '\n';
        std::cerr.precision(precision);
    }
#endif
    if (!FP_t::ks_profile) return;
    p.omitted_tide_samples += report.omitted_tide_samples;
    p.max_omitted_tide_ratio = std::max(p.max_omitted_tide_ratio, report.max_omitted_tide_ratio);
    p.screened_candidates += report.screened_candidates;
    ++p.attempts;
    if (success) { ++p.accepted; p.groups += report.groups; }
    else ++p.reasons[static_cast<int>(report.reason)];
    if (!success && report.resume_time > begin) ++p.prefixes;
    p.blocks += report.blocks;
    p.coupling_iterations += report.coupling_iterations;
    p.full_force_evaluations += report.full_force_evaluations;
    p.attempt_seconds += report.attempt_seconds;
    p.integration.step_calls += report.integration.step_calls;
    p.integration.correction_iterations += report.integration.correction_iterations;
    p.integration.force_evaluations += report.integration.force_evaluations;
    p.integration.accepted_steps += report.integration.accepted_steps;
    p.integration.rejected_steps += report.integration.rejected_steps;
    p.integration.reject_iteration += report.integration.reject_iteration;
    p.integration.reject_target += report.integration.reject_target;
    p.integration.reject_error += report.integration.reject_error;
}

inline void recordBinaryRescanProfile(
    const PS::U64 attempts,
    const PS::U64 skipped,
    const PS::U64 accepted)
{
    if (!FP_t::ks_profile || (attempts == 0 && skipped == 0)) return;
    std::lock_guard<std::mutex> lock(binaryProfileMutex());
    auto& profile = binaryProfileStorage();
    profile.rescan_attempts += attempts;
    profile.rescan_precheck_skips += skipped;
    profile.reactivations += accepted;
}

class BinaryAttemptTimer {
    bool enabled;
    std::chrono::steady_clock::time_point start;
public:
    BinaryAttemptTimer(): enabled(FP_t::ks_profile), start{} {
        if (enabled) start = std::chrono::steady_clock::now();
    }
    PS::F64 elapsed() const {
        if (!enabled) return 0;
        return std::chrono::duration<PS::F64>(std::chrono::steady_clock::now()-start).count();
    }
};

class BinaryFallbackTimer {
public:
    bool enabled;
    PS::F64 excluded_seconds = 0;
    std::chrono::steady_clock::time_point start;
    explicit BinaryFallbackTimer(bool active): enabled(false), start{} {
        if (active) activate();
    }
    void activate() {
        if (!FP_t::ks_profile || enabled) return;
        start = std::chrono::steady_clock::now();
        enabled = true;
    }
    ~BinaryFallbackTimer() {
        if (!enabled) return;
        const PS::F64 seconds = std::chrono::duration<PS::F64>(std::chrono::steady_clock::now()-start).count();
        std::lock_guard<std::mutex> lock(binaryProfileMutex());
        ++binaryProfileStorage().fallback_intervals;
        binaryProfileStorage().fallback_seconds += std::max(PS::F64(0), seconds-excluded_seconds);
    }
};

inline void writeBinaryProfile(std::ostream& out, PS::F64 time, bool header = false)
{
    if (header) {
        out << "time,attempts,accepted,groups,prefixes,blocks,coupling_iterations,full_force_evaluations,"
            << "ks_step_calls,ks_corrections,ks_force_evaluations,ks_accepted_steps,ks_rejected_steps,"
            << "fallback_intervals,rescan_attempts,reactivations,attempt_seconds,fallback_seconds";
        for (int i = 1; i < static_cast<int>(BinaryFailure::count); ++i)
            out << ",reject_" << binaryFailureName(static_cast<BinaryFailure>(i));
        out << ",screened_candidates,ks_reject_iteration,ks_reject_target,ks_reject_error,rescan_precheck_skips,omitted_tide_samples,max_omitted_tide_ratio\n";
    }
    const auto p = getBinaryProfile();
    out << time << ',' << p.attempts << ',' << p.accepted << ',' << p.groups << ',' << p.prefixes
        << ',' << p.blocks << ',' << p.coupling_iterations << ',' << p.full_force_evaluations
        << ',' << p.integration.step_calls << ',' << p.integration.correction_iterations
        << ',' << p.integration.force_evaluations << ',' << p.integration.accepted_steps
        << ',' << p.integration.rejected_steps << ',' << p.fallback_intervals
        << ',' << p.rescan_attempts << ',' << p.reactivations
        << ',' << p.attempt_seconds << ',' << p.fallback_seconds;
    for (int i = 1; i < static_cast<int>(BinaryFailure::count); ++i) out << ',' << p.reasons[i];
    out << ',' << p.screened_candidates << ',' << p.integration.reject_iteration
        << ',' << p.integration.reject_target << ',' << p.integration.reject_error << ',' << p.rescan_precheck_skips << ',' << p.omitted_tide_samples
        << ',' << p.max_omitted_tide_ratio << '\n';
    out.flush();
}

// avoid writing the same time twice
//
class BinaryProfileOutput {
    std::ostream& out;
    PS::F64 last_time;
public:
    explicit BinaryProfileOutput(std::ostream& stream,
        PS::F64 previous_time = std::numeric_limits<PS::F64>::quiet_NaN())
        : out(stream), last_time(previous_time) {}

    void write(PS::F64 time, bool header = false)
    {
        if (time == last_time) return;
        writeBinaryProfile(out, time, header);
        if (out) last_time = time;
    }
};
