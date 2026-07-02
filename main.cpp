#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef WITH_CPLEX
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#endif

using std::cerr;
using std::cout;
using std::endl;
using std::size_t;

namespace {

constexpr double EPS = 1e-9;

std::ofstream g_bgmh_trace;

double env_double(const char* name, double fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || raw[0] == '\0') return fallback;
    try {
        return std::stod(raw);
    } catch (...) {
        return fallback;
    }
}

int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || raw[0] == '\0') return fallback;
    try {
        return std::stoi(raw);
    } catch (...) {
        return fallback;
    }
}

struct ExperimentConfig {
    int n = 50;
    double theta = 0.5;
    double due_range = 0.8;
    double sigma = 3.0;
    std::string due_mode = "taur";
    int instance_id = 0;
    int p_min = 2;
    int p_max = 30;
    int e_min = 2;
    int e_max = 50;
    double alpha = 10.0;
    double beta = 1.0;
    long long seed_base = 260985;
    bool common_base_seed = false;
    double due_reference_sigma = -1.0;
    double recharge_rho = -1.0;
    std::string pe_correlation = "independent";
    std::string due_reference_correlation = "";
};

struct Instance {
    std::vector<std::pair<int, int>> jobs;
    std::vector<int> due_dates;
    int capacity = 0;
    double alpha = 5.0;
    double beta = 1.0;
    ExperimentConfig config;
};

void trace_bgmh_event(
    const Instance& inst,
    const std::string& phase,
    int round,
    double time,
    int objective
) {
    const char* path = std::getenv("BGMH_TRACE_FILE");
    if (path == nullptr || path[0] == '\0') return;
    if (!g_bgmh_trace.is_open()) {
        bool exists = static_cast<bool>(std::ifstream(path));
        g_bgmh_trace.open(path, std::ios::app);
        if (!g_bgmh_trace) return;
        if (!exists) {
            g_bgmh_trace
                << "n,tau,R,gamma,rep,phase,round,time,objective\n";
        }
    }
    const auto& cfg = inst.config;
    g_bgmh_trace << std::fixed << std::setprecision(6)
                 << cfg.n << ','
                 << cfg.theta << ','
                 << cfg.due_range << ','
                 << cfg.sigma << ','
                 << cfg.instance_id << ','
                 << phase << ','
                 << round << ','
                 << time << ','
                 << objective << '\n';
    g_bgmh_trace.flush();
}

struct TypedInstance {
    std::vector<int> p;
    std::vector<int> e;
    std::vector<int> count;
    std::vector<int> due_dates;
    int capacity = 0;
    double alpha = 5.0;
    double beta = 1.0;
    int n = 0;
    int k_types = 0;
};

struct SolveResult {
    std::string algorithm;
    int objective = 0;
    double best_bound = 0.0;
    bool certified = false;
    long long nodes = 0;
    double runtime = 0.0;
    int incumbent_initial = 0;
    double time_limit = 0.0;
    bool available = true;
    int recharges = -1;
    double makespan = -1.0;
    double total_recharge_time = -1.0;
    double battery_utilization = -1.0;
};

long long instance_seed(const ExperimentConfig& cfg) {
    long long scaled_theta = static_cast<long long>(std::llround(cfg.theta * 1000.0));
    long long scaled_range = static_cast<long long>(std::llround(cfg.due_range * 1000.0));
    long long scaled_sigma = static_cast<long long>(std::llround(cfg.sigma * 1000.0));
    long long scaled_alpha = static_cast<long long>(std::llround(cfg.alpha * 1000.0));
    long long scaled_beta = static_cast<long long>(std::llround(cfg.beta * 1000.0));
    long long parameter_component = cfg.common_base_seed
        ? 577LL
        : 37LL * scaled_theta + 101LL * scaled_range + scaled_sigma
            + 541LL * scaled_alpha + 547LL * scaled_beta;
    long long correlation_component = cfg.common_base_seed
        ? 0LL
        : static_cast<long long>(std::hash<std::string>{}(cfg.pe_correlation) % 1000003ULL);
    return cfg.seed_base + 1000003LL * cfg.instance_id + 10007LL * cfg.n
        + parameter_component
        + correlation_component
        + 503LL * cfg.p_min + 509LL * cfg.p_max
        + 521LL * cfg.e_min + 523LL * cfg.e_max;
}

double heuristic_makespan(
    const std::vector<std::pair<int, int>>& jobs,
    int capacity,
    double alpha,
    double beta
);

Instance generate_instance(const ExperimentConfig& cfg) {
    if (cfg.p_min <= 0 || cfg.p_max < cfg.p_min) {
        throw std::runtime_error("Invalid processing-time range");
    }
    if (cfg.e_min <= 0 || cfg.e_max < cfg.e_min) {
        throw std::runtime_error("Invalid energy range");
    }
    std::mt19937_64 rng(static_cast<std::uint64_t>(instance_seed(cfg)));
    std::uniform_int_distribution<int> p_dist(cfg.p_min, cfg.p_max);
    std::uniform_int_distribution<int> e_dist(cfg.e_min, cfg.e_max);

    Instance inst;
    inst.config = cfg;
    inst.alpha = cfg.alpha;
    inst.beta = cfg.beta;
    inst.jobs.reserve(cfg.n);
    std::vector<int> p_values;
    std::vector<int> e_values;
    p_values.reserve(cfg.n);
    e_values.reserve(cfg.n);
    for (int i = 0; i < cfg.n; ++i) {
        p_values.push_back(p_dist(rng));
        e_values.push_back(e_dist(rng));
    }
    auto paired_jobs = [&](const std::string& mode) {
        std::vector<int> p_copy = p_values;
        std::vector<int> e_copy = e_values;
        if (mode == "positive" || mode == "pos") {
            std::sort(p_copy.begin(), p_copy.end());
            std::sort(e_copy.begin(), e_copy.end());
        } else if (mode == "negative" || mode == "neg") {
            std::sort(p_copy.begin(), p_copy.end());
            std::sort(e_copy.begin(), e_copy.end(), std::greater<int>());
        } else if (mode == "independent" || mode == "ind") {
            std::sort(p_copy.begin(), p_copy.end());
            std::sort(e_copy.begin(), e_copy.end());
            std::mt19937_64 pair_rng(static_cast<std::uint64_t>(instance_seed(cfg) + 900001LL));
            std::shuffle(e_copy.begin(), e_copy.end(), pair_rng);
        } else {
            throw std::runtime_error("Unknown processing-energy correlation mode: " + mode);
        }
        std::vector<std::pair<int, int>> out;
        out.reserve(cfg.n);
        for (int i = 0; i < cfg.n; ++i) out.emplace_back(p_copy[i], e_copy[i]);
        return out;
    };
    inst.jobs = paired_jobs(cfg.pe_correlation);
    std::shuffle(inst.jobs.begin(), inst.jobs.end(), rng);

    int max_e = 0, min_p = std::numeric_limits<int>::max();
    int sum_p = 0, sum_e = 0;
    for (const auto& job : inst.jobs) {
        max_e = std::max(max_e, job.second);
        min_p = std::min(min_p, job.first);
        sum_p += job.first;
        sum_e += job.second;
    }
    inst.capacity = std::max(max_e, static_cast<int>(std::ceil(cfg.sigma * max_e)));
    if (cfg.recharge_rho > 0.0 && cfg.recharge_rho < 1.0) {
        inst.alpha = cfg.recharge_rho / (1.0 - cfg.recharge_rho) * cfg.beta * inst.capacity;
        inst.config.alpha = inst.alpha;
    }
    int d_min = 1;
    int d_max = 1;
    if (cfg.due_mode == "taur" || cfg.due_mode == "tau-r") {
        int due_capacity = inst.capacity;
        if (cfg.due_reference_sigma > 0.0) {
            due_capacity = std::max(
                max_e,
                static_cast<int>(std::ceil(cfg.due_reference_sigma * max_e))
            );
        }
        std::vector<std::pair<int, int>> reference_jobs = inst.jobs;
        if (!cfg.due_reference_correlation.empty()) {
            reference_jobs = paired_jobs(cfg.due_reference_correlation);
        }
        double horizon = heuristic_makespan(reference_jobs, due_capacity, inst.alpha, cfg.beta);
        d_min = static_cast<int>(std::ceil(horizon * (1.0 - cfg.theta - cfg.due_range / 2.0)));
        d_max = static_cast<int>(std::floor(horizon * (1.0 - cfg.theta + cfg.due_range / 2.0)));
    } else if (cfg.due_mode == "theta") {
        double full_recharge_time = inst.alpha + cfg.beta * inst.capacity;
        int min_blocks = static_cast<int>(std::ceil(static_cast<double>(sum_e) / inst.capacity));
        d_min = static_cast<int>(std::ceil(full_recharge_time + min_p));
        d_max = static_cast<int>(
            std::floor(cfg.theta * (full_recharge_time * min_blocks + sum_p))
        );
    } else {
        throw std::runtime_error("Unknown due-date mode: " + cfg.due_mode);
    }
    d_min = std::max(1, d_min);
    d_max = std::max(d_min, d_max);
    std::uniform_int_distribution<int> d_dist(d_min, d_max);
    inst.due_dates.reserve(cfg.n);
    for (int i = 0; i < cfg.n; ++i) inst.due_dates.push_back(d_dist(rng));
    std::sort(inst.due_dates.begin(), inst.due_dates.end());
    return inst;
}

TypedInstance to_typed_instance(const Instance& inst) {
    std::map<std::pair<int, int>, int> counts;
    for (const auto& job : inst.jobs) counts[job]++;
    TypedInstance t;
    for (const auto& kv : counts) {
        t.p.push_back(kv.first.first);
        t.e.push_back(kv.first.second);
        t.count.push_back(kv.second);
    }
    t.due_dates = inst.due_dates;
    std::sort(t.due_dates.begin(), t.due_dates.end());
    t.capacity = inst.capacity;
    t.alpha = inst.alpha;
    t.beta = inst.beta;
    t.n = static_cast<int>(inst.jobs.size());
    t.k_types = static_cast<int>(t.p.size());
    return t;
}

int sum_vec(const std::vector<int>& v) {
    return std::accumulate(v.begin(), v.end(), 0);
}

double heuristic_makespan(const std::vector<std::pair<int, int>>& jobs, int capacity, double alpha, double beta) {
    std::vector<std::pair<int, int>> sorted_jobs = jobs;
    std::sort(sorted_jobs.begin(), sorted_jobs.end());
    double best = std::numeric_limits<double>::infinity();

    for (int mode = 0; mode < 2; ++mode) {
        std::vector<int> block_e;
        std::vector<std::vector<std::pair<int, int>>> block_jobs;
        for (const auto& job : sorted_jobs) {
            int chosen = -1;
            if (mode == 0) {
                for (int b = 0; b < static_cast<int>(block_e.size()); ++b) {
                    if (block_e[b] + job.second <= capacity) {
                        chosen = b;
                        break;
                    }
                }
            } else {
                for (int b = 0; b < static_cast<int>(block_e.size()); ++b) {
                    if (block_e[b] + job.second <= capacity
                        && (chosen < 0 || block_e[b] > block_e[chosen])) {
                        chosen = b;
                    }
                }
            }
            if (chosen < 0) {
                block_e.push_back(job.second);
                block_jobs.push_back({job});
            } else {
                block_e[chosen] += job.second;
                block_jobs[chosen].push_back(job);
            }
        }

        double time = 0.0;
        for (int b = 0; b < static_cast<int>(block_jobs.size()); ++b) {
            time += alpha + beta * block_e[b];
            std::sort(block_jobs[b].begin(), block_jobs[b].end());
            for (const auto& job : block_jobs[b]) time += job.first;
        }
        best = std::min(best, time);
    }
    return best;
}

struct ScheduleEval {
    int late = 0;
    double makespan = 0.0;
};

struct ScheduleMetrics {
    int recharges = -1;
    double makespan = -1.0;
    double total_recharge_time = -1.0;
    double battery_utilization = -1.0;
};

std::vector<int> flatten_blocks(const std::vector<std::vector<int>>& blocks) {
    std::vector<int> sequence;
    for (const auto& block : blocks) {
        sequence.insert(sequence.end(), block.begin(), block.end());
    }
    return sequence;
}

int block_energy(const Instance& inst, const std::vector<int>& block) {
    int energy = 0;
    for (int j : block) energy += inst.jobs[j].second;
    return energy;
}

void sort_block_spt(const Instance& inst, std::vector<int>& block) {
    std::sort(block.begin(), block.end(), [&](int a, int b) {
        if (inst.jobs[a].first != inst.jobs[b].first) return inst.jobs[a].first < inst.jobs[b].first;
        if (inst.jobs[a].second != inst.jobs[b].second) return inst.jobs[a].second < inst.jobs[b].second;
        return a < b;
    });
}

void normalize_blocks(const Instance& inst, std::vector<std::vector<int>>& blocks) {
    for (auto& block : blocks) sort_block_spt(inst, block);
    blocks.erase(
        std::remove_if(blocks.begin(), blocks.end(), [](const auto& block) { return block.empty(); }),
        blocks.end()
    );
}

void validate_blocks(const Instance& inst, const std::vector<std::vector<int>>& blocks) {
    std::vector<int> seen(inst.jobs.size(), 0);
    for (const auto& block : blocks) {
        int energy = 0;
        for (int j : block) {
            if (j < 0 || j >= static_cast<int>(inst.jobs.size())) {
                throw std::runtime_error("Heuristic produced an invalid job index");
            }
            ++seen[j];
            energy += inst.jobs[j].second;
        }
        if (energy > inst.capacity) {
            throw std::runtime_error("Heuristic produced an over-capacity block");
        }
    }
    for (int c : seen) {
        if (c != 1) throw std::runtime_error("Heuristic produced a duplicate or missing job");
    }
}

double block_duration(const Instance& inst, const std::vector<int>& block) {
    int energy = 0;
    int processing = 0;
    for (int j : block) {
        processing += inst.jobs[j].first;
        energy += inst.jobs[j].second;
    }
    return inst.alpha + inst.beta * energy + processing;
}

void reorder_blocks(const Instance& inst, std::vector<std::vector<int>>& blocks, int mode) {
    if (mode == 0) return;
    std::stable_sort(blocks.begin(), blocks.end(), [&](const auto& a, const auto& b) {
        auto key = [&](const std::vector<int>& block) {
            int energy = 0;
            int processing = 0;
            int min_p = std::numeric_limits<int>::max();
            for (int j : block) {
                processing += inst.jobs[j].first;
                energy += inst.jobs[j].second;
                min_p = std::min(min_p, inst.jobs[j].first);
            }
            if (mode == 1) return inst.alpha + inst.beta * energy + processing;
            if (mode == 2) return static_cast<double>(processing) + inst.beta * energy;
            if (mode == 3) return static_cast<double>(energy);
            return static_cast<double>(min_p);
        };
        double ka = key(a);
        double kb = key(b);
        if (std::abs(ka - kb) > EPS) return ka < kb;
        return a.size() < b.size();
    });
}

ScheduleEval evaluate_blocks(const Instance& inst, const std::vector<std::vector<int>>& blocks) {
    double time = 0.0;
    int pos = 0;
    int late = 0;
    for (const auto& block : blocks) {
        int energy = block_energy(inst, block);
        if (energy > inst.capacity) return {static_cast<int>(inst.jobs.size()) + 1, 1e100};
        time += inst.alpha + inst.beta * energy;
        for (int j : block) {
            time += inst.jobs[j].first;
            if (time > inst.due_dates[pos] + EPS) ++late;
            ++pos;
        }
    }
    return {late, time};
}

ScheduleEval evaluate_blocks_full_recharge(
    const Instance& inst,
    const std::vector<std::vector<int>>& blocks
) {
    double time = 0.0;
    int pos = 0;
    int late = 0;
    const double full_recharge_time = inst.alpha + inst.beta * inst.capacity;
    for (const auto& block : blocks) {
        int energy = block_energy(inst, block);
        if (energy > inst.capacity) return {static_cast<int>(inst.jobs.size()) + 1, 1e100};
        time += full_recharge_time;
        for (int j : block) {
            time += inst.jobs[j].first;
            if (time > inst.due_dates[pos] + EPS) ++late;
            ++pos;
        }
    }
    return {late, time};
}

ScheduleMetrics schedule_metrics_from_blocks(
    const Instance& inst,
    const std::vector<std::vector<int>>& blocks,
    bool full_recharge = false
) {
    ScheduleMetrics metrics;
    metrics.recharges = static_cast<int>(blocks.size());
    metrics.total_recharge_time = 0.0;
    double utilization_sum = 0.0;
    for (const auto& block : blocks) {
        int energy = block_energy(inst, block);
        metrics.total_recharge_time += full_recharge
            ? inst.alpha + inst.beta * inst.capacity
            : inst.alpha + inst.beta * energy;
        utilization_sum += static_cast<double>(energy) / std::max(1, inst.capacity);
    }
    metrics.battery_utilization = blocks.empty()
        ? 0.0
        : utilization_sum / static_cast<double>(blocks.size());
    metrics.makespan = full_recharge
        ? evaluate_blocks_full_recharge(inst, blocks).makespan
        : evaluate_blocks(inst, blocks).makespan;
    return metrics;
}

void attach_schedule_metrics(
    SolveResult& res,
    const Instance& inst,
    const std::vector<std::vector<int>>& blocks,
    bool full_recharge = false
) {
    ScheduleMetrics metrics = schedule_metrics_from_blocks(inst, blocks, full_recharge);
    res.recharges = metrics.recharges;
    res.makespan = metrics.makespan;
    res.total_recharge_time = metrics.total_recharge_time;
    res.battery_utilization = metrics.battery_utilization;
}

bool better_schedule(const ScheduleEval& a, const ScheduleEval& b) {
    if (a.late != b.late) return a.late < b.late;
    return a.makespan < b.makespan - EPS;
}

std::vector<double> completion_times(const Instance& inst, const std::vector<std::vector<int>>& blocks) {
    std::vector<double> out;
    out.reserve(inst.jobs.size());
    double time = 0.0;
    for (const auto& block : blocks) {
        int energy = block_energy(inst, block);
        time += inst.alpha + inst.beta * energy;
        for (int j : block) {
            time += inst.jobs[j].first;
            out.push_back(time);
        }
    }
    return out;
}

std::vector<double> completion_times_full_recharge(
    const Instance& inst,
    const std::vector<std::vector<int>>& blocks
) {
    std::vector<double> out;
    out.reserve(inst.jobs.size());
    double time = 0.0;
    const double full_recharge_time = inst.alpha + inst.beta * inst.capacity;
    for (const auto& block : blocks) {
        time += full_recharge_time;
        for (int j : block) {
            time += inst.jobs[j].first;
            out.push_back(time);
        }
    }
    return out;
}

std::vector<std::vector<int>> best_blocks_for_fixed_sequence(const Instance& inst, const std::vector<int>& sequence) {
    const int n = static_cast<int>(sequence.size());
    const double inf = 1e100;
    std::vector<std::vector<double>> dp(n + 1, std::vector<double>(n + 1, inf));
    std::vector<std::vector<int>> parent_start(n + 1, std::vector<int>(n + 1, -1));
    std::vector<std::vector<int>> parent_late(n + 1, std::vector<int>(n + 1, -1));
    dp[0][0] = 0.0;

    for (int start = 0; start < n; ++start) {
        for (int late_so_far = 0; late_so_far <= start; ++late_so_far) {
            double base_time = dp[start][late_so_far];
            if (base_time >= inf / 2) continue;
            int energy = 0;
            for (int end = start; end < n; ++end) {
                int job = sequence[end];
                energy += inst.jobs[job].second;
                if (energy > inst.capacity) break;
                double time = base_time + inst.alpha + inst.beta * energy;
                int segment_late = 0;
                for (int pos = start; pos <= end; ++pos) {
                    time += inst.jobs[sequence[pos]].first;
                    if (time > inst.due_dates[pos] + EPS) ++segment_late;
                }
                int total_late = late_so_far + segment_late;
                if (time + EPS < dp[end + 1][total_late]) {
                    dp[end + 1][total_late] = time;
                    parent_start[end + 1][total_late] = start;
                    parent_late[end + 1][total_late] = late_so_far;
                }
            }
        }
    }

    int best_late = -1;
    for (int late = 0; late <= n; ++late) {
        if (dp[n][late] < inf / 2) {
            best_late = late;
            break;
        }
    }
    if (best_late < 0) throw std::runtime_error("No feasible segmentation for sequence");

    std::vector<std::vector<int>> blocks;
    int end = n;
    int late = best_late;
    while (end > 0) {
        int start = parent_start[end][late];
        int prev_late = parent_late[end][late];
        if (start < 0 || prev_late < 0) throw std::runtime_error("Invalid segmentation parent");
        blocks.emplace_back(sequence.begin() + start, sequence.begin() + end);
        end = start;
        late = prev_late;
    }
    std::reverse(blocks.begin(), blocks.end());
    normalize_blocks(inst, blocks);
    return blocks;
}

std::vector<std::vector<int>> best_blocks_for_fixed_sequence_full_recharge(
    const Instance& inst,
    const std::vector<int>& sequence
) {
    const int n = static_cast<int>(sequence.size());
    const double inf = 1e100;
    const double full_recharge_time = inst.alpha + inst.beta * inst.capacity;
    std::vector<std::vector<double>> dp(n + 1, std::vector<double>(n + 1, inf));
    std::vector<std::vector<int>> parent_start(n + 1, std::vector<int>(n + 1, -1));
    std::vector<std::vector<int>> parent_late(n + 1, std::vector<int>(n + 1, -1));
    dp[0][0] = 0.0;

    for (int start = 0; start < n; ++start) {
        for (int late_so_far = 0; late_so_far <= start; ++late_so_far) {
            double base_time = dp[start][late_so_far];
            if (base_time >= inf / 2) continue;
            int energy = 0;
            for (int end = start; end < n; ++end) {
                int job = sequence[end];
                energy += inst.jobs[job].second;
                if (energy > inst.capacity) break;
                double time = base_time + full_recharge_time;
                int segment_late = 0;
                for (int pos = start; pos <= end; ++pos) {
                    time += inst.jobs[sequence[pos]].first;
                    if (time > inst.due_dates[pos] + EPS) ++segment_late;
                }
                int total_late = late_so_far + segment_late;
                if (time + EPS < dp[end + 1][total_late]) {
                    dp[end + 1][total_late] = time;
                    parent_start[end + 1][total_late] = start;
                    parent_late[end + 1][total_late] = late_so_far;
                }
            }
        }
    }

    int best_late = -1;
    for (int late = 0; late <= n; ++late) {
        if (dp[n][late] < inf / 2) {
            best_late = late;
            break;
        }
    }
    if (best_late < 0) throw std::runtime_error("No feasible full-recharge segmentation for sequence");

    std::vector<std::vector<int>> blocks;
    int end = n;
    int late = best_late;
    while (end > 0) {
        int start = parent_start[end][late];
        int prev_late = parent_late[end][late];
        if (start < 0 || prev_late < 0) throw std::runtime_error("Invalid full-recharge segmentation parent");
        blocks.emplace_back(sequence.begin() + start, sequence.begin() + end);
        end = start;
        late = prev_late;
    }
    std::reverse(blocks.begin(), blocks.end());
    normalize_blocks(inst, blocks);
    return blocks;
}

std::vector<std::vector<int>> construct_blocks(
    const Instance& inst, const std::vector<int>& order, int pack_mode, int reorder_mode
) {
    std::vector<int> block_e;
    std::vector<std::vector<int>> blocks;
    for (int j : order) {
        int chosen = -1;
        if (pack_mode == 0) {
            if (!block_e.empty() && block_e.back() + inst.jobs[j].second <= inst.capacity) {
                chosen = static_cast<int>(block_e.size()) - 1;
            }
        } else {
            for (int b = 0; b < static_cast<int>(block_e.size()); ++b) {
                if (block_e[b] + inst.jobs[j].second > inst.capacity) continue;
                if (chosen < 0) {
                    chosen = b;
                } else if (pack_mode == 1) {
                    break;
                } else if (pack_mode == 2 && block_e[b] > block_e[chosen]) {
                    chosen = b;
                } else if (pack_mode == 3 && block_e[b] < block_e[chosen]) {
                    chosen = b;
                }
            }
        }
        if (chosen < 0) {
            block_e.push_back(inst.jobs[j].second);
            blocks.push_back({j});
        } else {
            block_e[chosen] += inst.jobs[j].second;
            blocks[chosen].push_back(j);
        }
    }
    normalize_blocks(inst, blocks);
    reorder_blocks(inst, blocks, reorder_mode);
    return blocks;
}

void improve_block_order(const Instance& inst, std::vector<std::vector<int>>& blocks) {
    ScheduleEval best_eval = evaluate_blocks(inst, blocks);
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i + 1 < static_cast<int>(blocks.size()); ++i) {
            std::swap(blocks[i], blocks[i + 1]);
            ScheduleEval cand = evaluate_blocks(inst, blocks);
            if (better_schedule(cand, best_eval)) {
                best_eval = cand;
                improved = true;
            } else {
                std::swap(blocks[i], blocks[i + 1]);
            }
        }
    }
}

void local_improve_blocks(const Instance& inst, std::vector<std::vector<int>>& blocks) {
    normalize_blocks(inst, blocks);
    ScheduleEval best_eval = evaluate_blocks(inst, blocks);

    for (int pass = 0; pass < 3; ++pass) {
        bool changed = false;

        std::vector<int> energies(blocks.size());
        for (int b = 0; b < static_cast<int>(blocks.size()); ++b) energies[b] = block_energy(inst, blocks[b]);

        std::vector<std::vector<int>> best_blocks;
        ScheduleEval move_eval = best_eval;
        for (int a = 0; a < static_cast<int>(blocks.size()); ++a) {
            for (int ia = 0; ia < static_cast<int>(blocks[a].size()); ++ia) {
                int job = blocks[a][ia];
                for (int b = 0; b < static_cast<int>(blocks.size()); ++b) {
                    if (a == b || energies[b] + inst.jobs[job].second > inst.capacity) continue;
                    auto cand = blocks;
                    cand[a].erase(cand[a].begin() + ia);
                    cand[b].push_back(job);
                    normalize_blocks(inst, cand);
                    improve_block_order(inst, cand);
                    ScheduleEval cand_eval = evaluate_blocks(inst, cand);
                    if (better_schedule(cand_eval, move_eval)) {
                        move_eval = cand_eval;
                        best_blocks = std::move(cand);
                    }
                }
            }
        }
        if (!best_blocks.empty()) {
            blocks = std::move(best_blocks);
            best_eval = move_eval;
            changed = true;
        }

        energies.assign(blocks.size(), 0);
        for (int b = 0; b < static_cast<int>(blocks.size()); ++b) energies[b] = block_energy(inst, blocks[b]);

        std::vector<std::vector<int>> swap_blocks;
        ScheduleEval swap_eval = best_eval;
        for (int a = 0; a < static_cast<int>(blocks.size()); ++a) {
            for (int b = a + 1; b < static_cast<int>(blocks.size()); ++b) {
                for (int ia = 0; ia < static_cast<int>(blocks[a].size()); ++ia) {
                    for (int ib = 0; ib < static_cast<int>(blocks[b].size()); ++ib) {
                        int ja = blocks[a][ia];
                        int jb = blocks[b][ib];
                        if (energies[a] - inst.jobs[ja].second + inst.jobs[jb].second > inst.capacity) continue;
                        if (energies[b] - inst.jobs[jb].second + inst.jobs[ja].second > inst.capacity) continue;
                        auto cand = blocks;
                        cand[a][ia] = jb;
                        cand[b][ib] = ja;
                        normalize_blocks(inst, cand);
                        improve_block_order(inst, cand);
                        ScheduleEval cand_eval = evaluate_blocks(inst, cand);
                        if (better_schedule(cand_eval, swap_eval)) {
                            swap_eval = cand_eval;
                            swap_blocks = std::move(cand);
                        }
                    }
                }
            }
        }
        if (!swap_blocks.empty()) {
            blocks = std::move(swap_blocks);
            best_eval = swap_eval;
            changed = true;
        }

        if (!changed) break;
    }
}

std::vector<std::vector<int>> block_packing_improvement_heuristic(
    const Instance& inst,
    bool use_sequence_resegmentation = true
) {
    const int n = static_cast<int>(inst.jobs.size());
    std::vector<int> base(n);
    std::iota(base.begin(), base.end(), 0);

    std::vector<std::vector<int>> orders;
    auto add_order = [&](auto cmp) {
        std::vector<int> order = base;
        std::stable_sort(order.begin(), order.end(), cmp);
        orders.push_back(std::move(order));
    };
    const double fac = (inst.alpha + inst.beta * inst.capacity) / inst.capacity;
    add_order([&](int a, int b) {
        if (inst.jobs[a].first != inst.jobs[b].first) return inst.jobs[a].first < inst.jobs[b].first;
        if (inst.jobs[a].second != inst.jobs[b].second) return inst.jobs[a].second < inst.jobs[b].second;
        return a < b;
    });
    add_order([&](int a, int b) {
        double ka = inst.jobs[a].first + fac * inst.jobs[a].second;
        double kb = inst.jobs[b].first + fac * inst.jobs[b].second;
        if (std::abs(ka - kb) > EPS) return ka < kb;
        return inst.jobs[a].first < inst.jobs[b].first;
    });
    add_order([&](int a, int b) {
        double ka = inst.jobs[a].first + inst.beta * inst.jobs[a].second;
        double kb = inst.jobs[b].first + inst.beta * inst.jobs[b].second;
        if (std::abs(ka - kb) > EPS) return ka < kb;
        return inst.jobs[a].second < inst.jobs[b].second;
    });
    add_order([&](int a, int b) {
        if (inst.jobs[a].second != inst.jobs[b].second) return inst.jobs[a].second < inst.jobs[b].second;
        return inst.jobs[a].first < inst.jobs[b].first;
    });
    add_order([&](int a, int b) {
        if (inst.jobs[a].second != inst.jobs[b].second) return inst.jobs[a].second > inst.jobs[b].second;
        return inst.jobs[a].first < inst.jobs[b].first;
    });
    add_order([&](int a, int b) {
        double ka = static_cast<double>(inst.jobs[a].first + inst.beta * inst.jobs[a].second)
            / std::max(1, inst.jobs[a].second);
        double kb = static_cast<double>(inst.jobs[b].first + inst.beta * inst.jobs[b].second)
            / std::max(1, inst.jobs[b].second);
        if (std::abs(ka - kb) > EPS) return ka < kb;
        return inst.jobs[a].first < inst.jobs[b].first;
    });

    struct Candidate {
        ScheduleEval eval;
        std::vector<std::vector<int>> blocks;
    };
    std::vector<Candidate> candidates;
    for (const auto& order : orders) {
        for (int pack_mode = 0; pack_mode < 4; ++pack_mode) {
            for (int reorder_mode = 0; reorder_mode < 5; ++reorder_mode) {
                auto blocks = construct_blocks(inst, order, pack_mode, reorder_mode);
                candidates.push_back({evaluate_blocks(inst, blocks), std::move(blocks)});
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
        return better_schedule(a.eval, b.eval);
    });

    const int improve_count = std::min(12, static_cast<int>(candidates.size()));
    for (int i = 0; i < improve_count; ++i) {
        local_improve_blocks(inst, candidates[i].blocks);
        candidates[i].eval = evaluate_blocks(inst, candidates[i].blocks);
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
        return better_schedule(a.eval, b.eval);
    });

    if (use_sequence_resegmentation) {
        const int resegment_count = std::min(
            static_cast<int>(candidates.size()),
            n >= 100 ? 10 : 24
        );
        for (int i = 0; i < resegment_count; ++i) {
            auto sequence = flatten_blocks(candidates[i].blocks);
            auto resegmented = best_blocks_for_fixed_sequence(inst, sequence);
            ScheduleEval reseg_eval = evaluate_blocks(inst, resegmented);
            if (better_schedule(reseg_eval, candidates[i].eval)) {
                candidates[i].blocks = std::move(resegmented);
                candidates[i].eval = reseg_eval;
            }
        }
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
            return better_schedule(a.eval, b.eval);
        });
    }
    validate_blocks(inst, candidates.front().blocks);
    return candidates.front().blocks;
}

std::vector<int> effective_load_order(const Instance& inst) {
    const int n = static_cast<int>(inst.jobs.size());
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    const double lambda = (inst.alpha + inst.beta * inst.capacity) / inst.capacity;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
        double wa = inst.jobs[a].first + lambda * inst.jobs[a].second;
        double wb = inst.jobs[b].first + lambda * inst.jobs[b].second;
        if (std::abs(wa - wb) > EPS) return wa < wb;
        if (inst.jobs[a].first != inst.jobs[b].first) return inst.jobs[a].first < inst.jobs[b].first;
        if (inst.jobs[a].second != inst.jobs[b].second) return inst.jobs[a].second < inst.jobs[b].second;
        return a < b;
    });
    return order;
}

std::vector<std::vector<int>> rigorous_start_blocks_with_indices(const Instance& inst) {
    return best_blocks_for_fixed_sequence(inst, effective_load_order(inst));
}

std::vector<std::vector<int>> heuristic_blocks_with_indices(const Instance& inst) {
    return rigorous_start_blocks_with_indices(inst);
}

std::vector<std::vector<int>> grasp_effective_load_construction(
    const Instance& inst,
    double construction_time_limit,
    std::uint64_t seed,
    double epsilon = 0.15
) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };

    const int n = static_cast<int>(inst.jobs.size());
    const double lambda = (inst.alpha + inst.beta * inst.capacity) / inst.capacity;
    std::vector<double> weight(n, 0.0);
    for (int j = 0; j < n; ++j) {
        weight[j] = inst.jobs[j].first + lambda * inst.jobs[j].second;
    }

    auto best_blocks = best_blocks_for_fixed_sequence(inst, effective_load_order(inst));
    ScheduleEval best_eval = evaluate_blocks(inst, best_blocks);

    std::mt19937_64 rng(seed);
    int iterations = 0;
    const int min_iterations = std::min(30, std::max(8, n / 10));
    const int max_iterations = std::max(80, std::min(2000, 120000 / std::max(1, n)));
    while ((iterations < min_iterations || elapsed() < construction_time_limit - EPS)
           && iterations < max_iterations) {
        std::vector<int> remaining(n);
        std::iota(remaining.begin(), remaining.end(), 0);
        std::vector<int> seq;
        seq.reserve(n);
        while (!remaining.empty()) {
            double best_w = std::numeric_limits<double>::infinity();
            for (int j : remaining) best_w = std::min(best_w, weight[j]);
            std::vector<int> rcl;
            for (int j : remaining) {
                if (weight[j] <= (1.0 + epsilon) * best_w + EPS) rcl.push_back(j);
            }
            std::uniform_int_distribution<int> pick(0, static_cast<int>(rcl.size()) - 1);
            int chosen = rcl[pick(rng)];
            seq.push_back(chosen);
            remaining.erase(std::find(remaining.begin(), remaining.end(), chosen));
        }
        auto cand_blocks = best_blocks_for_fixed_sequence(inst, seq);
        ScheduleEval cand_eval = evaluate_blocks(inst, cand_blocks);
        if (better_schedule(cand_eval, best_eval)) {
            best_blocks = std::move(cand_blocks);
            best_eval = cand_eval;
        }
        ++iterations;
    }
    validate_blocks(inst, best_blocks);
    return best_blocks;
}

std::vector<int> randomized_effective_load_sequence(
    const Instance& inst,
    std::mt19937_64& rng,
    double epsilon
) {
    const int n = static_cast<int>(inst.jobs.size());
    const double lambda = (inst.alpha + inst.beta * inst.capacity) / inst.capacity;
    std::vector<double> weight(n, 0.0);
    for (int j = 0; j < n; ++j) {
        weight[j] = inst.jobs[j].first + lambda * inst.jobs[j].second;
    }
    std::vector<int> remaining(n);
    std::iota(remaining.begin(), remaining.end(), 0);
    std::vector<int> seq;
    seq.reserve(n);
    while (!remaining.empty()) {
        double best_w = std::numeric_limits<double>::infinity();
        for (int j : remaining) best_w = std::min(best_w, weight[j]);
        std::vector<int> rcl;
        for (int j : remaining) {
            if (weight[j] <= (1.0 + epsilon) * best_w + EPS) rcl.push_back(j);
        }
        std::uniform_int_distribution<int> pick(0, static_cast<int>(rcl.size()) - 1);
        int chosen = rcl[pick(rng)];
        seq.push_back(chosen);
        remaining.erase(std::find(remaining.begin(), remaining.end(), chosen));
    }
    return seq;
}

struct SequenceSolution {
    std::vector<int> sequence;
    std::vector<std::vector<int>> blocks;
    ScheduleEval eval;
};

SequenceSolution evaluate_sequence_solution(
    const Instance& inst,
    const std::vector<int>& sequence,
    bool full_recharge = false
) {
    auto blocks = full_recharge
        ? best_blocks_for_fixed_sequence_full_recharge(inst, sequence)
        : best_blocks_for_fixed_sequence(inst, sequence);
    auto eval = full_recharge
        ? evaluate_blocks_full_recharge(inst, blocks)
        : evaluate_blocks(inst, blocks);
    return {flatten_blocks(blocks), std::move(blocks), eval};
}

int relaxed_binding_position(const Instance& inst) {
    const int n = static_cast<int>(inst.jobs.size());
    auto order = effective_load_order(inst);
    const double lambda = (inst.alpha + inst.beta * inst.capacity) / inst.capacity;
    double relaxed_time = 0.0;
    for (int i = 0; i < n; ++i) {
        int job = order[i];
        relaxed_time += inst.jobs[job].first + lambda * inst.jobs[job].second;
        if (relaxed_time > inst.due_dates[i] + EPS) return i;
    }
    return n - 1;
}

std::vector<int> bound_guided_positions(
    const Instance& inst,
    const std::vector<std::vector<int>>& blocks,
    bool full_recharge = false
) {
    const int n = static_cast<int>(inst.jobs.size());
    const double h_factor = env_double("GRASP_LS_RADIUS_FACTOR", 0.75);
    const int h = std::max(
        2,
        static_cast<int>(std::ceil(h_factor * std::sqrt(static_cast<double>(n))))
    );
    const int critical = relaxed_binding_position(inst);
    std::vector<double> c = full_recharge
        ? completion_times_full_recharge(inst, blocks)
        : completion_times(inst, blocks);
    std::vector<int> ranked(n);
    std::iota(ranked.begin(), ranked.end(), 0);
    std::stable_sort(ranked.begin(), ranked.end(), [&](int a, int b) {
        double slack_a = inst.due_dates[a] - c[a];
        double slack_b = inst.due_dates[b] - c[b];
        bool late_a = slack_a < -EPS;
        bool late_b = slack_b < -EPS;
        if (late_a != late_b) return late_a > late_b;
        if (std::abs(slack_a - slack_b) > EPS) return slack_a < slack_b;
        return std::abs(a - critical) < std::abs(b - critical);
    });

    std::set<int> selected;
    for (int pos = std::max(0, critical - h); pos <= std::min(n - 1, critical + h); ++pos) {
        selected.insert(pos);
    }
    for (int i = 0; i < std::min(n, 2 * h); ++i) selected.insert(ranked[i]);
    std::vector<int> out(selected.begin(), selected.end());
    return out;
}

std::vector<int> insertion_sequence(const std::vector<int>& seq, int from, int to) {
    std::vector<int> out = seq;
    int job = out[from];
    out.erase(out.begin() + from);
    if (to > from) --to;
    to = std::max(0, std::min(to, static_cast<int>(out.size())));
    out.insert(out.begin() + to, job);
    return out;
}

long long bound_guided_local_search(
    const Instance& inst,
    SequenceSolution& sol,
    int root_lb,
    double time_limit,
    bool full_recharge = false
) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    const int n = static_cast<int>(inst.jobs.size());
    long long evaluations = 0;
    bool improved = true;
    while (improved && elapsed() < time_limit - EPS && sol.eval.late > root_lb) {
        improved = false;
        auto positions = bound_guided_positions(inst, sol.blocks, full_recharge);
        std::set<std::pair<int, int>> insertion_moves;
        std::set<int> adjacent_moves;
        for (int from : positions) {
            for (int to : positions) {
                if (from != to) insertion_moves.insert({from, to});
            }
            if (from > 0) adjacent_moves.insert(from - 1);
            if (from + 1 < n) adjacent_moves.insert(from);
        }

        SequenceSolution best_move = sol;
        for (auto [from, to] : insertion_moves) {
            if (elapsed() >= time_limit - EPS) break;
            auto cand_seq = insertion_sequence(sol.sequence, from, to);
            auto cand = evaluate_sequence_solution(inst, cand_seq, full_recharge);
            ++evaluations;
            if (better_schedule(cand.eval, best_move.eval)) {
                best_move = std::move(cand);
            }
        }
        for (int pos : adjacent_moves) {
            if (elapsed() >= time_limit - EPS) break;
            auto cand_seq = sol.sequence;
            std::swap(cand_seq[pos], cand_seq[pos + 1]);
            auto cand = evaluate_sequence_solution(inst, cand_seq, full_recharge);
            ++evaluations;
            if (better_schedule(cand.eval, best_move.eval)) {
                best_move = std::move(cand);
            }
        }
        if (better_schedule(best_move.eval, sol.eval)) {
            sol = std::move(best_move);
            improved = true;
        }
    }
    return evaluations;
}

bool same_sequence(const std::vector<int>& a, const std::vector<int>& b) {
    return a == b;
}

void update_elite_pool(std::vector<SequenceSolution>& elite, SequenceSolution cand, int elite_size) {
    for (const auto& sol : elite) {
        if (same_sequence(sol.sequence, cand.sequence)) return;
    }
    elite.push_back(std::move(cand));
    std::stable_sort(elite.begin(), elite.end(), [](const auto& a, const auto& b) {
        return better_schedule(a.eval, b.eval);
    });
    if (static_cast<int>(elite.size()) > elite_size) elite.resize(elite_size);
}

long long path_relink(
    const Instance& inst,
    const SequenceSolution& start_sol,
    const SequenceSolution& guide_sol,
    SequenceSolution& best_on_path,
    double time_limit,
    bool full_recharge = false
) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    std::vector<int> current = start_sol.sequence;
    const std::vector<int>& guide = guide_sol.sequence;
    const int n = static_cast<int>(current.size());
    std::vector<int> pos(n, 0);
    for (int i = 0; i < n; ++i) pos[current[i]] = i;
    best_on_path = start_sol;
    long long evaluations = 0;

    for (int i = 0; i < n && elapsed() < time_limit - EPS; ++i) {
        if (current[i] == guide[i]) continue;
        int wanted = guide[i];
        int from = pos[wanted];
        int moved = current[from];
        current.erase(current.begin() + from);
        current.insert(current.begin() + i, moved);
        int left = std::min(i, from);
        int right = std::max(i, from);
        for (int k = left; k <= right; ++k) pos[current[k]] = k;
        auto cand = evaluate_sequence_solution(inst, current, full_recharge);
        ++evaluations;
        if (better_schedule(cand.eval, best_on_path.eval)) {
            best_on_path = std::move(cand);
        }
    }
    return evaluations;
}

bool target_set_grasp_rebuild(
    const Instance& inst,
    std::vector<std::vector<int>>& best_blocks,
    ScheduleEval& best_eval,
    int root_lb,
    double time_limit,
    std::uint64_t seed
) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    const int n = static_cast<int>(inst.jobs.size());
    if (best_eval.late <= root_lb || best_eval.late <= 0 || time_limit <= EPS) return false;
    const int target_late = std::max(root_lb, best_eval.late - 1);
    const int target_early = n - target_late;
    const double lambda = (inst.alpha + inst.beta * inst.capacity) / inst.capacity;
    std::vector<double> weight(n, 0.0);
    for (int j = 0; j < n; ++j) {
        weight[j] = inst.jobs[j].first + lambda * inst.jobs[j].second;
    }
    auto by_weight_asc = [&](int a, int b) {
        if (std::abs(weight[a] - weight[b]) > EPS) return weight[a] < weight[b];
        if (inst.jobs[a].first != inst.jobs[b].first) return inst.jobs[a].first < inst.jobs[b].first;
        if (inst.jobs[a].second != inst.jobs[b].second) return inst.jobs[a].second < inst.jobs[b].second;
        return a < b;
    };

    bool improved = false;
    std::mt19937_64 rng(seed);
    int trials = 0;
    const int max_trials = std::max(100, std::min(2000, 120000 / std::max(1, n)));
    while (elapsed() < time_limit - EPS && trials < max_trials && best_eval.late > root_lb) {
        std::vector<int> remaining(n);
        std::iota(remaining.begin(), remaining.end(), 0);
        std::vector<int> early;
        early.reserve(target_early);
        const double epsilon = env_double("GRASP_EPS", 0.15);
        while (static_cast<int>(early.size()) < target_early) {
            double best_w = std::numeric_limits<double>::infinity();
            for (int j : remaining) best_w = std::min(best_w, weight[j]);
            std::vector<int> rcl;
            for (int j : remaining) {
                if (weight[j] <= (1.0 + epsilon) * best_w + EPS) rcl.push_back(j);
            }
            std::uniform_int_distribution<int> pick(0, static_cast<int>(rcl.size()) - 1);
            int chosen = rcl[pick(rng)];
            early.push_back(chosen);
            remaining.erase(std::find(remaining.begin(), remaining.end(), chosen));
        }
        std::stable_sort(early.begin(), early.end(), by_weight_asc);
        std::stable_sort(remaining.begin(), remaining.end(), by_weight_asc);
        std::vector<int> seq = early;
        seq.insert(seq.end(), remaining.begin(), remaining.end());
        auto cand_blocks = best_blocks_for_fixed_sequence(inst, seq);
        ScheduleEval cand_eval = evaluate_blocks(inst, cand_blocks);
        if (better_schedule(cand_eval, best_eval)) {
            best_blocks = std::move(cand_blocks);
            best_eval = cand_eval;
            improved = true;
        }
        ++trials;
    }
    return improved;
}

std::vector<double> expanded_values(
    const TypedInstance& t, const std::vector<int>& rem, const std::vector<double>& type_values
) {
    std::vector<double> values;
    values.reserve(sum_vec(rem));
    for (int i = 0; i < t.k_types; ++i) {
        for (int c = 0; c < rem[i]; ++c) values.push_back(type_values[i]);
    }
    return values;
}

int lb_fractional(const TypedInstance& t, const std::vector<int>& rem, int sc, double t0) {
    std::vector<double> type_values(t.k_types);
    const double fac = (t.alpha + t.beta * t.capacity) / t.capacity;
    for (int i = 0; i < t.k_types; ++i) type_values[i] = t.p[i] + t.e[i] * fac;
    auto values = expanded_values(t, rem, type_values);
    std::sort(values.begin(), values.end());
    double time = t0;
    int late = 0;
    for (int r = 0; r < static_cast<int>(values.size()); ++r) {
        time += values[r];
        if (time > t.due_dates[sc + r] + EPS) ++late;
    }
    return late;
}

int lb_single_setup(const TypedInstance& t, const std::vector<int>& rem, int sc, double t0) {
    std::vector<double> type_values(t.k_types);
    for (int i = 0; i < t.k_types; ++i) type_values[i] = t.p[i] + t.beta * t.e[i];
    auto values = expanded_values(t, rem, type_values);
    std::sort(values.begin(), values.end());
    double time = t0 + t.alpha;
    int late = 0;
    for (int r = 0; r < static_cast<int>(values.size()); ++r) {
        time += values[r];
        if (time > t.due_dates[sc + r] + EPS) ++late;
    }
    return late;
}

int bin_packing_lower_bound(const TypedInstance& t, const std::vector<int>& rem) {
    int total_e = 0;
    std::vector<int> energies;
    energies.reserve(sum_vec(rem));
    for (int i = 0; i < t.k_types; ++i) {
        total_e += rem[i] * t.e[i];
        for (int c = 0; c < rem[i]; ++c) energies.push_back(t.e[i]);
    }
    if (total_e == 0) return 0;
    const int cap = t.capacity;
    int lb = static_cast<int>(std::ceil(static_cast<double>(total_e) / cap));
    std::sort(energies.rbegin(), energies.rend());

    int big = 0;
    for (int e : energies) if (2 * e > cap) ++big;
    lb = std::max(lb, big);

    int clique = energies.empty() ? 0 : 1;
    for (int r = 2; r <= static_cast<int>(energies.size()); ++r) {
        if (energies[r - 2] + energies[r - 1] > cap) clique = r;
        else break;
    }
    lb = std::max(lb, clique);

    std::vector<int> unique_e = energies;
    std::sort(unique_e.begin(), unique_e.end());
    unique_e.erase(std::unique(unique_e.begin(), unique_e.end()), unique_e.end());
    for (int u : unique_e) {
        if (!(0 < 2 * u && 2 * u <= cap)) continue;
        int transformed = 0;
        for (int e : energies) {
            if (e > cap - u) transformed += cap;
            else if (u <= e && e <= cap - u) transformed += e;
        }
        lb = std::max(lb, static_cast<int>(std::ceil(static_cast<double>(transformed) / cap)));
    }
    return lb;
}

std::pair<int, bool> tail_values(const TypedInstance& t, const std::vector<int>& rem, int sc, double t0) {
    int m = sum_vec(rem);
    if (m == 0) return {0, false};
    int rem_p = 0, rem_e = 0;
    std::vector<double> values;
    values.reserve(m);
    for (int i = 0; i < t.k_types; ++i) {
        rem_p += rem[i] * t.p[i];
        rem_e += rem[i] * t.e[i];
        for (int c = 0; c < rem[i]; ++c) values.push_back(t.p[i] + t.beta * t.e[i]);
    }
    int q_min = bin_packing_lower_bound(t, rem);
    double makespan_lb = t0 + rem_p + t.beta * rem_e + t.alpha * q_min;
    std::sort(values.rbegin(), values.rend());

    int forced = 0;
    double cum = 0.0;
    bool last_forced = false;
    for (int tail_count = 0; tail_count < m; ++tail_count) {
        int pos = sc + m - 1 - tail_count;
        double completion_lb = makespan_lb - cum - t.alpha * tail_count;
        if (completion_lb > t.due_dates[pos] + EPS) {
            ++forced;
            if (tail_count == 0) last_forced = true;
        } else {
            break;
        }
        cum += values[tail_count];
    }
    return {forced, last_forced};
}

int lb_tail(const TypedInstance& t, const std::vector<int>& rem, int sc, double t0) {
    return tail_values(t, rem, sc, t0).first;
}

int lb_position_count(const TypedInstance& t, const std::vector<int>& rem, int sc, double t0) {
    int m = sum_vec(rem);
    if (m == 0) return 0;

    std::vector<double> work_values;
    std::vector<int> energies;
    work_values.reserve(m);
    energies.reserve(m);
    for (int i = 0; i < t.k_types; ++i) {
        for (int c = 0; c < rem[i]; ++c) {
            work_values.push_back(t.p[i] + t.beta * t.e[i]);
            energies.push_back(t.e[i]);
        }
    }
    std::sort(work_values.begin(), work_values.end());
    std::sort(energies.begin(), energies.end());

    int late = 0;
    double prefix_work = 0.0;
    int prefix_energy = 0;
    for (int r = 1; r <= m; ++r) {
        prefix_work += work_values[r - 1];
        prefix_energy += energies[r - 1];
        double completion_lb = t0 + prefix_work
            + t.alpha * std::ceil(static_cast<double>(prefix_energy) / t.capacity);
        if (completion_lb > t.due_dates[sc + r - 1] + EPS) ++late;
    }
    return late;
}

int lb_position_dp(
    const TypedInstance& t, const std::vector<int>& rem, int sc, double t0,
    int max_jobs = 35, int max_energy_sum = 500
) {
    int m = sum_vec(rem);
    if (m == 0) return 0;
    int total_e = 0;
    std::vector<std::pair<int, double>> jobs;
    jobs.reserve(m);
    for (int i = 0; i < t.k_types; ++i) {
        total_e += rem[i] * t.e[i];
        for (int c = 0; c < rem[i]; ++c) jobs.push_back({t.e[i], t.p[i] + t.beta * t.e[i]});
    }
    if (m > max_jobs || total_e > max_energy_sum) return 0;

    const double inf = 1e100;
    std::vector<std::vector<double>> dp(m + 1, std::vector<double>(total_e + 1, inf));
    dp[0][0] = 0.0;
    for (auto [e, w] : jobs) {
        for (int r = m - 1; r >= 0; --r) {
            for (int g = 0; g + e <= total_e; ++g) {
                double val = dp[r][g];
                if (val < inf) dp[r + 1][g + e] = std::min(dp[r + 1][g + e], val + w);
            }
        }
    }

    int late = 0;
    for (int r = 1; r <= m; ++r) {
        double best = inf;
        for (int g = 1; g <= total_e; ++g) {
            if (dp[r][g] < inf) {
                best = std::min(best, dp[r][g] + t.alpha * std::ceil(static_cast<double>(g) / t.capacity));
            }
        }
        if (t0 + best > t.due_dates[sc + r - 1] + EPS) ++late;
    }
    return late;
}

int lb_forced_position_union(const TypedInstance& t, const std::vector<int>& rem, int sc, double t0) {
    int m = sum_vec(rem);
    if (m == 0) return 0;

    std::vector<char> forced(m, 0);

    std::vector<double> fractional_values;
    std::vector<double> work_values;
    std::vector<int> energies;
    fractional_values.reserve(m);
    work_values.reserve(m);
    energies.reserve(m);

    const double fac = (t.alpha + t.beta * t.capacity) / t.capacity;
    int rem_p = 0;
    int rem_e = 0;
    for (int i = 0; i < t.k_types; ++i) {
        rem_p += rem[i] * t.p[i];
        rem_e += rem[i] * t.e[i];
        for (int c = 0; c < rem[i]; ++c) {
            fractional_values.push_back(t.p[i] + t.e[i] * fac);
            work_values.push_back(t.p[i] + t.beta * t.e[i]);
            energies.push_back(t.e[i]);
        }
    }

    std::sort(fractional_values.begin(), fractional_values.end());
    double fractional_time = t0;
    for (int r = 0; r < m; ++r) {
        fractional_time += fractional_values[r];
        if (fractional_time > t.due_dates[sc + r] + EPS) forced[r] = 1;
    }

    std::sort(work_values.begin(), work_values.end());
    double single_setup_time = t0 + t.alpha;
    for (int r = 0; r < m; ++r) {
        single_setup_time += work_values[r];
        if (single_setup_time > t.due_dates[sc + r] + EPS) forced[r] = 1;
    }

    std::sort(energies.begin(), energies.end());
    double prefix_work = 0.0;
    int prefix_energy = 0;
    for (int r = 0; r < m; ++r) {
        prefix_work += work_values[r];
        prefix_energy += energies[r];
        double completion_lb = t0 + prefix_work
            + t.alpha * std::ceil(static_cast<double>(prefix_energy) / t.capacity);
        if (completion_lb > t.due_dates[sc + r] + EPS) forced[r] = 1;
    }

    int q_min = bin_packing_lower_bound(t, rem);
    double makespan_lb = t0 + rem_p + t.beta * rem_e + t.alpha * q_min;
    std::sort(work_values.rbegin(), work_values.rend());
    double tail_deduction = 0.0;
    for (int tail_count = 0; tail_count < m; ++tail_count) {
        int pos = m - 1 - tail_count;
        double completion_lb = makespan_lb - tail_deduction - t.alpha * tail_count;
        if (completion_lb > t.due_dates[sc + pos] + EPS) forced[pos] = 1;
        tail_deduction += work_values[tail_count];
    }

    return static_cast<int>(std::count(forced.begin(), forced.end(), 1));
}

int node_lower_bound(
    const TypedInstance& t, const std::vector<int>& rem, int sc, double t0,
    bool use_position_dp = false
) {
    int lb = std::max({lb_fractional(t, rem, sc, t0),
                       lb_single_setup(t, rem, sc, t0),
                       lb_tail(t, rem, sc, t0),
                       lb_position_count(t, rem, sc, t0),
                       lb_forced_position_union(t, rem, sc, t0)});
    if (use_position_dp) {
        if (sc == 0 && std::abs(t0) <= EPS) {
            lb = std::max(
                lb,
                lb_position_dp(
                    t, rem, sc, t0,
                    env_int("ROOT_POSITION_DP_MAX_JOBS", 220),
                    env_int("ROOT_POSITION_DP_MAX_ENERGY", 12000)
                )
            );
        } else {
            lb = std::max(lb, lb_position_dp(t, rem, sc, t0));
        }
    }
    return lb;
}

SolveResult solve_constructive_heuristic(const Instance& inst, double time_limit) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto blocks = rigorous_start_blocks_with_indices(inst);
    ScheduleEval eval = evaluate_blocks(inst, blocks);
    TypedInstance typed = to_typed_instance(inst);
    int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);
    double runtime = std::chrono::duration<double>(Clock::now() - start).count();
    SolveResult res{
        "Relaxation-guided DP start", eval.late, static_cast<double>(root_lb),
        eval.late <= root_lb, 0, std::min(runtime, time_limit),
        static_cast<int>(inst.jobs.size()), time_limit, true
    };
    attach_schedule_metrics(res, inst, blocks);
    return res;
}

SolveResult solve_constructive_heuristic_no_resegmentation(const Instance& inst, double time_limit) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto blocks = block_packing_improvement_heuristic(inst, false);
    ScheduleEval eval = evaluate_blocks(inst, blocks);
    TypedInstance typed = to_typed_instance(inst);
    int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);
    double runtime = std::chrono::duration<double>(Clock::now() - start).count();
    SolveResult res{
        "H0 without sequence DP", eval.late, static_cast<double>(root_lb),
        eval.late <= root_lb, 0, std::min(runtime, time_limit),
        static_cast<int>(inst.jobs.size()), time_limit, true
    };
    attach_schedule_metrics(res, inst, blocks);
    return res;
}

struct WindowPolishResult {
    bool feasible = false;
    std::vector<std::vector<int>> blocks;
    long long nodes = 0;
};

#ifdef WITH_CPLEX
SolveResult solve_cplex_position(
    const Instance& inst,
    double time_limit,
    bool use_heuristic_start,
    bool use_strengthening,
    const std::string& algorithm_name
);

SolveResult solve_cplex_full_recharge_position(
    const Instance& inst,
    double time_limit,
    bool use_heuristic_start,
    bool use_strengthening,
    const std::string& algorithm_name
);

WindowPolishResult polish_fixed_window_cplex(
    const Instance& inst,
    const std::vector<std::vector<int>>& current_blocks,
    const std::vector<char>& free_position,
    double time_limit,
    bool use_strengthening,
    int tardy_lb,
    int recharge_lb
) {
    const int n = static_cast<int>(inst.jobs.size());
    const int E = inst.capacity;
    std::vector<int> current_sequence = flatten_blocks(current_blocks);
    std::vector<int> free_positions;
    free_positions.reserve(n);
    for (int pos = 0; pos < n; ++pos) {
        if (free_position[pos]) free_positions.push_back(pos);
    }
    if (free_positions.empty()) return {};
    const int m = static_cast<int>(free_positions.size());
    std::vector<int> free_jobs;
    free_jobs.reserve(m);
    std::unordered_map<int, int> free_job_to_local;
    for (int local = 0; local < m; ++local) {
        int job = current_sequence[free_positions[local]];
        free_jobs.push_back(job);
        free_job_to_local[job] = local;
    }

    int sum_p = 0;
    for (const auto& job : inst.jobs) sum_p += job.first;
    const double horizon_ub = sum_p + n * (inst.alpha + inst.beta * E);

    IloEnv env;
    try {
        IloModel model(env);
        std::vector<IloBoolVarArray> x;
        x.reserve(m);
        for (int j = 0; j < m; ++j) {
            x.emplace_back(env, m);
            for (int i = 0; i < m; ++i) {
                std::ostringstream name;
                name << "xw_" << free_jobs[j] << "_" << free_positions[i];
                x[j][i].setName(name.str().c_str());
            }
        }
        IloBoolVarArray z(env, n);
        IloBoolVarArray u(env, n);
        IloNumVarArray r(env, n, 0.0, E);
        IloNumVarArray g(env, n, 0.0, E);
        IloNumVarArray ctime(env, n, 0.0, horizon_ub);

        for (int j = 0; j < m; ++j) {
            IloExpr expr(env);
            for (int i = 0; i < m; ++i) expr += x[j][i];
            model.add(expr == 1);
            expr.end();
        }
        for (int i = 0; i < m; ++i) {
            IloExpr expr(env);
            for (int j = 0; j < m; ++j) expr += x[j][i];
            model.add(expr == 1);
            expr.end();
        }

        for (int i = 0; i < n; ++i) {
            IloExpr pexpr(env), eexpr(env);
            int free_local_pos = -1;
            if (free_position[i]) {
                auto it = std::lower_bound(free_positions.begin(), free_positions.end(), i);
                free_local_pos = static_cast<int>(it - free_positions.begin());
            }
            if (free_local_pos >= 0) {
                for (int j = 0; j < m; ++j) {
                    int job = free_jobs[j];
                    pexpr += inst.jobs[job].first * x[j][free_local_pos];
                    eexpr += inst.jobs[job].second * x[j][free_local_pos];
                }
            } else {
                int job = current_sequence[i];
                pexpr += inst.jobs[job].first;
                eexpr += inst.jobs[job].second;
            }

            if (i == 0) {
                model.add(ctime[i] == inst.alpha * z[i] + inst.beta * r[i] + pexpr);
                model.add(g[i] == r[i] - eexpr);
                model.add(r[i] <= E * z[i]);
            } else {
                model.add(ctime[i] == ctime[i - 1] + inst.alpha * z[i] + inst.beta * r[i] + pexpr);
                model.add(g[i] == g[i - 1] + r[i] - eexpr);
                model.add(g[i - 1] + r[i] <= E);
                model.add(r[i] <= E * z[i]);
            }
            double m_i = std::max(0.0, horizon_ub - inst.due_dates[i]);
            model.add(ctime[i] <= inst.due_dates[i] + m_i * u[i]);
            pexpr.end();
            eexpr.end();
        }

        IloExpr tardy_count(env);
        for (int i = 0; i < n; ++i) tardy_count += u[i];
        if (use_strengthening) {
            model.add(tardy_count >= tardy_lb);

            IloExpr recharge_count(env);
            for (int i = 0; i < n; ++i) recharge_count += z[i];
            model.add(recharge_count >= recharge_lb);
            recharge_count.end();

            model.add(g[n - 1] == 0);
        }

        IloExpr obj(env);
        obj += 100000.0 * tardy_count;
        obj += ctime[n - 1];
        model.add(IloMinimize(env, obj));
        obj.end();
        tardy_count.end();

        IloCplex cplex(model);
        IloNumVarArray start_vars(env);
        IloNumArray start_vals(env);
        std::vector<double> z_val(n, 0.0), u_val(n, 0.0), r_val(n, 0.0), g_val(n, 0.0), c_val(n, 0.0);
        std::vector<std::vector<double>> x_val(m, std::vector<double>(m, 0.0));
        for (int local_pos = 0; local_pos < m; ++local_pos) {
            int job = current_sequence[free_positions[local_pos]];
            auto it = free_job_to_local.find(job);
            if (it != free_job_to_local.end()) x_val[it->second][local_pos] = 1.0;
        }
        double time = 0.0;
        int pos = 0;
        for (const auto& block : current_blocks) {
            int energy = block_energy(inst, block);
            int residual = energy;
            time += inst.alpha + inst.beta * energy;
            for (int local = 0; local < static_cast<int>(block.size()); ++local) {
                int job = block[local];
                if (local == 0) {
                    z_val[pos] = 1.0;
                    r_val[pos] = energy;
                }
                residual -= inst.jobs[job].second;
                g_val[pos] = residual;
                time += inst.jobs[job].first;
                c_val[pos] = time;
                u_val[pos] = time > inst.due_dates[pos] + EPS ? 1.0 : 0.0;
                ++pos;
            }
        }
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < m; ++i) {
                start_vars.add(x[j][i]);
                start_vals.add(x_val[j][i]);
            }
        }
        for (int i = 0; i < n; ++i) {
            start_vars.add(z[i]); start_vals.add(z_val[i]);
            start_vars.add(u[i]); start_vals.add(u_val[i]);
            start_vars.add(r[i]); start_vals.add(r_val[i]);
            start_vars.add(g[i]); start_vals.add(g_val[i]);
            start_vars.add(ctime[i]); start_vals.add(c_val[i]);
        }
        cplex.addMIPStart(start_vars, start_vals, IloCplex::MIPStartAuto, "current_schedule");
        start_vars.end();
        start_vals.end();

        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::TiLim, std::max(0.05, time_limit));
        cplex.setParam(IloCplex::Threads, 1);
        bool ok = cplex.solve();
        WindowPolishResult out;
        out.nodes = static_cast<long long>(cplex.getNnodes());
        if (ok) {
            std::vector<int> sequence(n, -1);
            for (int i = 0; i < n; ++i) sequence[i] = current_sequence[i];
            for (int local_pos = 0; local_pos < m; ++local_pos) {
                int pos = free_positions[local_pos];
                for (int j = 0; j < m; ++j) {
                    if (cplex.getValue(x[j][local_pos]) > 0.5) {
                        sequence[pos] = free_jobs[j];
                        break;
                    }
                }
            }
            out.blocks = best_blocks_for_fixed_sequence(inst, sequence);
            out.feasible = true;
        }
        env.end();
        return out;
    } catch (const IloException& ex) {
        std::string msg = ex.getMessage();
        env.end();
        throw std::runtime_error("CPLEX exception in window polish: " + msg);
    }
}

WindowPolishResult polish_fixed_window_cplex(
    const Instance& inst,
    const std::vector<std::vector<int>>& current_blocks,
    int left,
    int right,
    double time_limit,
    bool use_strengthening,
    int tardy_lb,
    int recharge_lb
) {
    std::vector<char> free_position(inst.jobs.size(), 0);
    for (int pos = left; pos <= right; ++pos) free_position[pos] = 1;
    return polish_fixed_window_cplex(
        inst, current_blocks, free_position, time_limit,
        use_strengthening, tardy_lb, recharge_lb
    );
}

WindowPolishResult polish_local_branching_cplex(
    const Instance& inst,
    const std::vector<std::vector<int>>& current_blocks,
    int radius,
    double time_limit,
    bool use_strengthening,
    int tardy_lb,
    int recharge_lb
) {
    const int n = static_cast<int>(inst.jobs.size());
    const int E = inst.capacity;
    std::vector<int> current_sequence = flatten_blocks(current_blocks);
    ScheduleEval current_eval = evaluate_blocks(inst, current_blocks);
    radius = std::max(1, std::min(radius, n));

    int sum_p = 0;
    for (const auto& job : inst.jobs) sum_p += job.first;
    const double horizon_ub = sum_p + n * (inst.alpha + inst.beta * E);

    IloEnv env;
    try {
        IloModel model(env);
        std::vector<IloBoolVarArray> x;
        x.reserve(n);
        for (int j = 0; j < n; ++j) {
            x.emplace_back(env, n);
            for (int i = 0; i < n; ++i) {
                std::ostringstream name;
                name << "xlb_" << j << "_" << i;
                x[j][i].setName(name.str().c_str());
            }
        }
        IloBoolVarArray z(env, n);
        IloBoolVarArray u(env, n);
        IloNumVarArray r(env, n, 0.0, E);
        IloNumVarArray g(env, n, 0.0, E);
        IloNumVarArray ctime(env, n, 0.0, horizon_ub);

        for (int j = 0; j < n; ++j) {
            IloExpr expr(env);
            for (int i = 0; i < n; ++i) expr += x[j][i];
            model.add(expr == 1);
            expr.end();
        }
        for (int i = 0; i < n; ++i) {
            IloExpr expr(env);
            for (int j = 0; j < n; ++j) expr += x[j][i];
            model.add(expr == 1);
            expr.end();
        }

        for (int i = 0; i < n; ++i) {
            IloExpr pexpr(env), eexpr(env);
            for (int j = 0; j < n; ++j) {
                pexpr += inst.jobs[j].first * x[j][i];
                eexpr += inst.jobs[j].second * x[j][i];
            }
            if (i == 0) {
                model.add(ctime[i] == inst.alpha * z[i] + inst.beta * r[i] + pexpr);
                model.add(g[i] == r[i] - eexpr);
                model.add(r[i] <= E * z[i]);
            } else {
                model.add(ctime[i] == ctime[i - 1] + inst.alpha * z[i] + inst.beta * r[i] + pexpr);
                model.add(g[i] == g[i - 1] + r[i] - eexpr);
                model.add(g[i - 1] + r[i] <= E);
                model.add(r[i] <= E * z[i]);
            }
            double m_i = std::max(0.0, horizon_ub - inst.due_dates[i]);
            model.add(ctime[i] <= inst.due_dates[i] + m_i * u[i]);
            pexpr.end();
            eexpr.end();
        }

        IloExpr tardy_count(env);
        for (int i = 0; i < n; ++i) tardy_count += u[i];
        model.add(tardy_count <= current_eval.late);
        if (use_strengthening) {
            model.add(tardy_count >= tardy_lb);
            IloExpr recharge_count(env);
            for (int i = 0; i < n; ++i) recharge_count += z[i];
            model.add(recharge_count >= recharge_lb);
            recharge_count.end();
            model.add(g[n - 1] == 0);
        }

        IloExpr keep(env);
        for (int i = 0; i < n; ++i) keep += x[current_sequence[i]][i];
        model.add(keep >= n - radius);
        keep.end();

        IloExpr obj(env);
        obj += 100000.0 * tardy_count;
        obj += ctime[n - 1];
        model.add(IloMinimize(env, obj));
        obj.end();
        tardy_count.end();

        IloNumVarArray start_vars(env);
        IloNumArray start_vals(env);
        std::vector<double> z_val(n, 0.0), u_val(n, 0.0), r_val(n, 0.0), g_val(n, 0.0), c_val(n, 0.0);
        std::vector<std::vector<double>> x_val(n, std::vector<double>(n, 0.0));
        double time = 0.0;
        int pos = 0;
        for (const auto& block : current_blocks) {
            int energy = block_energy(inst, block);
            int residual = energy;
            time += inst.alpha + inst.beta * energy;
            for (int local = 0; local < static_cast<int>(block.size()); ++local) {
                int job = block[local];
                x_val[job][pos] = 1.0;
                if (local == 0) {
                    z_val[pos] = 1.0;
                    r_val[pos] = energy;
                }
                residual -= inst.jobs[job].second;
                g_val[pos] = residual;
                time += inst.jobs[job].first;
                c_val[pos] = time;
                u_val[pos] = time > inst.due_dates[pos] + EPS ? 1.0 : 0.0;
                ++pos;
            }
        }
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                start_vars.add(x[j][i]);
                start_vals.add(x_val[j][i]);
            }
        }
        for (int i = 0; i < n; ++i) {
            start_vars.add(z[i]); start_vals.add(z_val[i]);
            start_vars.add(u[i]); start_vals.add(u_val[i]);
            start_vars.add(r[i]); start_vals.add(r_val[i]);
            start_vars.add(g[i]); start_vals.add(g_val[i]);
            start_vars.add(ctime[i]); start_vals.add(c_val[i]);
        }

        IloCplex cplex(model);
        cplex.addMIPStart(start_vars, start_vals, IloCplex::MIPStartAuto, "current_schedule");
        start_vars.end();
        start_vals.end();
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::TiLim, std::max(0.05, time_limit));
        cplex.setParam(IloCplex::Threads, 1);
        bool ok = cplex.solve();
        WindowPolishResult out;
        out.nodes = static_cast<long long>(cplex.getNnodes());
        if (ok) {
            std::vector<int> sequence(n, -1);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (cplex.getValue(x[j][i]) > 0.5) {
                        sequence[i] = j;
                        break;
                    }
                }
            }
            if (std::find(sequence.begin(), sequence.end(), -1) == sequence.end()) {
                out.blocks = best_blocks_for_fixed_sequence(inst, sequence);
                out.feasible = true;
            }
        }
        env.end();
        return out;
    } catch (const IloException& ex) {
        std::string msg = ex.getMessage();
        env.end();
        throw std::runtime_error("CPLEX exception in local branching polish: " + msg);
    }
}

std::vector<std::pair<int, int>> critical_windows(
    const Instance& inst,
    const std::vector<std::vector<int>>& blocks,
    int round
) {
    const int n = static_cast<int>(inst.jobs.size());
    std::vector<double> c = completion_times(inst, blocks);
    std::vector<int> centers;

    int first_late = -1;
    int max_tardy = 0;
    double max_tardiness = -1e100;
    int tightest = 0;
    double min_slack = 1e100;
    for (int i = 0; i < n; ++i) {
        double slack = inst.due_dates[i] - c[i];
        if (first_late < 0 && slack < -EPS) first_late = i;
        if (-slack > max_tardiness) {
            max_tardiness = -slack;
            max_tardy = i;
        }
        if (slack < min_slack) {
            min_slack = slack;
            tightest = i;
        }
    }
    if (first_late >= 0) {
        centers.push_back(first_late);
        centers.push_back(std::max(0, first_late - 4));
        centers.push_back(std::min(n - 1, first_late + 4));
    }
    centers.push_back(tightest);
    centers.push_back(max_tardy);
    centers.push_back(n / 3);
    centers.push_back(n / 2);
    centers.push_back(2 * n / 3);

    std::vector<int> sizes;
    if (n <= 70) sizes = {10, 14, 18, 24, 30};
    else if (n <= 140) sizes = {12, 18, 24, 32};
    else sizes = {16, 24, 32, 44};
    std::rotate(sizes.begin(), sizes.begin() + (round % static_cast<int>(sizes.size())), sizes.end());

    std::set<std::pair<int, int>> seen;
    std::vector<std::pair<int, int>> windows;
    for (int center : centers) {
        for (int size : sizes) {
            int left = std::max(0, center - size / 2);
            int right = std::min(n - 1, left + size - 1);
            left = std::max(0, right - size + 1);
            if (right - left + 1 < 2) continue;
            if (seen.insert({left, right}).second) windows.push_back({left, right});
        }
    }
    return windows;
}

std::vector<char> alns_destroy_positions(
    const Instance& inst,
    const std::vector<std::vector<int>>& blocks,
    int q,
    int op,
    std::mt19937_64& rng
) {
    const int n = static_cast<int>(inst.jobs.size());
    std::vector<int> seq = flatten_blocks(blocks);
    std::vector<char> free_position(n, 0);
    q = std::max(2, std::min(q, n));

    auto mark_ranked = [&](std::vector<int> ranked) {
        int marked = 0;
        for (int pos : ranked) {
            if (pos < 0 || pos >= n || free_position[pos]) continue;
            free_position[pos] = 1;
            if (++marked >= q) break;
        }
        if (marked < q) {
            std::vector<int> all(n);
            std::iota(all.begin(), all.end(), 0);
            std::shuffle(all.begin(), all.end(), rng);
            for (int pos : all) {
                if (free_position[pos]) continue;
                free_position[pos] = 1;
                if (++marked >= q) break;
            }
        }
    };

    if (op == 0) {
        auto c = completion_times(inst, blocks);
        std::vector<int> lb_order = effective_load_order(inst);
        const double lambda = (inst.alpha + inst.beta * inst.capacity) / inst.capacity;
        double relaxed_time = 0.0;
        int critical = n / 2;
        for (int i = 0; i < n; ++i) {
            int job = lb_order[i];
            relaxed_time += inst.jobs[job].first + lambda * inst.jobs[job].second;
            if (relaxed_time > inst.due_dates[i] + EPS) {
                critical = i;
                break;
            }
        }
        std::vector<int> ranked(n);
        std::iota(ranked.begin(), ranked.end(), 0);
        std::stable_sort(ranked.begin(), ranked.end(), [&](int a, int b) {
            double slack_a = inst.due_dates[a] - c[a];
            double slack_b = inst.due_dates[b] - c[b];
            double ka = std::abs(a - critical) + 0.05 * slack_a;
            double kb = std::abs(b - critical) + 0.05 * slack_b;
            if (std::abs(ka - kb) > EPS) return ka < kb;
            return slack_a < slack_b;
        });
        mark_ranked(std::move(ranked));
    } else if (op == 1) {
        std::uniform_int_distribution<int> seed_dist(0, n - 1);
        int seed_pos = seed_dist(rng);
        int seed_job = seq[seed_pos];
        double p_scale = 1.0;
        double e_scale = 1.0;
        for (const auto& job : inst.jobs) {
            p_scale = std::max(p_scale, static_cast<double>(job.first));
            e_scale = std::max(e_scale, static_cast<double>(job.second));
        }
        std::vector<int> ranked(n);
        std::iota(ranked.begin(), ranked.end(), 0);
        std::stable_sort(ranked.begin(), ranked.end(), [&](int a, int b) {
            int ja = seq[a];
            int jb = seq[b];
            double da = std::abs(inst.jobs[ja].first - inst.jobs[seed_job].first) / p_scale
                + std::abs(inst.jobs[ja].second - inst.jobs[seed_job].second) / e_scale
                + 0.02 * std::abs(a - seed_pos);
            double db = std::abs(inst.jobs[jb].first - inst.jobs[seed_job].first) / p_scale
                + std::abs(inst.jobs[jb].second - inst.jobs[seed_job].second) / e_scale
                + 0.02 * std::abs(b - seed_pos);
            if (std::abs(da - db) > EPS) return da < db;
            return a < b;
        });
        mark_ranked(std::move(ranked));
    } else {
        std::vector<int> ranked(n);
        std::iota(ranked.begin(), ranked.end(), 0);
        std::shuffle(ranked.begin(), ranked.end(), rng);
        mark_ranked(std::move(ranked));
    }
    return free_position;
}

SolveResult solve_grasp_lns_matheuristic(const Instance& inst, double time_limit) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };

    const int n = static_cast<int>(inst.jobs.size());
    std::uint64_t seed = static_cast<std::uint64_t>(
        instance_seed(inst.config) ^ 0x9e3779b97f4a7c15ULL
    );
    auto current_blocks = grasp_effective_load_construction(
        inst, std::min(2.0, 0.15 * time_limit), seed
    );
    auto best_blocks = current_blocks;
    ScheduleEval current_eval = evaluate_blocks(inst, current_blocks);
    ScheduleEval best_eval = current_eval;
    const int initial = best_eval.late;

    TypedInstance typed = to_typed_instance(inst);
    int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);
    int recharge_lb = bin_packing_lower_bound(typed, typed.count);
    if (best_eval.late <= root_lb) {
        return {
            "GRASP-ALNS exact-repair matheuristic", best_eval.late,
            static_cast<double>(root_lb), true, 0, elapsed(),
            initial, time_limit, true
        };
    }

    double rebuild_budget = std::min(time_limit - elapsed(), std::max(0.0, 0.35 * time_limit));
    if (rebuild_budget > 0.05) {
        target_set_grasp_rebuild(
            inst, best_blocks, best_eval, root_lb, rebuild_budget,
            seed ^ 0x94d049bb133111ebULL
        );
        current_blocks = best_blocks;
        current_eval = best_eval;
    }

    std::mt19937_64 rng(seed ^ 0xd1b54a32d192ed03ULL);
    std::vector<double> weights = {1.0, 1.0, 1.0};
    const double reaction = 0.20;
    int q_min = std::max(4, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n)))));
    int q_max = std::max(q_min, static_cast<int>(std::ceil(0.25 * n)));
    int q = std::max(q_min, static_cast<int>(std::ceil(0.12 * n)));
    int no_improve = 0;
    long long neighborhoods = 0;
    long long nodes = 0;

    while (elapsed() < time_limit - 0.10 && best_eval.late > root_lb) {
        std::discrete_distribution<int> choose_op(weights.begin(), weights.end());
        int op = choose_op(rng);
        std::vector<char> free_position = alns_destroy_positions(inst, current_blocks, q, op, rng);
        double remaining = time_limit - elapsed();
        double sub_time = std::min(remaining, std::max(0.25, 0.06 * time_limit));
        auto cand = polish_fixed_window_cplex(
            inst, current_blocks, free_position, sub_time, true, root_lb, recharge_lb
        );
        ++neighborhoods;
        nodes += cand.nodes;
        double reward = 0.05;
        if (cand.feasible) {
            ScheduleEval cand_eval = evaluate_blocks(inst, cand.blocks);
            if (better_schedule(cand_eval, current_eval)) {
                current_blocks = std::move(cand.blocks);
                current_eval = cand_eval;
                reward = 1.0;
                if (better_schedule(current_eval, best_eval)) {
                    best_blocks = current_blocks;
                    best_eval = current_eval;
                    reward = 5.0;
                    no_improve = 0;
                    q = std::max(q_min, static_cast<int>(std::floor(0.85 * q)));
                }
            } else {
                ++no_improve;
            }
        } else {
            ++no_improve;
        }
        weights[op] = (1.0 - reaction) * weights[op] + reaction * reward;
        if (no_improve >= 3) {
            q = std::min(q_max, std::max(q + 1, static_cast<int>(std::ceil(1.20 * q))));
            no_improve = 0;
        }
    }

    return {
        "GRASP-ALNS exact-repair matheuristic", best_eval.late,
        static_cast<double>(root_lb), best_eval.late <= root_lb,
        std::max(neighborhoods, nodes), elapsed(), initial, time_limit, true
    };
}

SolveResult solve_grasp_local_branching_matheuristic(const Instance& inst, double time_limit) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };

    const int n = static_cast<int>(inst.jobs.size());
    std::uint64_t seed = static_cast<std::uint64_t>(
        instance_seed(inst.config) ^ 0x517cc1b727220a95ULL
    );
    auto best_blocks = grasp_effective_load_construction(
        inst, std::min(2.0, 0.15 * time_limit), seed
    );
    ScheduleEval best_eval = evaluate_blocks(inst, best_blocks);
    const int initial = best_eval.late;

    TypedInstance typed = to_typed_instance(inst);
    int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);
    int recharge_lb = bin_packing_lower_bound(typed, typed.count);
    long long neighborhoods = 0;
    long long nodes = 0;

    if (best_eval.late > root_lb) {
        double rebuild_budget = std::min(time_limit - elapsed(), std::max(0.0, 0.40 * time_limit));
        if (rebuild_budget > 0.05) {
            target_set_grasp_rebuild(
                inst, best_blocks, best_eval, root_lb, rebuild_budget,
                seed ^ 0x94d049bb133111ebULL
            );
        }
    }

    std::vector<double> radius_fraction = {0.04, 0.06, 0.08, 0.12, 0.16};
    int sweep = 0;
    while (elapsed() < time_limit - 0.10 && best_eval.late > root_lb) {
        bool improved = false;
        for (double frac : radius_fraction) {
            double remaining = time_limit - elapsed();
            if (remaining <= 0.10 || best_eval.late <= root_lb) break;
            int radius = std::max(2, static_cast<int>(std::ceil(frac * n)));
            if (sweep > 0) radius = std::min(n, static_cast<int>(std::ceil(radius * (1.0 + 0.25 * sweep))));
            double sub_time = std::min(remaining, std::max(0.5, 0.12 * time_limit));
            auto cand = polish_local_branching_cplex(
                inst, best_blocks, radius, sub_time, true, root_lb, recharge_lb
            );
            ++neighborhoods;
            nodes += cand.nodes;
            if (!cand.feasible) continue;
            ScheduleEval cand_eval = evaluate_blocks(inst, cand.blocks);
            if (better_schedule(cand_eval, best_eval)) {
                best_blocks = std::move(cand.blocks);
                best_eval = cand_eval;
                improved = true;
                break;
            }
        }
        if (!improved) ++sweep;
        if (sweep >= 3) break;
    }

    return {
        "GRASP local-branching matheuristic", best_eval.late,
        static_cast<double>(root_lb), best_eval.late <= root_lb,
        std::max(neighborhoods, nodes), elapsed(), initial, time_limit, true
    };
}

SolveResult solve_grasp_dp_matheuristic(const Instance& inst, double time_limit) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };

    std::uint64_t seed = static_cast<std::uint64_t>(
        instance_seed(inst.config) ^ 0xbf58476d1ce4e5b9ULL
    );
    auto best_blocks = grasp_effective_load_construction(
        inst, std::min(time_limit, std::max(1.0, 0.25 * time_limit)), seed, 0.20
    );
    ScheduleEval best_eval = evaluate_blocks(inst, best_blocks);
    const int initial = best_eval.late;

    TypedInstance typed = to_typed_instance(inst);
    int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);
    long long neighborhoods = 0;
    int no_improve = 0;
    while (elapsed() < time_limit - EPS && best_eval.late > root_lb) {
        double slice = std::min(time_limit - elapsed(), std::max(0.25, 0.15 * time_limit));
        ScheduleEval before = best_eval;
        bool improved = target_set_grasp_rebuild(
            inst, best_blocks, best_eval, root_lb, slice,
            seed ^ static_cast<std::uint64_t>(0x9e3779b97f4a7c15ULL + neighborhoods * 1315423911ULL)
        );
        ++neighborhoods;
        if (improved && better_schedule(best_eval, before)) {
            no_improve = 0;
        } else if (++no_improve >= 50) {
            break;
        }
    }

    return {
        "GRASP-DP matheuristic", best_eval.late,
        static_cast<double>(root_lb), best_eval.late <= root_lb,
        neighborhoods, elapsed(), initial, time_limit, true
    };
}

SolveResult solve_grasp_pr_matheuristic(
    const Instance& inst,
    double time_limit,
    bool long_search = false,
    bool full_recharge = false
) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };

    std::uint64_t seed = static_cast<std::uint64_t>(
        instance_seed(inst.config) ^ 0x94d049bb133111ebULL
    );
    std::mt19937_64 rng(seed);
    TypedInstance typed = to_typed_instance(inst);
    int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);

    auto base = evaluate_sequence_solution(inst, effective_load_order(inst), full_recharge);
    const int initial = base.eval.late;
    long long evaluations = 0;
    double initial_ls_budget = std::min(time_limit - elapsed(), std::max(0.25, 0.08 * time_limit));
    evaluations += bound_guided_local_search(inst, base, root_lb, initial_ls_budget, full_recharge);

    SequenceSolution best = base;
    std::vector<SequenceSolution> elite;
    const int elite_size = env_int("GRASP_ELITE_SIZE", 10);
    update_elite_pool(elite, base, elite_size);

    int no_improve = 0;
    while (elapsed() < time_limit - EPS && best.eval.late > root_lb) {
        double remaining = time_limit - elapsed();
        if (remaining <= EPS) break;
        bool improved_iteration = false;
        const double epsilon = env_double("GRASP_EPS", 0.20);
        auto seq = randomized_effective_load_sequence(inst, rng, epsilon);
        auto cand = evaluate_sequence_solution(inst, seq, full_recharge);
        ++evaluations;

        double ls_budget = std::min(
            time_limit - elapsed(),
            std::max(0.15, env_double("GRASP_LS_BUDGET_FRAC", 0.05) * time_limit)
        );
        evaluations += bound_guided_local_search(inst, cand, root_lb, ls_budget, full_recharge);

        if (better_schedule(cand.eval, best.eval)) {
            best = cand;
            improved_iteration = true;
        }
        update_elite_pool(elite, cand, elite_size);

        if (elite.size() >= 2 && elapsed() < time_limit - EPS) {
            std::uniform_int_distribution<int> pick(0, static_cast<int>(elite.size()) - 1);
            int idx = pick(rng);
            if (same_sequence(elite[idx].sequence, cand.sequence)) {
                idx = (idx + 1) % static_cast<int>(elite.size());
            }
            SequenceSolution path_best;
            double pr_budget = std::min(
                time_limit - elapsed(),
                std::max(0.15, env_double("GRASP_PR_BUDGET_FRAC", 0.04) * time_limit)
            );
            evaluations += path_relink(inst, cand, elite[idx], path_best, pr_budget, full_recharge);
            double path_ls_budget = std::min(
                time_limit - elapsed(),
                std::max(0.10, env_double("GRASP_PATH_LS_BUDGET_FRAC", 0.03) * time_limit)
            );
            evaluations += bound_guided_local_search(inst, path_best, root_lb, path_ls_budget, full_recharge);
            if (better_schedule(path_best.eval, best.eval)) {
                best = path_best;
                improved_iteration = true;
            }
            update_elite_pool(elite, std::move(path_best), elite_size);
        }

        if (improved_iteration) {
            no_improve = 0;
        } else {
            ++no_improve;
            int stop_no_improve = env_int("GRASP_STOP_NO_IMPROVE", 25);
            if (no_improve >= stop_no_improve) break;
        }
    }

    SolveResult res{
        full_recharge ? (long_search ? "GRASP-PR full-recharge heuristic long"
                                      : "GRASP-PR full-recharge heuristic")
                      : (long_search ? "GRASP-PR heuristic long"
                                     : "GRASP-PR heuristic"),
        best.eval.late, static_cast<double>(root_lb), best.eval.late <= root_lb,
        evaluations, elapsed(), initial, time_limit, true
    };
    attach_schedule_metrics(res, inst, best.blocks, full_recharge);
    return res;
}

bool global_late_set_rebuild_search(
    const Instance& inst,
    std::vector<std::vector<int>>& best_blocks,
    ScheduleEval& best_eval,
    int root_lb,
    const std::chrono::steady_clock::time_point& start,
    double time_limit,
    int round
) {
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    };
    const int n = static_cast<int>(inst.jobs.size());
    if (n < 100 || best_eval.late <= root_lb) return false;
    const double fac = (inst.alpha + inst.beta * inst.capacity) / inst.capacity;

    std::vector<int> current_sequence = flatten_blocks(best_blocks);
    std::vector<int> current_pos(n, 0);
    for (int pos = 0; pos < n; ++pos) current_pos[current_sequence[pos]] = pos;

    auto contribution = [&](int job) {
        return inst.jobs[job].first + fac * inst.jobs[job].second;
    };
    auto work = [&](int job) {
        return inst.jobs[job].first + inst.beta * inst.jobs[job].second;
    };

    std::vector<int> jobs(n);
    std::iota(jobs.begin(), jobs.end(), 0);

    std::vector<std::vector<int>> late_rankings;
    auto add_late_ranking = [&](auto cmp) {
        std::vector<int> order = jobs;
        std::stable_sort(order.begin(), order.end(), cmp);
        late_rankings.push_back(std::move(order));
    };
    add_late_ranking([&](int a, int b) {
        double ca = contribution(a), cb = contribution(b);
        if (std::abs(ca - cb) > EPS) return ca > cb;
        return current_pos[a] < current_pos[b];
    });
    add_late_ranking([&](int a, int b) {
        if (inst.jobs[a].first != inst.jobs[b].first) return inst.jobs[a].first > inst.jobs[b].first;
        return inst.jobs[a].second > inst.jobs[b].second;
    });
    add_late_ranking([&](int a, int b) {
        if (inst.jobs[a].second != inst.jobs[b].second) return inst.jobs[a].second > inst.jobs[b].second;
        return inst.jobs[a].first > inst.jobs[b].first;
    });
    add_late_ranking([&](int a, int b) {
        double ca = work(a), cb = work(b);
        if (std::abs(ca - cb) > EPS) return ca > cb;
        return inst.jobs[a].second > inst.jobs[b].second;
    });
    add_late_ranking([&](int a, int b) {
        return current_pos[a] > current_pos[b];
    });

    enum class EarlyOrder { Current, SPT, Contribution, Work, Energy };
    std::vector<EarlyOrder> early_orders{
        EarlyOrder::Current, EarlyOrder::SPT, EarlyOrder::Contribution,
        EarlyOrder::Work, EarlyOrder::Energy
    };

    auto sort_early = [&](std::vector<int>& early, EarlyOrder mode) {
        std::stable_sort(early.begin(), early.end(), [&](int a, int b) {
            if (mode == EarlyOrder::Current) return current_pos[a] < current_pos[b];
            if (mode == EarlyOrder::SPT) {
                if (inst.jobs[a].first != inst.jobs[b].first) return inst.jobs[a].first < inst.jobs[b].first;
                return inst.jobs[a].second < inst.jobs[b].second;
            }
            if (mode == EarlyOrder::Contribution) {
                double ca = contribution(a), cb = contribution(b);
                if (std::abs(ca - cb) > EPS) return ca < cb;
                return inst.jobs[a].first < inst.jobs[b].first;
            }
            if (mode == EarlyOrder::Work) {
                double ca = work(a), cb = work(b);
                if (std::abs(ca - cb) > EPS) return ca < cb;
                return inst.jobs[a].first < inst.jobs[b].first;
            }
            if (inst.jobs[a].second != inst.jobs[b].second) return inst.jobs[a].second < inst.jobs[b].second;
            return inst.jobs[a].first < inst.jobs[b].first;
        });
    };

    auto sort_late_tail = [&](std::vector<int>& late) {
        std::stable_sort(late.begin(), late.end(), [&](int a, int b) {
            double ca = contribution(a), cb = contribution(b);
            if (std::abs(ca - cb) > EPS) return ca > cb;
            return current_pos[a] > current_pos[b];
        });
    };

    auto evaluate_late_set = [&](const std::vector<int>& late_jobs, EarlyOrder mode) {
        std::vector<char> is_late(n, 0);
        for (int job : late_jobs) is_late[job] = 1;
        std::vector<int> early;
        early.reserve(n - static_cast<int>(late_jobs.size()));
        for (int job : jobs) if (!is_late[job]) early.push_back(job);
        std::vector<int> late = late_jobs;
        sort_early(early, mode);
        sort_late_tail(late);
        early.insert(early.end(), late.begin(), late.end());
        auto cand_blocks = best_blocks_for_fixed_sequence(inst, early);
        ScheduleEval cand_eval = evaluate_blocks(inst, cand_blocks);
        if (better_schedule(cand_eval, best_eval)) {
            best_blocks = std::move(cand_blocks);
            best_eval = cand_eval;
            trace_bgmh_event(inst, "global", round, elapsed(), best_eval.late);
            return true;
        }
        return false;
    };

    bool improved = false;
    int max_depth = std::min(best_eval.late - root_lb, 5);
    for (int depth = 1; depth <= max_depth; ++depth) {
        int target_late = best_eval.late - depth;
        if (target_late < root_lb || target_late <= 0) break;
        for (const auto& ranking : late_rankings) {
            if (elapsed() >= time_limit - 0.05) return improved;
            std::vector<int> late_jobs(ranking.begin(), ranking.begin() + std::min(target_late, n));
            for (EarlyOrder mode : early_orders) {
                if (elapsed() >= time_limit - 0.05) return improved;
                if (evaluate_late_set(late_jobs, mode)) {
                    improved = true;
                    if (best_eval.late <= root_lb) return true;
                }
            }
        }
    }

    std::uint64_t seed = 7809847782465536322ULL ^ static_cast<std::uint64_t>(n * 131 + best_eval.late);
    for (const auto& job : inst.jobs) {
        seed ^= static_cast<std::uint64_t>(job.first * 911 + job.second * 3571);
        seed *= 1099511628211ULL;
    }
    std::mt19937_64 rng(seed);
    int target_late = std::max(root_lb, best_eval.late - 1);
    std::vector<int> pool = late_rankings.front();
    int pool_size = std::min(n, std::max(target_late + 20, std::min(n, 2 * target_late)));
    pool.resize(pool_size);
    std::uniform_int_distribution<int> pick_dist(0, pool_size - 1);
    for (int attempt = 0; attempt < 80; ++attempt) {
        if (elapsed() >= time_limit - 0.05) return improved;
        std::set<int> selected;
        while (static_cast<int>(selected.size()) < target_late) {
            selected.insert(pool[pick_dist(rng)]);
        }
        std::vector<int> late_jobs(selected.begin(), selected.end());
        EarlyOrder mode = early_orders[attempt % early_orders.size()];
        if (evaluate_late_set(late_jobs, mode)) {
            improved = true;
            target_late = std::max(root_lb, best_eval.late - 1);
            if (best_eval.late <= root_lb) return true;
        }
    }
    return improved;
}

struct BGMHOptions {
    bool use_global_rebuild = true;
    bool use_windows = true;
    bool strengthen_windows = true;
    bool interleave_global_rebuild = false;
    std::string algorithm_name = "Bound-guided matheuristic";
};

SolveResult solve_bound_guided_matheuristic(
    const Instance& inst,
    double time_limit,
    const BGMHOptions& options = BGMHOptions{}
) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };

    auto best_blocks = heuristic_blocks_with_indices(inst);
    best_blocks = best_blocks_for_fixed_sequence(inst, flatten_blocks(best_blocks));
    ScheduleEval best_eval = evaluate_blocks(inst, best_blocks);
    const int initial = best_eval.late;
    trace_bgmh_event(inst, "initial", 0, elapsed(), best_eval.late);
    TypedInstance typed = to_typed_instance(inst);
    int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);
    int recharge_lb = bin_packing_lower_bound(typed, typed.count);
    long long neighborhoods = 0;
    long long nodes = 0;
    double best_bound_report = static_cast<double>(root_lb);
    const double local_time_limit = time_limit;

    for (int round = 0; elapsed() < local_time_limit - EPS; ++round) {
        if (best_eval.late <= root_lb) break;
        bool improved = false;
        if (round == 0 || options.interleave_global_rebuild) {
            double rebuild_budget = 0.0;
            if (inst.jobs.size() >= 100) rebuild_budget = (round == 0 ? 25.0 : 1.5);
            double rebuild_deadline = std::min(local_time_limit, elapsed() + rebuild_budget);
            if (options.use_global_rebuild && rebuild_deadline > elapsed() + EPS) {
                while (elapsed() < rebuild_deadline - EPS
                       && best_eval.late > root_lb
                       && global_late_set_rebuild_search(inst, best_blocks, best_eval, root_lb, start, rebuild_deadline, round)) {
                    improved = true;
                }
            }
        }
        if (improved) continue;
        if (options.use_windows) {
            auto windows = critical_windows(inst, best_blocks, round);
            for (auto [left, right] : windows) {
                double remaining = local_time_limit - elapsed();
                if (remaining <= 0.05) break;
                double per_window = inst.jobs.size() <= 70 ? 1.0 : (inst.jobs.size() <= 140 ? 1.8 : 0.8);
                if (round >= 2) per_window *= 1.5;
                if (round >= 4) per_window *= 1.5;
                per_window = std::min(per_window, remaining);
                auto cand = polish_fixed_window_cplex(
                    inst, best_blocks, left, right, per_window,
                    options.strengthen_windows, root_lb, recharge_lb
                );
                ++neighborhoods;
                nodes += cand.nodes;
                if (!cand.feasible) continue;
                ScheduleEval cand_eval = evaluate_blocks(inst, cand.blocks);
                if (better_schedule(cand_eval, best_eval)) {
                    best_blocks = std::move(cand.blocks);
                    best_eval = cand_eval;
                    trace_bgmh_event(inst, "window", round, elapsed(), best_eval.late);
                    improved = true;
                    break;
                }
            }
        }
        if (!improved && round >= 5) break;
    }

    trace_bgmh_event(inst, "final", -1, elapsed(), best_eval.late);
    return {
        options.algorithm_name, best_eval.late, best_bound_report,
        best_eval.late <= best_bound_report + EPS, std::max(neighborhoods, nodes), elapsed(),
        initial, time_limit, true
    };
}

SolveResult solve_cplex_full_recharge_position(
    const Instance& inst,
    double time_limit,
    bool use_heuristic_start,
    bool use_strengthening,
    const std::string& algorithm_name
) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    const int n = static_cast<int>(inst.jobs.size());
    const int E = inst.capacity;
    const double full_recharge_time = inst.alpha + inst.beta * E;
    int sum_p = 0;
    for (const auto& job : inst.jobs) sum_p += job.first;
    const double horizon_ub = sum_p + n * full_recharge_time;

    IloEnv env;
    try {
        IloModel model(env);
        std::vector<IloBoolVarArray> x;
        x.reserve(n);
        for (int j = 0; j < n; ++j) {
            x.emplace_back(env, n);
            for (int i = 0; i < n; ++i) {
                std::ostringstream name;
                name << "xf_" << j << "_" << i;
                x[j][i].setName(name.str().c_str());
            }
        }
        IloBoolVarArray z(env, n);
        IloBoolVarArray u(env, n);
        IloNumVarArray g(env, n, 0.0, E);
        IloNumVarArray ctime(env, n, 0.0, horizon_ub);

        for (int j = 0; j < n; ++j) {
            IloExpr expr(env);
            for (int i = 0; i < n; ++i) expr += x[j][i];
            model.add(expr == 1);
            expr.end();
        }
        for (int i = 0; i < n; ++i) {
            IloExpr expr(env);
            for (int j = 0; j < n; ++j) expr += x[j][i];
            model.add(expr == 1);
            expr.end();
        }

        model.add(z[0] == 1);
        for (int i = 0; i < n; ++i) {
            IloExpr pexpr(env), eexpr(env);
            for (int j = 0; j < n; ++j) {
                pexpr += inst.jobs[j].first * x[j][i];
                eexpr += inst.jobs[j].second * x[j][i];
            }

            if (i == 0) {
                model.add(ctime[i] == full_recharge_time + pexpr);
                model.add(g[i] == E - eexpr);
            } else {
                model.add(ctime[i] == ctime[i - 1] + full_recharge_time * z[i] + pexpr);
                model.add(g[i] <= g[i - 1] - eexpr + 2 * E * z[i]);
                model.add(g[i] >= g[i - 1] - eexpr - 2 * E * z[i]);
                model.add(g[i] <= E - eexpr + 2 * E * (1 - z[i]));
                model.add(g[i] >= E - eexpr - 2 * E * (1 - z[i]));
            }
            double m_i = std::max(0.0, horizon_ub - inst.due_dates[i]);
            model.add(ctime[i] <= inst.due_dates[i] + m_i * u[i]);
            pexpr.end();
            eexpr.end();
        }

        IloExpr tardy_count(env);
        for (int i = 0; i < n; ++i) tardy_count += u[i];
        model.add(IloMinimize(env, tardy_count));

        int mip_start_obj = n;
        bool has_mip_start = false;
        IloNumVarArray start_vars(env);
        IloNumArray start_vals(env);
        if (use_heuristic_start) {
            auto blocks = best_blocks_for_fixed_sequence_full_recharge(inst, effective_load_order(inst));
            mip_start_obj = evaluate_blocks_full_recharge(inst, blocks).late;
            if (use_strengthening) model.add(tardy_count <= mip_start_obj);
            std::vector<int> sequence;
            sequence.reserve(n);
            for (const auto& block : blocks) {
                for (int j : block) sequence.push_back(j);
            }
            std::map<std::pair<int, int>, std::vector<int>> group_positions;
            for (int pos = 0; pos < n; ++pos) group_positions[inst.jobs[sequence[pos]]].push_back(pos);
            for (const auto& kv : group_positions) {
                std::vector<int> jobs;
                for (int j = 0; j < n; ++j) {
                    if (inst.jobs[j] == kv.first) jobs.push_back(j);
                }
                const auto& positions = kv.second;
                for (int idx = 0; idx < static_cast<int>(positions.size()); ++idx) {
                    sequence[positions[idx]] = jobs[idx];
                }
            }

            std::vector<double> z_val(n, 0.0), u_val(n, 0.0), g_val(n, 0.0), c_val(n, 0.0);
            std::vector<std::vector<double>> x_val(n, std::vector<double>(n, 0.0));
            double time = 0.0;
            int pos = 0;
            for (const auto& block : blocks) {
                int residual = E;
                time += full_recharge_time;
                for (int local = 0; local < static_cast<int>(block.size()); ++local) {
                    int j = sequence[pos];
                    x_val[j][pos] = 1.0;
                    if (local == 0) z_val[pos] = 1.0;
                    residual -= inst.jobs[j].second;
                    g_val[pos] = residual;
                    time += inst.jobs[j].first;
                    c_val[pos] = time;
                    u_val[pos] = time > inst.due_dates[pos] + EPS ? 1.0 : 0.0;
                    ++pos;
                }
            }

            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    start_vars.add(x[j][i]);
                    start_vals.add(x_val[j][i]);
                }
            }
            for (int i = 0; i < n; ++i) {
                start_vars.add(z[i]); start_vals.add(z_val[i]);
                start_vars.add(u[i]); start_vals.add(u_val[i]);
                start_vars.add(g[i]); start_vals.add(g_val[i]);
                start_vars.add(ctime[i]); start_vals.add(c_val[i]);
            }
            has_mip_start = true;
        }
        if (use_strengthening) {
            TypedInstance typed = to_typed_instance(inst);
            int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);
            model.add(tardy_count >= root_lb);

            IloExpr recharge_count(env);
            for (int i = 0; i < n; ++i) recharge_count += z[i];
            model.add(recharge_count >= bin_packing_lower_bound(typed, typed.count));
            recharge_count.end();
        }
        tardy_count.end();

        IloCplex cplex(model);
        if (has_mip_start) {
            cplex.addMIPStart(start_vars, start_vals, IloCplex::MIPStartAuto, "full_recharge_start");
        }
        start_vars.end();
        start_vals.end();
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::TiLim, time_limit);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.0);
        cplex.setParam(IloCplex::EpAGap, 0.0);
        bool ok = cplex.solve();
        double runtime = std::chrono::duration<double>(Clock::now() - start).count();
        SolveResult res;
        res.algorithm = algorithm_name;
        res.available = true;
        res.time_limit = time_limit;
        res.runtime = runtime;
        if (ok) {
            res.objective = static_cast<int>(std::llround(cplex.getObjValue()));
            res.best_bound = cplex.getBestObjValue();
            auto status = cplex.getStatus();
            res.certified = (status == IloAlgorithm::Optimal);
            if (res.certified) res.best_bound = static_cast<double>(res.objective);
            res.nodes = static_cast<long long>(cplex.getNnodes());
            res.incumbent_initial = use_heuristic_start ? mip_start_obj : n;
            std::vector<int> sequence(n, -1);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (cplex.getValue(x[j][i]) > 0.5) {
                        sequence[i] = j;
                        break;
                    }
                }
            }
            if (std::find(sequence.begin(), sequence.end(), -1) == sequence.end()) {
                std::vector<std::vector<int>> blocks;
                for (int i = 0; i < n; ++i) {
                    if (i == 0 || cplex.getValue(z[i]) > 0.5 || blocks.empty()) {
                        blocks.push_back({});
                    }
                    blocks.back().push_back(sequence[i]);
                }
                attach_schedule_metrics(res, inst, blocks, true);
            }
        } else {
            res.objective = n;
            res.best_bound = 0.0;
            res.certified = false;
            res.nodes = static_cast<long long>(cplex.getNnodes());
            res.incumbent_initial = use_heuristic_start ? mip_start_obj : n;
        }
        env.end();
        return res;
    } catch (const IloException& ex) {
        std::string msg = ex.getMessage();
        env.end();
        throw std::runtime_error("CPLEX full-recharge exception: " + msg);
    }
}

SolveResult solve_cplex_position(
    const Instance& inst,
    double time_limit,
    bool use_heuristic_start,
    bool use_strengthening,
    const std::string& algorithm_name
) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    const int n = static_cast<int>(inst.jobs.size());
    const int E = inst.capacity;
    int sum_p = 0;
    for (const auto& job : inst.jobs) sum_p += job.first;
    const double horizon_ub = sum_p + n * (inst.alpha + inst.beta * E);

    IloEnv env;
    try {
        IloModel model(env);
        std::vector<IloBoolVarArray> x;
        x.reserve(n);
        for (int j = 0; j < n; ++j) {
            x.emplace_back(env, n);
            for (int i = 0; i < n; ++i) {
                std::ostringstream name;
                name << "x_" << j << "_" << i;
                x[j][i].setName(name.str().c_str());
            }
        }
        IloBoolVarArray z(env, n);
        IloBoolVarArray u(env, n);
        IloNumVarArray r(env, n, 0.0, E);
        IloNumVarArray g(env, n, 0.0, E);
        IloNumVarArray ctime(env, n, 0.0, horizon_ub);

        for (int j = 0; j < n; ++j) {
            IloExpr expr(env);
            for (int i = 0; i < n; ++i) expr += x[j][i];
            model.add(expr == 1);
            expr.end();
        }
        for (int i = 0; i < n; ++i) {
            IloExpr expr(env);
            for (int j = 0; j < n; ++j) expr += x[j][i];
            model.add(expr == 1);
            expr.end();
        }

        for (int i = 0; i < n; ++i) {
            IloExpr pexpr(env), eexpr(env);
            for (int j = 0; j < n; ++j) {
                pexpr += inst.jobs[j].first * x[j][i];
                eexpr += inst.jobs[j].second * x[j][i];
            }

            if (i == 0) {
                model.add(ctime[i] == inst.alpha * z[i] + inst.beta * r[i] + pexpr);
                model.add(g[i] == r[i] - eexpr);
                model.add(r[i] <= E * z[i]);
            } else {
                model.add(ctime[i] == ctime[i - 1] + inst.alpha * z[i] + inst.beta * r[i] + pexpr);
                model.add(g[i] == g[i - 1] + r[i] - eexpr);
                model.add(g[i - 1] + r[i] <= E);
                model.add(r[i] <= E * z[i]);
            }
            double m_i = std::max(0.0, horizon_ub - inst.due_dates[i]);
            model.add(ctime[i] <= inst.due_dates[i] + m_i * u[i]);
            pexpr.end();
            eexpr.end();
        }

        IloExpr tardy_count(env);
        for (int i = 0; i < n; ++i) tardy_count += u[i];
        model.add(IloMinimize(env, tardy_count));

        int mip_start_obj = n;
        bool has_mip_start = false;
        IloNumVarArray start_vars(env);
        IloNumArray start_vals(env);
        if (use_heuristic_start) {
            auto blocks = rigorous_start_blocks_with_indices(inst);
            mip_start_obj = evaluate_blocks(inst, blocks).late;
            if (use_strengthening) model.add(tardy_count <= mip_start_obj);
            std::vector<int> sequence;
            sequence.reserve(n);
            for (const auto& block : blocks) {
                for (int j : block) sequence.push_back(j);
            }
            std::map<std::pair<int, int>, std::vector<int>> group_positions;
            for (int pos = 0; pos < n; ++pos) group_positions[inst.jobs[sequence[pos]]].push_back(pos);
            for (const auto& kv : group_positions) {
                std::vector<int> jobs;
                for (int j = 0; j < n; ++j) {
                    if (inst.jobs[j] == kv.first) jobs.push_back(j);
                }
                const auto& positions = kv.second;
                for (int idx = 0; idx < static_cast<int>(positions.size()); ++idx) {
                    sequence[positions[idx]] = jobs[idx];
                }
            }

            std::vector<double> z_val(n, 0.0), u_val(n, 0.0), r_val(n, 0.0), g_val(n, 0.0), c_val(n, 0.0);
            std::vector<std::vector<double>> x_val(n, std::vector<double>(n, 0.0));

            double time = 0.0;
            int pos = 0;
            for (const auto& block : blocks) {
                int block_e = 0;
                for (int j : block) block_e += inst.jobs[j].second;
                int residual = block_e;
                time += inst.alpha + inst.beta * block_e;
                for (int local = 0; local < static_cast<int>(block.size()); ++local) {
                    int j = sequence[pos];
                    x_val[j][pos] = 1.0;
                    if (local == 0) {
                        z_val[pos] = 1.0;
                        r_val[pos] = block_e;
                    }
                    residual -= inst.jobs[j].second;
                    g_val[pos] = residual;
                    time += inst.jobs[j].first;
                    c_val[pos] = time;
                    u_val[pos] = time > inst.due_dates[pos] + EPS ? 1.0 : 0.0;
                    ++pos;
                }
            }

            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    start_vars.add(x[j][i]);
                    start_vals.add(x_val[j][i]);
                }
            }
            for (int i = 0; i < n; ++i) {
                start_vars.add(z[i]); start_vals.add(z_val[i]);
                start_vars.add(u[i]); start_vals.add(u_val[i]);
                start_vars.add(r[i]); start_vals.add(r_val[i]);
                start_vars.add(g[i]); start_vals.add(g_val[i]);
                start_vars.add(ctime[i]); start_vals.add(c_val[i]);
            }
            has_mip_start = true;
        }
        if (use_strengthening) {
            TypedInstance typed = to_typed_instance(inst);
            int root_lb = node_lower_bound(typed, typed.count, 0, 0.0, true);
            model.add(tardy_count >= root_lb);

            IloExpr recharge_count(env);
            for (int i = 0; i < n; ++i) recharge_count += z[i];
            model.add(recharge_count >= bin_packing_lower_bound(typed, typed.count));
            recharge_count.end();

            model.add(g[n - 1] == 0);
        }
        tardy_count.end();

        IloCplex cplex(model);
        if (has_mip_start) {
            cplex.addMIPStart(start_vars, start_vals, IloCplex::MIPStartAuto, "relaxation_dp_start");
        }
        start_vars.end();
        start_vals.end();
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::TiLim, time_limit);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.0);
        cplex.setParam(IloCplex::EpAGap, 0.0);
        bool ok = cplex.solve();
        double runtime = std::chrono::duration<double>(Clock::now() - start).count();
        SolveResult res;
        res.algorithm = algorithm_name;
        res.available = true;
        res.time_limit = time_limit;
        res.runtime = runtime;
        if (ok) {
            res.objective = static_cast<int>(std::llround(cplex.getObjValue()));
            res.best_bound = cplex.getBestObjValue();
            auto status = cplex.getStatus();
            res.certified = (status == IloAlgorithm::Optimal);
            if (res.certified) res.best_bound = static_cast<double>(res.objective);
            res.nodes = static_cast<long long>(cplex.getNnodes());
            res.incumbent_initial = use_heuristic_start ? mip_start_obj : n;
            std::vector<int> sequence(n, -1);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (cplex.getValue(x[j][i]) > 0.5) {
                        sequence[i] = j;
                        break;
                    }
                }
            }
            if (std::find(sequence.begin(), sequence.end(), -1) == sequence.end()) {
                std::vector<std::vector<int>> blocks;
                for (int i = 0; i < n; ++i) {
                    if (i == 0 || cplex.getValue(z[i]) > 0.5 || blocks.empty()) {
                        blocks.push_back({});
                    }
                    blocks.back().push_back(sequence[i]);
                }
                attach_schedule_metrics(res, inst, blocks);
            }
        } else {
            res.objective = n;
            res.best_bound = 0.0;
            res.certified = false;
            res.nodes = static_cast<long long>(cplex.getNnodes());
            res.incumbent_initial = use_heuristic_start ? mip_start_obj : n;
        }
        env.end();
        return res;
    } catch (const IloException& ex) {
        std::string msg = ex.getMessage();
        env.end();
        throw std::runtime_error("CPLEX exception: " + msg);
    }
}

SolveResult solve_milp_cplex(const Instance& inst, double time_limit) {
    return solve_cplex_position(inst, time_limit, false, false, "Plain CPLEX MILP");
}

SolveResult solve_strengthened_milp_cplex(const Instance& inst, double time_limit) {
    return solve_cplex_position(inst, time_limit, false, true, "Strengthened CPLEX MILP");
}

SolveResult solve_milp_with_heuristic_cplex(const Instance& inst, double time_limit) {
    return solve_cplex_position(
        inst, time_limit, true, true, "Strengthened CPLEX MILP heuristic-start"
    );
}

SolveResult solve_full_recharge_milp_with_heuristic_cplex(const Instance& inst, double time_limit) {
    return solve_cplex_full_recharge_position(
        inst, time_limit, true, true, "Full-recharge strengthened CPLEX MILP heuristic-start"
    );
}

#else
SolveResult solve_milp_cplex(const Instance& inst, double time_limit) {
    SolveResult res;
    res.algorithm = "CPLEX MILP";
    res.objective = static_cast<int>(inst.jobs.size());
    res.certified = false;
    res.nodes = 0;
    res.runtime = 0.0;
    res.incumbent_initial = static_cast<int>(inst.jobs.size());
    res.time_limit = time_limit;
    res.available = false;
    return res;
}

SolveResult solve_milp_with_heuristic_cplex(const Instance& inst, double time_limit) {
    SolveResult res = solve_milp_cplex(inst, time_limit);
    res.algorithm = "Strengthened CPLEX MILP heuristic-start";
    return res;
}

SolveResult solve_strengthened_milp_cplex(const Instance& inst, double time_limit) {
    SolveResult res = solve_milp_cplex(inst, time_limit);
    res.algorithm = "Strengthened CPLEX MILP";
    return res;
}

SolveResult solve_full_recharge_milp_with_heuristic_cplex(const Instance& inst, double time_limit) {
    SolveResult res = solve_milp_cplex(inst, time_limit);
    res.algorithm = "Full-recharge strengthened CPLEX MILP heuristic-start";
    return res;
}

SolveResult solve_bound_guided_matheuristic(const Instance& inst, double time_limit) {
    return solve_constructive_heuristic(inst, time_limit);
}

#endif

void print_result_csv_row(
    std::ostream& os, const ExperimentConfig& cfg, const Instance& inst,
    const TypedInstance& typed, const SolveResult& res
) {
    const auto& actual_cfg = inst.config;
    int sum_p = 0, sum_e = 0;
    for (const auto& job : inst.jobs) {
        sum_p += job.first;
        sum_e += job.second;
    }
    os << std::fixed << std::setprecision(6);
    os << cfg.n << ','
       << cfg.theta << ','
       << cfg.due_range << ','
       << cfg.sigma << ','
       << cfg.instance_id << ','
       << inst.capacity << ','
       << typed.k_types << ','
       << sum_p << ','
       << sum_e << ','
       << '"' << res.algorithm << '"' << ','
       << res.objective << ','
       << res.best_bound << ','
       << (res.certified ? 1 : 0) << ','
       << res.nodes << ','
       << res.runtime << ','
       << res.incumbent_initial << ','
       << res.time_limit << ','
       << (res.available ? 1 : 0) << ','
       << res.recharges << ','
       << res.makespan << ','
       << res.total_recharge_time << ','
       << res.battery_utilization << ','
       << actual_cfg.recharge_rho << ','
       << actual_cfg.alpha << ','
       << actual_cfg.beta << ','
       << '"' << actual_cfg.pe_correlation << '"'
       << '\n';
}

void write_instance_json(std::ostream& os, const Instance& inst) {
    const auto& cfg = inst.config;
    int sum_p = 0, sum_e = 0;
    int max_e = 0;
    for (const auto& job : inst.jobs) {
        sum_p += job.first;
        sum_e += job.second;
        max_e = std::max(max_e, job.second);
    }
    int due_capacity = inst.capacity;
    if (cfg.due_reference_sigma > 0.0) {
        due_capacity = std::max(
            max_e,
            static_cast<int>(std::ceil(cfg.due_reference_sigma * max_e))
        );
    }
    auto reference_jobs_for_export = [&]() {
        if (cfg.due_reference_correlation.empty()) return inst.jobs;
        std::vector<int> p_values;
        std::vector<int> e_values;
        p_values.reserve(inst.jobs.size());
        e_values.reserve(inst.jobs.size());
        for (const auto& job : inst.jobs) {
            p_values.push_back(job.first);
            e_values.push_back(job.second);
        }
        if (cfg.due_reference_correlation == "positive" || cfg.due_reference_correlation == "pos") {
            std::sort(p_values.begin(), p_values.end());
            std::sort(e_values.begin(), e_values.end());
        } else if (cfg.due_reference_correlation == "negative" || cfg.due_reference_correlation == "neg") {
            std::sort(p_values.begin(), p_values.end());
            std::sort(e_values.begin(), e_values.end(), std::greater<int>());
        } else if (cfg.due_reference_correlation == "independent" || cfg.due_reference_correlation == "ind") {
            std::sort(p_values.begin(), p_values.end());
            std::sort(e_values.begin(), e_values.end());
            std::mt19937_64 pair_rng(static_cast<std::uint64_t>(instance_seed(cfg) + 900001LL));
            std::shuffle(e_values.begin(), e_values.end(), pair_rng);
        } else {
            throw std::runtime_error("Unknown due-reference correlation mode: " + cfg.due_reference_correlation);
        }
        std::vector<std::pair<int, int>> out;
        out.reserve(inst.jobs.size());
        for (size_t i = 0; i < inst.jobs.size(); ++i) out.emplace_back(p_values[i], e_values[i]);
        return out;
    };
    double due_horizon = heuristic_makespan(reference_jobs_for_export(), due_capacity, inst.alpha, inst.beta);
    os << std::fixed << std::setprecision(6);
    os << "{\n";
    os << "  \"problem\": \"1|ren-e,T=alpha+beta r,GDD|sum U_j\",\n";
    os << "  \"n\": " << cfg.n << ",\n";
    os << "  \"tau\": " << cfg.theta << ",\n";
    os << "  \"R\": " << cfg.due_range << ",\n";
    os << "  \"gamma\": " << cfg.sigma << ",\n";
    os << "  \"rep\": " << cfg.instance_id << ",\n";
    os << "  \"seed\": " << instance_seed(cfg) << ",\n";
    os << "  \"p_range\": [" << cfg.p_min << ", " << cfg.p_max << "],\n";
    os << "  \"e_range\": [" << cfg.e_min << ", " << cfg.e_max << "],\n";
    os << "  \"processing_energy_correlation\": \"" << cfg.pe_correlation << "\",\n";
    os << "  \"due_reference_correlation\": \"" << cfg.due_reference_correlation << "\",\n";
    os << "  \"alpha\": " << inst.alpha << ",\n";
    os << "  \"beta\": " << inst.beta << ",\n";
    os << "  \"rho\": " << cfg.recharge_rho << ",\n";
    os << "  \"capacity\": " << inst.capacity << ",\n";
    if (cfg.due_reference_sigma > 0.0) {
        os << "  \"due_reference_gamma\": " << cfg.due_reference_sigma << ",\n";
        os << "  \"due_reference_capacity\": " << due_capacity << ",\n";
    }
    os << "  \"sum_p\": " << sum_p << ",\n";
    os << "  \"sum_e\": " << sum_e << ",\n";
    os << "  \"due_date_generation\": {\n";
    os << "    \"mode\": \"" << cfg.due_mode << "\",\n";
    os << "    \"reference_horizon\": " << due_horizon << ",\n";
    os << "    \"lower\": "
       << static_cast<int>(std::ceil(due_horizon * (1.0 - cfg.theta - cfg.due_range / 2.0))) << ",\n";
    os << "    \"upper\": "
       << static_cast<int>(std::floor(due_horizon * (1.0 - cfg.theta + cfg.due_range / 2.0))) << "\n";
    os << "  },\n";
    os << "  \"jobs\": [\n";
    for (size_t i = 0; i < inst.jobs.size(); ++i) {
        os << "    {\"id\": " << (i + 1)
           << ", \"p\": " << inst.jobs[i].first
           << ", \"e\": " << inst.jobs[i].second << "}";
        os << (i + 1 == inst.jobs.size() ? "\n" : ",\n");
    }
    os << "  ],\n";
    os << "  \"due_dates\": [";
    for (size_t i = 0; i < inst.due_dates.size(); ++i) {
        if (i > 0) os << ", ";
        os << inst.due_dates[i];
    }
    os << "]\n";
    os << "}\n";
}

struct Args {
    std::string mode = "part1";
    std::string solver = "matheur";
    std::string out;
    ExperimentConfig cfg;
    std::vector<int> n_values{50, 100, 150};
    std::vector<double> theta_values{0.5, 0.6};
    std::vector<double> R_values{0.8};
    std::vector<double> sigma_values{3.0, 5.0};
    int replications = 10;
    double time_limit = 30.0;
};

std::vector<std::string> split_csv(const std::string& text) {
    std::vector<std::string> out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

std::vector<int> parse_int_list(const std::string& text) {
    std::vector<int> out;
    for (const auto& item : split_csv(text)) out.push_back(std::stoi(item));
    if (out.empty()) throw std::runtime_error("Empty integer list");
    return out;
}

std::vector<double> parse_double_list(const std::string& text) {
    std::vector<double> out;
    for (const auto& item : split_csv(text)) out.push_back(std::stod(item));
    if (out.empty()) throw std::runtime_error("Empty numeric list");
    return out;
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value after " + name);
            return argv[++i];
        };
        if (key == "--mode") args.mode = need(key);
        else if (key == "--solver") args.solver = need(key);
        else if (key == "--out") args.out = need(key);
        else if (key == "--n") args.cfg.n = std::stoi(need(key));
        else if (key == "--theta" || key == "--tau") args.cfg.theta = std::stod(need(key));
        else if (key == "--R" || key == "--due-range") {
            args.cfg.due_range = std::stod(need(key));
            args.R_values = {args.cfg.due_range};
        }
        else if (key == "--due-mode") args.cfg.due_mode = need(key);
        else if (key == "--sigma" || key == "--gamma") args.cfg.sigma = std::stod(need(key));
        else if (key == "--p-min") args.cfg.p_min = std::stoi(need(key));
        else if (key == "--p-max") args.cfg.p_max = std::stoi(need(key));
        else if (key == "--e-min") args.cfg.e_min = std::stoi(need(key));
        else if (key == "--e-max") args.cfg.e_max = std::stoi(need(key));
        else if (key == "--alpha") args.cfg.alpha = std::stod(need(key));
        else if (key == "--beta") args.cfg.beta = std::stod(need(key));
        else if (key == "--rho") args.cfg.recharge_rho = std::stod(need(key));
        else if (key == "--pe-correlation" || key == "--correlation") args.cfg.pe_correlation = need(key);
        else if (key == "--due-reference-correlation") args.cfg.due_reference_correlation = need(key);
        else if (key == "--common-base-seeds") args.cfg.common_base_seed = true;
        else if (key == "--due-reference-gamma") args.cfg.due_reference_sigma = std::stod(need(key));
        else if (key == "--n-values") args.n_values = parse_int_list(need(key));
        else if (key == "--theta-values" || key == "--tau-values") args.theta_values = parse_double_list(need(key));
        else if (key == "--R-values" || key == "--due-range-values") args.R_values = parse_double_list(need(key));
        else if (key == "--sigma-values" || key == "--gamma-values") args.sigma_values = parse_double_list(need(key));
        else if (key == "--rep") args.cfg.instance_id = std::stoi(need(key));
        else if (key == "--reps") args.replications = std::stoi(need(key));
        else if (key == "--time") args.time_limit = std::stod(need(key));
        else if (key == "--help") {
            cout << "Usage: caie_part1_final --mode one|part1|export-instance --solver heuristic|bgmh|grasp-dp|grasp-pr|grasp-pr-full|grasp-lns|grasp-lb|milp|milp-s|milp-heur|milp-full-heur [options]\n"
                 << "Diagnostic solvers: h0-raw, bgmh-no-cplex, bgmh-no-rebuild, bgmh-plain-window, bgmh-interleave-rebuild, grasp-dp, grasp-pr, grasp-pr-long, grasp-lns, grasp-lb\n"
                 << "Generation options: --due-mode theta|taur --p-min --p-max --e-min --e-max --alpha --beta --rho --pe-correlation independent|positive|negative --theta/--tau --sigma/--gamma\n"
                 << "Grid options: --n-values --theta-values/--tau-values --sigma-values/--gamma-values --reps\n"
                 << "Export option: --mode export-instance --out path/to/instance.json\n"
                 << "Seed option: --common-base-seeds; sensitivity options: --due-reference-gamma --due-reference-correlation\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + key);
        }
    }
    return args;
}

std::vector<SolveResult> solve_requested(const Instance& inst, const Args& args) {
    std::vector<SolveResult> results;
    if (args.solver == "heuristic" || args.solver == "cbh") {
        results.push_back(solve_constructive_heuristic(inst, args.time_limit));
    }
    if (args.solver == "heuristic-raw" || args.solver == "h0-raw") {
        results.push_back(solve_constructive_heuristic_no_resegmentation(inst, args.time_limit));
    }
    if (args.solver == "matheur" || args.solver == "bgmh") {
        results.push_back(solve_bound_guided_matheuristic(inst, args.time_limit));
    }
    if (args.solver == "grasp-lns" || args.solver == "alns") {
#ifdef WITH_CPLEX
        results.push_back(solve_grasp_lns_matheuristic(inst, args.time_limit));
#else
        throw std::runtime_error("grasp-lns requires a CPLEX build");
#endif
    }
    if (args.solver == "grasp-dp") {
        results.push_back(solve_grasp_dp_matheuristic(inst, args.time_limit));
    }
    if (args.solver == "grasp-pr") {
        results.push_back(solve_grasp_pr_matheuristic(inst, args.time_limit));
    }
    if (args.solver == "grasp-pr-full" || args.solver == "full-grasp-pr") {
        results.push_back(solve_grasp_pr_matheuristic(inst, args.time_limit, false, true));
    }
    if (args.solver == "grasp-pr-long") {
        results.push_back(solve_grasp_pr_matheuristic(inst, args.time_limit, true));
    }
    if (args.solver == "grasp-lb" || args.solver == "local-branching") {
#ifdef WITH_CPLEX
        results.push_back(solve_grasp_local_branching_matheuristic(inst, args.time_limit));
#else
        throw std::runtime_error("grasp-lb requires a CPLEX build");
#endif
    }
    if (args.solver == "bgmh-no-rebuild") {
        BGMHOptions options;
        options.use_global_rebuild = false;
        options.algorithm_name = "BGMH without global rebuild";
        results.push_back(solve_bound_guided_matheuristic(inst, args.time_limit, options));
    }
    if (args.solver == "bgmh-no-cplex") {
        BGMHOptions options;
        options.use_windows = false;
        options.algorithm_name = "BGMH without CPLEX neighborhoods";
        results.push_back(solve_bound_guided_matheuristic(inst, args.time_limit, options));
    }
    if (args.solver == "bgmh-plain-window") {
        BGMHOptions options;
        options.strengthen_windows = false;
        options.algorithm_name = "BGMH with plain windows";
        results.push_back(solve_bound_guided_matheuristic(inst, args.time_limit, options));
    }
    if (args.solver == "bgmh-interleave-rebuild") {
        BGMHOptions options;
        options.interleave_global_rebuild = true;
        options.algorithm_name = "BGMH interleaved rebuild";
        results.push_back(solve_bound_guided_matheuristic(inst, args.time_limit, options));
    }
    if (args.solver == "milp") {
        results.push_back(solve_milp_cplex(inst, args.time_limit));
    }
    if (args.solver == "milp-s" || args.solver == "milp-strengthened") {
        results.push_back(solve_strengthened_milp_cplex(inst, args.time_limit));
    }
    if (args.solver == "milp-heur") {
        results.push_back(solve_milp_with_heuristic_cplex(inst, args.time_limit));
    }
    if (args.solver == "milp-full" || args.solver == "full-recharge" ||
        args.solver == "milp-full-heur" || args.solver == "full-milp-heur") {
        results.push_back(solve_full_recharge_milp_with_heuristic_cplex(inst, args.time_limit));
    }
    if (results.empty()) {
        throw std::runtime_error(
            "unknown solver"
        );
    }
    return results;
}

int run_one(const Args& args) {
    Instance inst = generate_instance(args.cfg);
    TypedInstance typed = to_typed_instance(inst);
    auto results = solve_requested(inst, args);
    cout << "n,theta,R,sigma,rep,capacity,actual_types,sum_p,sum_e,algorithm,objective,best_bound,certified,nodes,runtime,incumbent_initial,time_limit,available,recharges,makespan,total_recharge_time,battery_utilization,rho,alpha,beta,pe_correlation\n";
    for (const auto& res : results) print_result_csv_row(cout, args.cfg, inst, typed, res);
    return 0;
}

int run_part1(const Args& args) {
    std::ostream* out = &cout;
    std::ofstream file;
    if (!args.out.empty()) {
        file.open(args.out);
        if (!file) throw std::runtime_error("Cannot open output file: " + args.out);
        out = &file;
    }
    *out << "n,theta,R,sigma,rep,capacity,actual_types,sum_p,sum_e,algorithm,objective,best_bound,certified,nodes,runtime,incumbent_initial,time_limit,available,recharges,makespan,total_recharge_time,battery_utilization,rho,alpha,beta,pe_correlation\n";
    for (int n : args.n_values) {
        for (double theta : args.theta_values) {
            for (double R : args.R_values) {
                for (double sigma : args.sigma_values) {
                    for (int rep = 0; rep < args.replications; ++rep) {
                        ExperimentConfig cfg = args.cfg;
                        cfg.n = n;
                        cfg.theta = theta;
                        cfg.due_range = R;
                        cfg.sigma = sigma;
                        cfg.instance_id = rep;
                        Instance inst = generate_instance(cfg);
                        TypedInstance typed = to_typed_instance(inst);
                        for (const auto& res : solve_requested(inst, args)) {
                            print_result_csv_row(*out, cfg, inst, typed, res);
                        }
                        out->flush();
                    }
                }
            }
        }
    }
    return 0;
}

int export_instance(const Args& args) {
    Instance inst = generate_instance(args.cfg);
    if (args.out.empty()) {
        write_instance_json(cout, inst);
        return 0;
    }
    std::ofstream file(args.out);
    if (!file) throw std::runtime_error("Cannot open output file: " + args.out);
    write_instance_json(file, inst);
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        Args args = parse_args(argc, argv);
        if (args.mode == "one") return run_one(args);
        if (args.mode == "part1") return run_part1(args);
        if (args.mode == "export-instance") return export_instance(args);
        throw std::runtime_error("mode must be one, part1, or export-instance");
    } catch (const std::exception& ex) {
        cerr << "Error: " << ex.what() << endl;
        return 1;
    }
}
