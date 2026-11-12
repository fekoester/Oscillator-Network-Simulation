// Networks of SNIPER (active rotator) oscillators with diffusive coupling and delayed feedback.
//
//   x_i' = x_i (1 - r_i^2) + y_i (x_i - b_i) + eps * sum_j A_ij (x_j - x_i) + kappa (x_i(t - tau) - x_i(t))
//   y_i' = y_i (1 - r_i^2) - x_i (x_i - b_i) + eps * sum_j A_ij (y_j - y_i) + kappa (y_i(t - tau) - y_i(t))
//
// In polar coordinates the uncoupled oscillator is r' = r (1 - r^2), phi' = b - r cos(phi): a saddle-node on
// invariant circle (SNIPER) bifurcation at |b| = 1, excitable below it, oscillating with period
// 2 pi / sqrt(b^2 - 1) above it. Integration: classical RK4; the delayed state at off-grid times comes from
// cubic Hermite interpolation of the stored trajectory and its derivatives, so the scheme stays high-order.
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <functional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace oscnet {

constexpr double kPi = 3.14159265358979323846;

enum class Coupling { Diffusive, Legacy };   // Legacy: eps * sum_j [(x_j - x_i) + (y_j - y_i)] added to x' and y'
enum class Feedback { Pyragas, Legacy };     // Legacy: kappa * x_i * x_i(t - tau)  (and the same for y)

// Sparse coupling: neighbours[i] = list of (j, weight A_ij)
struct Graph {
    int n = 0;
    std::vector<std::vector<std::pair<int, double>>> neighbours;

    explicit Graph(int n_ = 0) : n(n_), neighbours((size_t)n_) {}
    void add(int i, int j, double w = 1.0) { neighbours.at((size_t)i).push_back({j, w}); }
    size_t edges() const {
        size_t e = 0;
        for (const auto &v : neighbours) e += v.size();
        return e;
    }
    // divide every row by its number of neighbours (mean-field coupling)
    void normalise() {
        for (auto &row : neighbours)
            for (auto &e : row) e.second /= (double)row.size();
    }

    static Graph allToAll(int n) {
        Graph g(n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j) g.add(i, j);
        return g;
    }
    // each node coupled to its k nearest neighbours on each side
    static Graph ring(int n, int k) {
        Graph g(n);
        for (int i = 0; i < n; i++) {
            std::vector<int> seen;
            auto add = [&](int j) {
                if (j == i || std::find(seen.begin(), seen.end(), j) != seen.end()) return;
                seen.push_back(j), g.add(i, j);
            };
            for (int d = 1; d <= k; d++) add((i + d) % n), add(((i - d) % n + n) % n);
        }
        return g;
    }
    // undirected Erdos-Renyi graph
    static Graph random(int n, double p, uint64_t seed) {
        Graph g(n);
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (u(rng) < p) g.add(i, j), g.add(j, i);
        return g;
    }
    // Rows of a square matrix: either digits without separators ("0110") or numbers separated by spaces,
    // tabs or commas. Nonzero entries become couplings with that weight.
    static Graph fromFile(const std::string &path) {
        std::ifstream f(path);
        if (!f) throw std::runtime_error("cannot open coupling matrix " + path);
        std::vector<std::vector<double>> rows;
        std::string line;
        while (std::getline(f, line)) {
            for (char &c : line)
                if (c == ',' || c == ';' || c == '\t' || c == '\r') c = ' ';
            if (line.find_first_not_of(' ') == std::string::npos) continue;
            std::vector<double> row;
            const auto a = line.find_first_not_of(' '), z = line.find_last_not_of(' ');
            const std::string core = line.substr(a, z - a + 1);
            if (core.find(' ') == std::string::npos && core.find_first_not_of("0123456789") == std::string::npos) {
                for (char c : core) row.push_back(c - '0');  // "0110": one digit per entry
            } else {
                std::istringstream ss(line);
                double v;
                while (ss >> v) row.push_back(v);
            }
            rows.push_back(row);
        }
        const int n = (int)rows.size();
        Graph g(n);
        for (int i = 0; i < n; i++) {
            if ((int)rows[i].size() != n)
                throw std::runtime_error(path + ": row " + std::to_string(i + 1) + " has " +
                                         std::to_string(rows[i].size()) + " entries, expected " + std::to_string(n));
            for (int j = 0; j < n; j++)
                if (rows[i][j] != 0.0) g.add(i, j, rows[i][j]);
        }
        return g;
    }
};

struct Params {
    std::vector<double> b;  // per oscillator
    double eps = 0.1;
    double tau = 0.0, kappa = 0.0;
    double dt = 0.001;
    Coupling coupling = Coupling::Diffusive;
    Feedback feedback = Feedback::Pyragas;
};

class Network {
public:
    Network(Graph graph, Params params, std::vector<double> x0, std::vector<double> y0)
        : g(std::move(graph)), p(std::move(params)), n(g.n), x(std::move(x0)), y(std::move(y0)) {
        if (n < 1) throw std::invalid_argument("empty network");
        if ((int)p.b.size() == 1) p.b.assign((size_t)n, p.b[0]);
        if ((int)p.b.size() != n || (int)x.size() != n || (int)y.size() != n)
            throw std::invalid_argument("b, x0 and y0 need one value per oscillator");
        if (!(p.dt > 0)) throw std::invalid_argument("dt must be positive");
        if (p.tau < 0) throw std::invalid_argument("tau must be >= 0");
        delayed_ = p.tau > 0 && p.kappa != 0.0;
        if (delayed_ && p.tau < p.dt) throw std::invalid_argument("tau must be 0 or at least dt");
        if (delayed_) {
            lag_ = p.tau / p.dt;
            hist_ = (size_t)std::ceil(lag_) + 3;
            hx_.assign(hist_ * n, 0.0), hy_.assign(hist_ * n, 0.0);
            hfx_.assign(hist_ * n, 0.0), hfy_.assign(hist_ * n, 0.0);
            lfx_.assign(hist_ * n, 0.0), lfy_.assign(hist_ * n, 0.0);
            // default history: constant initial state for t <= 0
            prefill([&](double, int i, double &xi, double &yi, double &fxi, double &fyi) {
                xi = x[(size_t)i], yi = y[(size_t)i], fxi = 0, fyi = 0;
            });
        }
        theta_.resize((size_t)n);
        for (int i = 0; i < n; i++) theta_[(size_t)i] = std::atan2(y[(size_t)i], x[(size_t)i]);
        theta0_ = theta_;
    }

    // Replace the history on the grid t = -k dt (k = 0 .. buffer length) by a given function that returns the
    // state and its time derivative of oscillator i at time t <= 0.
    using HistoryFn = std::function<void(double t, int i, double &x, double &y, double &fx, double &fy)>;
    void prefill(const HistoryFn &fn) {
        if (!delayed_) return;
        for (size_t k = 0; k < hist_; k++) {
            const long s = step_ - (long)k;
            const size_t slot = slotOf(s);
            for (int i = 0; i < n; i++) {
                fn(-(double)k * p.dt, i, hx_[slot * n + i], hy_[slot * n + i], hfx_[slot * n + i], hfy_[slot * n + i]);
                lfx_[slot * n + i] = hfx_[slot * n + i], lfy_[slot * n + i] = hfy_[slot * n + i];
            }
        }
    }

    // derivatives at the current state (and delayed state at stage offset c in units of dt)
    void rhs(const std::vector<double> &X, const std::vector<double> &Y, double c, std::vector<double> &FX,
             std::vector<double> &FY) const {
        FX.resize((size_t)n), FY.resize((size_t)n);
        for (int i = 0; i < n; i++) {
            const double xi = X[(size_t)i], yi = Y[(size_t)i], r2 = xi * xi + yi * yi, bi = p.b[(size_t)i];
            double fx = xi * (1 - r2) + yi * (xi - bi);
            double fy = yi * (1 - r2) - xi * (xi - bi);
            if (p.eps != 0.0) {
                double cx = 0, cy = 0;
                for (const auto &[j, w] : g.neighbours[(size_t)i]) {
                    cx += w * (X[(size_t)j] - xi);
                    cy += w * (Y[(size_t)j] - yi);
                }
                if (p.coupling == Coupling::Diffusive) fx += p.eps * cx, fy += p.eps * cy;
                else fx += p.eps * (cx + cy), fy += p.eps * (cx + cy);
            }
            if (delayed_) {
                double xd, yd;
                delayedState(i, c, xd, yd);
                if (p.feedback == Feedback::Pyragas) fx += p.kappa * (xd - xi), fy += p.kappa * (yd - yi);
                else fx += p.kappa * xi * xd, fy += p.kappa * yi * yd;
            }
            FX[(size_t)i] = fx, FY[(size_t)i] = fy;
        }
    }

    void step() {
        const double h = p.dt;
        rhs(x, y, 0.0, k1x_, k1y_);
        if (delayed_) {  // store the current point and its derivative before it is needed (c = 1 when tau == dt)
            const size_t slot = slotOf(step_);
            for (int i = 0; i < n; i++) {
                hx_[slot * n + i] = x[(size_t)i], hy_[slot * n + i] = y[(size_t)i];
                hfx_[slot * n + i] = k1x_[(size_t)i], hfy_[slot * n + i] = k1y_[(size_t)i];
                // the derivative jumps at t = 0 (history -> DDE); keep the history's one for interpolating t < 0
                if (step_ != 0) lfx_[slot * n + i] = k1x_[(size_t)i], lfy_[slot * n + i] = k1y_[(size_t)i];
            }
        }
        auto axpy = [&](const std::vector<double> &kx, const std::vector<double> &ky, double a) {
            tx_.resize((size_t)n), ty_.resize((size_t)n);
            for (int i = 0; i < n; i++) tx_[(size_t)i] = x[(size_t)i] + a * kx[(size_t)i], ty_[(size_t)i] = y[(size_t)i] + a * ky[(size_t)i];
        };
        axpy(k1x_, k1y_, h / 2), rhs(tx_, ty_, 0.5, k2x_, k2y_);
        axpy(k2x_, k2y_, h / 2), rhs(tx_, ty_, 0.5, k3x_, k3y_);
        axpy(k3x_, k3y_, h), rhs(tx_, ty_, 1.0, k4x_, k4y_);
        for (int i = 0; i < n; i++) {
            const size_t q = (size_t)i;
            x[q] += h / 6 * (k1x_[q] + 2 * k2x_[q] + 2 * k3x_[q] + k4x_[q]);
            y[q] += h / 6 * (k1y_[q] + 2 * k2y_[q] + 2 * k3y_[q] + k4y_[q]);
            // unwrapped phase for frequencies
            theta_[q] += std::remainder(std::atan2(y[q], x[q]) - theta_[q], 2 * kPi);
        }
        step_++;
    }

    double time() const { return step_ * p.dt; }
    long steps() const { return step_; }

    // Kuramoto order parameter R = |<exp(i phi_j)>|
    double orderParameter() const {
        std::complex<double> z = 0;
        for (int i = 0; i < n; i++) z += std::polar(1.0, std::atan2(y[(size_t)i], x[(size_t)i]));
        return std::abs(z) / n;
    }
    // mean angular frequency of each oscillator since the start (or since resetFrequencies())
    std::vector<double> frequencies() const {
        std::vector<double> w((size_t)n);
        const double T = (step_ - freqStart_) * p.dt;
        for (int i = 0; i < n; i++) w[(size_t)i] = T > 0 ? (theta_[(size_t)i] - theta0_[(size_t)i]) / T : 0.0;
        return w;
    }
    void resetFrequencies() { theta0_ = theta_, freqStart_ = step_; }

    Graph g;
    Params p;
    int n;
    std::vector<double> x, y;

private:
    size_t slotOf(long s) const { return (size_t)(((s % (long)hist_) + (long)hist_) % (long)hist_); }

    // state at time t_now + c dt - tau by cubic Hermite interpolation on the stored grid
    void delayedState(int i, double c, double &xd, double &yd) const {
        const double s = (double)step_ + c - lag_;
        const double fl = std::floor(s);
        long m = (long)fl;
        double th = s - fl;
        if (th < 1e-12) th = 0.0;
        const size_t a = slotOf(m) * n + i, b2 = slotOf(m + 1) * n + i;
        const double h = p.dt;
        const double t2 = th * th, t3 = t2 * th;
        const double h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + th, h01 = -2 * t3 + 3 * t2, h11 = t3 - t2;
        if (th == 0.0) {
            xd = hx_[a], yd = hy_[a];
            return;
        }
        // right-sided derivative at the left end, left-sided at the right end of the interval
        xd = h00 * hx_[a] + h10 * h * hfx_[a] + h01 * hx_[b2] + h11 * h * lfx_[b2];
        yd = h00 * hy_[a] + h10 * h * hfy_[a] + h01 * hy_[b2] + h11 * h * lfy_[b2];
    }

    bool delayed_ = false;
    double lag_ = 0;
    size_t hist_ = 0;
    long step_ = 0, freqStart_ = 0;
    std::vector<double> hx_, hy_, hfx_, hfy_, lfx_, lfy_;
    std::vector<double> k1x_, k1y_, k2x_, k2y_, k3x_, k3y_, k4x_, k4y_, tx_, ty_;
    std::vector<double> theta_, theta0_;
};

// period of the uncoupled oscillator on the unit circle (b > 1)
inline double period(double b) { return 2 * kPi / std::sqrt(b * b - 1); }

}  // namespace oscnet
