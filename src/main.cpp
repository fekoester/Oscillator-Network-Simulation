// oscnet: simulate a network of SNIPER oscillators with coupling and delayed feedback.
//
//   oscnet --n 100 --b 9 --eps 0.1 --topology ring:2 --tau 2 --kappa 0.5 --t 50 --out run.csv
//   oscnet --n 50 --b 2 --b-spread 0.2 --topology all --normalize --sweep eps:0:1:21 --t 200
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "oscnet/network.hpp"

using namespace oscnet;

namespace {

struct Options {
    int n = 100;
    double b = 9.0, bSpread = 0.0, eps = 0.1, tau = 0.0, kappa = 0.0, dt = 0.001, t = 50.0;
    std::string topology = "all", out, sweep, matrix;
    bool normalize = false, writeY = false;
    long every = 10;
    uint64_t seed = 1;
    Coupling coupling = Coupling::Diffusive;
    Feedback feedback = Feedback::Pyragas;
    std::string init = "square";
};

void usage() {
    std::cout <<
        "usage: oscnet [--n 100] [--b 9] [--b-spread 0] [--eps 0.1] [--topology all|ring:K|random:P|file:PATH]\n"
        "              [--normalize] [--tau 0] [--kappa 0] [--coupling diffusive|legacy] [--feedback pyragas|legacy]\n"
        "              [--dt 0.001] [--t 50] [--seed 1] [--init square|circle]\n"
        "              [--out run.csv] [--every 10] [--write-y] [--sweep PARAM:FROM:TO:COUNT]\n"
        "Writes t, R (order parameter) and x_i (and y_i with --write-y) every --every steps.\n"
        "--sweep runs one simulation per value of eps, kappa, tau or b-spread and prints the mean R over the\n"
        "second half and the spread of the oscillators' mean frequencies.\n"
        "--init square: x, y ~ U(0, 1) as in the original; circle: random phases on the unit circle.\n";
}

Graph makeGraph(const Options &o) {
    Graph g;
    const std::string &t = o.topology;
    if (t == "all") g = Graph::allToAll(o.n);
    else if (t.rfind("ring:", 0) == 0) g = Graph::ring(o.n, std::stoi(t.substr(5)));
    else if (t.rfind("random:", 0) == 0) g = Graph::random(o.n, std::stod(t.substr(7)), o.seed + 7);
    else if (t.rfind("file:", 0) == 0) g = Graph::fromFile(t.substr(5));
    else throw std::invalid_argument("unknown topology " + t);
    if (o.normalize) g.normalise();
    return g;
}

Network build(const Options &o) {
    Graph g = makeGraph(o);
    const int n = g.n;
    std::mt19937_64 rng(o.seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Params p;
    p.eps = o.eps, p.tau = o.tau, p.kappa = o.kappa, p.dt = o.dt, p.coupling = o.coupling, p.feedback = o.feedback;
    std::vector<double> x((size_t)n), y((size_t)n);
    for (int i = 0; i < n; i++) {
        if (o.init == "circle") {
            const double phi = 2 * kPi * u(rng);
            x[(size_t)i] = std::cos(phi), y[(size_t)i] = std::sin(phi);
        } else {
            x[(size_t)i] = u(rng), y[(size_t)i] = u(rng);
        }
    }
    p.b.resize((size_t)n);
    for (int i = 0; i < n; i++) p.b[(size_t)i] = o.b + o.bSpread * (2 * u(rng) - 1);
    return Network(std::move(g), p, x, y);
}

// mean R over the second half, and the standard deviation of the mean frequencies over the second half
std::pair<double, double> summarise(const Options &o) {
    Network net = build(o);
    const long total = (long)std::llround(o.t / o.dt), half = total / 2;
    double rsum = 0;
    long rcount = 0;
    for (long s = 1; s <= total; s++) {
        net.step();
        if (s == half) net.resetFrequencies();
        if (s > half && s % 10 == 0) rsum += net.orderParameter(), rcount++;
    }
    const auto w = net.frequencies();
    double mu = 0, var = 0;
    for (double v : w) mu += v;
    mu /= (double)w.size();
    for (double v : w) var += (v - mu) * (v - mu);
    return {rcount ? rsum / rcount : net.orderParameter(), std::sqrt(var / (double)w.size())};
}

}  // namespace

int main(int argc, char **argv) {
    Options o;
    try {
        for (int i = 1; i < argc; i++) {
            std::string a = argv[i];
            auto next = [&]() -> std::string {
                if (i + 1 >= argc) throw std::invalid_argument(a + " needs a value");
                return argv[++i];
            };
            if (a == "-h" || a == "--help") return usage(), 0;
            else if (a == "--n") o.n = std::stoi(next());
            else if (a == "--b") o.b = std::stod(next());
            else if (a == "--b-spread") o.bSpread = std::stod(next());
            else if (a == "--eps") o.eps = std::stod(next());
            else if (a == "--tau") o.tau = std::stod(next());
            else if (a == "--kappa") o.kappa = std::stod(next());
            else if (a == "--dt") o.dt = std::stod(next());
            else if (a == "--t") o.t = std::stod(next());
            else if (a == "--topology") o.topology = next();
            else if (a == "--normalize") o.normalize = true;
            else if (a == "--out") o.out = next();
            else if (a == "--every") o.every = std::max(1L, std::stol(next()));
            else if (a == "--write-y") o.writeY = true;
            else if (a == "--seed") o.seed = std::stoull(next());
            else if (a == "--init") o.init = next();
            else if (a == "--sweep") o.sweep = next();
            else if (a == "--coupling") {
                const std::string v = next();
                if (v == "diffusive") o.coupling = Coupling::Diffusive;
                else if (v == "legacy") o.coupling = Coupling::Legacy;
                else throw std::invalid_argument("unknown coupling " + v);
            } else if (a == "--feedback") {
                const std::string v = next();
                if (v == "pyragas") o.feedback = Feedback::Pyragas;
                else if (v == "legacy") o.feedback = Feedback::Legacy;
                else throw std::invalid_argument("unknown feedback " + v);
            } else throw std::invalid_argument("unknown argument " + a);
        }
        if (o.init != "square" && o.init != "circle") throw std::invalid_argument("--init must be square or circle");

        if (!o.sweep.empty()) {
            std::vector<std::string> f;
            std::stringstream ss(o.sweep);
            std::string item;
            while (std::getline(ss, item, ':')) f.push_back(item);
            if (f.size() != 4) throw std::invalid_argument("--sweep expects PARAM:FROM:TO:COUNT");
            const double from = std::stod(f[1]), to = std::stod(f[2]);
            const int count = std::stoi(f[3]);
            std::printf("%s,R_mean,frequency_std\n", f[0].c_str());
            for (int k = 0; k < count; k++) {
                const double v = count == 1 ? from : from + (to - from) * k / (count - 1);
                Options q = o;
                if (f[0] == "eps") q.eps = v;
                else if (f[0] == "kappa") q.kappa = v;
                else if (f[0] == "tau") q.tau = v;
                else if (f[0] == "b-spread") q.bSpread = v;
                else throw std::invalid_argument("can only sweep eps, kappa, tau or b-spread");
                const auto [r, ws] = summarise(q);
                std::printf("%g,%.6f,%.6f\n", v, r, ws);
                std::fflush(stdout);
            }
            return 0;
        }

        Network net = build(o);
        std::ofstream out;
        if (!o.out.empty()) {
            out.open(o.out);
            if (!out) throw std::runtime_error("cannot write " + o.out);
            out << "t,R";
            for (int i = 0; i < net.n; i++) out << ",x" << i;
            if (o.writeY)
                for (int i = 0; i < net.n; i++) out << ",y" << i;
            out << '\n';
        }
        auto row = [&]() {
            out << net.time() << ',' << net.orderParameter();
            for (double v : net.x) out << ',' << v;
            if (o.writeY)
                for (double v : net.y) out << ',' << v;
            out << '\n';
        };
        if (out) row();
        const long total = (long)std::llround(o.t / o.dt);
        for (long s = 1; s <= total; s++) {
            net.step();
            if (out && s % o.every == 0) row();
        }
        const auto w = net.frequencies();
        double wmin = w[0], wmax = w[0];
        for (double v : w) wmin = std::min(wmin, v), wmax = std::max(wmax, v);
        std::printf("%d oscillators, %zu couplings, t = %g: R = %.4f, mean frequencies in [%.4f, %.4f]\n", net.n,
                    net.g.edges(), net.time(), net.orderParameter(), wmin, wmax);
    } catch (const std::exception &e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        usage();
        return 2;
    }
    return 0;
}
