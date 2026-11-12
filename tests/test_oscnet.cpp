#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cstdio>
#include <fstream>
#include <tuple>

#include "oscnet/network.hpp"

using namespace oscnet;

static Network single(double b, double dt, double x0 = 1.0, double y0 = 0.0) {
    Params p;
    p.b = {b}, p.eps = 0.0, p.dt = dt;
    return Network(Graph(1), p, {x0}, {y0});
}

// time between successive upward crossings of y = 0 at x > 0 (linear interpolation)
static double measuredPeriod(Network &net, double tmax) {
    std::vector<double> crossings;
    double prevY = net.y[0], prevT = net.time();
    while (net.time() < tmax) {
        net.step();
        const double yy = net.y[0];
        if (prevY < 0 && yy >= 0 && net.x[0] > 0)
            crossings.push_back(prevT + (net.time() - prevT) * (-prevY) / (yy - prevY));
        prevY = yy, prevT = net.time();
    }
    REQUIRE(crossings.size() >= 3);
    return (crossings.back() - crossings[1]) / double(crossings.size() - 2);
}

TEST_CASE("an uncoupled oscillator has the analytic period 2 pi / sqrt(b^2 - 1)") {
    for (double b : {1.2, 2.0, 9.0}) {
        Network net = single(b, 1e-3);
        CHECK(measuredPeriod(net, 12 * period(b)) == doctest::Approx(period(b)).epsilon(1e-6));
        CHECK(std::hypot(net.x[0], net.y[0]) == doctest::Approx(1.0).epsilon(1e-9));  // limit cycle r = 1
    }
}

TEST_CASE("below the SNIPER bifurcation the oscillator settles at the stable node") {
    // fixed points on r = 1 with cos(phi) = b; the stable one has phi = -acos(b) (phi' = b - cos(phi))
    const double b = 0.6;
    Network net = single(b, 1e-3, 0.0, 1.0);
    for (int i = 0; i < 40000; i++) net.step();
    CHECK(net.x[0] == doctest::Approx(b).epsilon(1e-6));
    CHECK(net.y[0] == doctest::Approx(-std::sqrt(1 - b * b)).epsilon(1e-6));
}

TEST_CASE("RK4 converges with order 4") {
    auto endState = [](double dt) {
        Network net = single(3.0, dt, 0.3, 0.2);
        const long steps = std::lround(2.0 / dt);
        for (long i = 0; i < steps; i++) net.step();
        return std::make_pair(net.x[0], net.y[0]);
    };
    const auto ref = endState(1e-4);
    const auto a = endState(0.02), b = endState(0.01);
    const double ea = std::hypot(a.first - ref.first, a.second - ref.second);
    const double eb = std::hypot(b.first - ref.first, b.second - ref.second);
    CHECK(std::log2(ea / eb) == doctest::Approx(4.0).epsilon(0.1));
}

static std::pair<double, double> delayedOrder(double tau, double dt1, double dt2, double dtRef) {
    auto endState = [tau](double dt) {
        Params p;
        p.b = {2.0}, p.eps = 0, p.dt = dt, p.tau = tau, p.kappa = 0.8;
        // smooth history: the exact uncoupled orbit through (1, 0) at t = 0
        Network net(Graph(1), p, {1.0}, {0.0});
        net.prefill([&](double t, int, double &x, double &y, double &fx, double &fy) {
            // on r = 1, phi' = b - cos(phi) has the closed form tan(phi/2) = sqrt((b-1)/(b+1)) tan(w t / 2)
            const double b = 2.0, w = std::sqrt(b * b - 1);
            const double phi = 2 * std::atan(std::sqrt((b - 1) / (b + 1)) * std::tan(w * t / 2));
            x = std::cos(phi), y = std::sin(phi);
            const double dphi = b - std::cos(phi);
            fx = -std::sin(phi) * dphi, fy = std::cos(phi) * dphi;
        });
        const long steps = std::lround(3.0 / dt);
        for (long i = 0; i < steps; i++) net.step();
        return std::make_pair(net.x[0], net.y[0]);
    };
    const auto ref = endState(dtRef);
    const auto a = endState(dt1), b = endState(dt2);
    return {std::hypot(a.first - ref.first, a.second - ref.second), std::hypot(b.first - ref.first, b.second - ref.second)};
}

TEST_CASE("delayed RK4 + Hermite interpolation: order 4 with tau on the grid, small errors off it") {
    // with tau a multiple of dt the derivative breakpoints t = 0, tau, 2 tau fall on grid points
    auto [ea, eb] = delayedOrder(0.75, 0.0125, 0.00625, 0.0125 / 64);
    CHECK(std::log2(ea / eb) > 3.7);
    // otherwise the kink of x(t - tau) at t = tau falls inside a step, which costs one low-order local error
    // (its size depends on where in the step the kink lies, so there is no clean order)
    std::tie(ea, eb) = delayedOrder(0.73, 0.008, 0.004, 0.008 / 64);
    CHECK(ea < 1e-6);
    CHECK(eb < 1e-6);
}

TEST_CASE("Pyragas feedback with tau = T is non-invasive on the periodic orbit") {
    const double b = 2.0, T = period(b), dt = T / 2000;  // tau exactly on the grid
    auto orbit = [&](double t, double &x, double &y, double &fx, double &fy) {
        const double w = std::sqrt(b * b - 1);
        // phi(t) with phi(0) = 0; tan has branch jumps, so unwrap via t modulo the period
        double tt = std::fmod(t, T);
        if (tt < 0) tt += T;
        const double phi = 2 * std::atan(std::sqrt((b - 1) / (b + 1)) * std::tan(w * tt / 2));
        x = std::cos(phi), y = std::sin(phi);
        const double dphi = b - std::cos(phi);
        fx = -std::sin(phi) * dphi, fy = std::cos(phi) * dphi;
    };
    Params p;
    p.b = {b}, p.eps = 0, p.dt = dt, p.tau = T, p.kappa = 2.0;
    Network controlled(Graph(1), p, {1.0}, {0.0});
    controlled.prefill([&](double t, int, double &x, double &y, double &fx, double &fy) { orbit(t, x, y, fx, fy); });
    Network free = single(b, dt);
    for (int i = 0; i < 6000; i++) controlled.step(), free.step();
    CHECK(controlled.x[0] == doctest::Approx(free.x[0]).epsilon(1e-8));
    CHECK(controlled.y[0] == doctest::Approx(free.y[0]).epsilon(1e-8));

    // with tau != T the feedback acts and the orbit changes
    p.tau = 0.7 * T;
    Network off(Graph(1), p, {1.0}, {0.0});
    off.prefill([&](double t, int, double &x, double &y, double &fx, double &fy) { orbit(t, x, y, fx, fy); });
    Network free2 = single(b, dt);
    for (int i = 0; i < 6000; i++) off.step(), free2.step();
    CHECK(std::hypot(off.x[0] - free2.x[0], off.y[0] - free2.y[0]) > 1e-3);
}

TEST_CASE("diffusive coupling synchronises identical oscillators; uncoupled ones stay apart") {
    for (double eps : {0.0, 0.5}) {
        Params p;
        p.b = {3.0}, p.eps = eps, p.dt = 1e-3;
        Graph g = Graph::allToAll(10);
        g.normalise();
        std::vector<double> x, y;
        for (int i = 0; i < 10; i++) x.push_back(std::cos(0.6 * i)), y.push_back(std::sin(0.6 * i));
        Network net(g, p, x, y);
        const double r0 = net.orderParameter();
        for (int i = 0; i < 30000; i++) net.step();
        if (eps > 0) CHECK(net.orderParameter() > 0.999);
        else CHECK(net.orderParameter() < 0.9);
        CHECK(r0 < 0.5);
    }
}

TEST_CASE("legacy coupling and feedback reproduce the 2017 right-hand side") {
    Params p;
    p.b = {9.0}, p.eps = 0.1, p.dt = 1e-3, p.coupling = Coupling::Legacy;
    Graph g(3);
    g.add(0, 1), g.add(0, 2);
    Network net(g, p, {0.2, 0.5, -0.3}, {0.1, 0.4, 0.9});
    std::vector<double> fx, fy;
    net.rhs(net.x, net.y, 0.0, fx, fy);
    const double x = 0.2, y = 0.1, r2 = x * x + y * y;
    const double c = 0.1 * ((0.5 - x) + (0.4 - y) + (-0.3 - x) + (0.9 - y));
    CHECK(fx[0] == doctest::Approx(x * (1 - r2) + y * (x - 9) + c));
    CHECK(fy[0] == doctest::Approx(y * (1 - r2) - x * (x - 9) + c));

    // legacy feedback kappa * x * x(t - tau), with a constant history equal to the initial state
    Params q;
    q.b = {9.0}, q.eps = 0, q.dt = 1e-3, q.tau = 0.5, q.kappa = 1.0, q.feedback = Feedback::Legacy;
    Network fb(Graph(1), q, {x}, {y});
    fb.rhs(fb.x, fb.y, 0.0, fx, fy);
    CHECK(fx[0] == doctest::Approx(x * (1 - r2) + y * (x - 9) + x * x));
    CHECK(fy[0] == doctest::Approx(y * (1 - r2) - x * (x - 9) + y * y));
}

TEST_CASE("graphs") {
    CHECK(Graph::allToAll(5).edges() == 20u);
    CHECK(Graph::ring(10, 2).edges() == 40u);
    CHECK(Graph::ring(3, 2).edges() == 6u);  // no duplicate neighbours on tiny rings
    const Graph r = Graph::random(200, 0.1, 3);
    CHECK(r.edges() > 3000u);
    CHECK(r.edges() < 4900u);
    for (int i = 0; i < 200; i++)
        for (const auto &[j, w] : r.neighbours[(size_t)i]) {
            bool back = false;
            for (const auto &e : r.neighbours[(size_t)j]) back |= e.first == i;
            CHECK(back);  // undirected
            CHECK(w == 1.0);
        }
    Graph n = Graph::ring(6, 1);
    n.normalise();
    CHECK(n.neighbours[0][0].second == doctest::Approx(0.5));

    {
        std::ofstream f("m_digits.txt");
        f << "011\n101\n110\n";
    }
    {
        std::ofstream f("m_numbers.txt");
        f << "0, 0.5, 0\n0 0 2\n1\t0\t0\n";
    }
    const Graph a = Graph::fromFile("m_digits.txt"), b = Graph::fromFile("m_numbers.txt");
    CHECK(a.n == 3);
    CHECK(a.edges() == 6u);
    CHECK(b.edges() == 3u);
    CHECK(b.neighbours[1][0].first == 2);
    CHECK(b.neighbours[1][0].second == 2.0);
    {
        std::ofstream f("m_bad.txt");
        f << "01\n0\n";
    }
    CHECK_THROWS(Graph::fromFile("m_bad.txt"));
    CHECK_THROWS(Graph::fromFile("does_not_exist.txt"));
    std::remove("m_digits.txt"), std::remove("m_numbers.txt"), std::remove("m_bad.txt");
}

TEST_CASE("invalid parameters are rejected") {
    Params p;
    p.b = {2.0, 2.0};
    CHECK_THROWS_AS(Network(Graph(3), p, {0, 0, 0}, {0, 0, 0}), std::invalid_argument);
    p.b = {2.0};
    p.tau = 0.0005, p.kappa = 1.0, p.dt = 0.001;
    CHECK_THROWS_AS(Network(Graph(1), p, {1}, {0}), std::invalid_argument);
}
