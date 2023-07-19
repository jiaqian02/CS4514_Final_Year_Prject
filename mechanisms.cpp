#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#include <vector>
#include <cstdio>
#include <iostream>

struct frac_t {
	int a, b;
	frac_t(int a = 0, int b = 1) : a(a), b(b) {
		assert(b != 0);
		int g = std::__gcd(a, b);
		a /= g;
		b /= g;
		if (b < 0) {
			a = -a;
			b = -b;
		}
	}
	frac_t operator +(const frac_t &f) const {
		return frac_t(a * f.b + b * f.a, b * f.b);
	}
	frac_t operator *(const frac_t &f) const {
		return frac_t(a * f.a, b * f.b);
	}
	bool operator <(const frac_t &f) const {
		return a * f.b < b * f.a;
	}
	bool operator ==(const frac_t &f) const {
		return a * f.b == b * f.a;
	}
	double value() {
		return (double)a / b;
	}
};

int n, m;
std::vector<std::vector<frac_t> > G;

/**
 * @brief Compute the total cost of group $G_g$ at location $x$
 *
 * @param[in] g - the group index
 * @param[in] x - the location
 *
 * @return the total cost of group $G_g$ at location $x$
 */
double h(int g, double x) {
	double ans = 0;
	for (int i = 0; i < G[g].size(); i++)
		ans += 1 - fabs(G[g][i].value() - x);
	return ans;
}

/**
 * @brief Compute the maximum total group cost at location $x$
 *
 * @param[in] x - the location
 *
 * @return the maximum total group cost at location $x$
 */
double max_h(double x) {
	double ans = 0;
	for (int i = 1; i <= m; i++)
		ans = std::max(ans, h(i, x));
	return ans;
}

/**
 * @brief Compute the optimal location
 *
 * @return the optimal location
 */
double find_optimal() {
    double ans = n, ansloc = 0;
	for (double i = 0; i <= 1; i += 0.00001) {
		if (max_h(i) < ans) {
			ans = max_h(i);
			ansloc = i;
		}
	}
    return ansloc;
}

/**
 * @brief Compute the approximation ratio
 *
 * @param[in] p - the probability of putting the facility at location 0
 *
 * @return void
 */
void approx_ratio(double p) {
    double exp_cost = p * max_h(0) + (1 - p) * max_h(1);
    double opt_loc = find_optimal(), opt = max_h(opt_loc);

    if (opt <= 1e-10 && exp_cost >= 1e-10)
        std::cout << "The approximaton ratio is INF." << std::endl;
    else if (opt <= 1e-10 && exp_cost <= 1e-10) {
		std::cout << "The apprximation ratio is 1." << std::endl;
	} else {
        std::cout << "The approximaton ratio is " << exp_cost / opt << "." << std::endl;
    }
}

/**
 * @brief Run the PEPM algorithm and compute the approximation ratio
 *
 * @return void
 */
void PEPM() {
	std::cin >> n >> m;
	G.resize(m + 1);
	for (int i = 1; i <= n; i++) {
		int a, b, g;
		std::cin >> a >> b >> g;
		assert(a <= b && 0 <= g && g <= m);
		G[g].push_back(frac_t(a, b));
	}
	for (int i = 1; i <= m; i++)
		std::sort(G[i].begin(), G[i].end());

    std::vector<std::pair<int, int> > V;
    V.push_back(std::make_pair(0, 0));
    for (int i = 1; i <= m; i++) {
        int n_1 = 0, n_2 = 0;
        for (size_t j = 0; j < G[i].size(); j++) {
            if (G[i][j] < frac_t(1, 2))
                ++n_1;
            else
                ++n_2;
        }
        V.push_back(std::make_pair(n_1, n_2));
    }
    int a = 0, b = 0, c = 0, d = 0;
	double p = 0;
    for (int i = 1; i <= m; i++) {
        a = std::max(a, V[i].first + V[i].second * 2);
        b = std::max(b, V[i].first);
        c = std::max(c, V[i].first * 2 + V[i].second);
        d = std::max(d, V[i].second);
    }
    if (b == 0)
        p = 1;
    else if (d == 0)
        p = 0;
    else {
        p = ((double)a / b - 1) / ((double)a / b + (double) c / d - 2);
    }

	std::cout << "PEPM puts the facility at 0 (resp. 1) with probability " << p << " (resp. " << 1 - p << ")." << std::endl;
	approx_ratio(p);
}

/**
 * @brief Run the LGRV algorithm and compute the approximation ratio
 *
 * @return void
 */
void LGRV() {
	std::cin >> n >> m;
	G.resize(m + 1);
	for (int i = 1; i <= n; i++) {
		int a, b, g;
		std::cin >> a >> b >> g;
		assert(a <= b && 0 <= g && g <= m);
		G[g].push_back(frac_t(a, b));
	}
	for (int i = 1; i <= m; i++)
		std::sort(G[i].begin(), G[i].end());

    std::vector<std::pair<int, int> > V;
    V.push_back(std::make_pair(0, 0));
    for (int i = 1; i <= m; i++) {
        int n_1 = 0, n_2 = 0;
        for (size_t j = 0; j < G[i].size(); j++) {
            if (G[i][j] < frac_t(1, 2))
                ++n_1;
            else
                ++n_2;
        }
        V.push_back(std::make_pair(n_1, n_2));
    }
	int cur_max = 0;
	double p = 0;
    for (int i = 1; i <= m; i++) {
		if (V[i].first > cur_max) {
			cur_max = V[i].first, p = 0;
		}
		if (V[i].second > cur_max) {
			cur_max = V[i].second, p = 1;
		}
    }

	std::cout << "LGRV puts the facility at 0 (resp. 1) with probability " << p << " (resp. " << 1 - p << ")." << std::endl;
	approx_ratio(p);
}

int main() {
    PEPM();

    return 0;
}