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
 * @brief g = 0: Plot the maximum total group cost at different locations 
 *		g > 0: Plot the total cost of group $G_g$ at different locations
 *
 * @param[in] g - the group index, 0 for the maximum total group cost
 *
 * @return void
 */
void plot(int g) {
    if (g == 0) {
        std::vector<double> X, Y;
        for (double i = 0; i <= 1; i += 0.00001) {
            X.push_back(i);
            Y.push_back(max_h(i));
        }
        plt::plot(X, Y, {{"label", "Maximum total group cost"}});
    } else {
        if (G[g].size() == 0)
            return;
        std::vector<double> X_g, Y_g;
        if (G[g][0].a != 0) {
            X_g.push_back(0);
            Y_g.push_back(h(g, 0));
        }
        for (int j = 0; j < G[g].size(); j++) {
            if (j + 1 < G[g].size() && G[g][j] == G[g][j + 1])
                j++;
            X_g.push_back(G[g][j].value());
            Y_g.push_back(h(g, G[g][j].value()));
        }
        if (G[g][G[g].size() - 1].a != G[g][G[g].size() - 1].b) {
            X_g.push_back(1);
            Y_g.push_back(h(g, 1));
        }
        plt::plot(X_g, Y_g, {{"label", "Group $" + std::to_string(g) + "$'s total cost"}, {"ls", "--"}});
    }
}

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

void plot_optimal() {
	double ansloc = find_optimal(), ans = max_h(ansloc);
	plt::plot({ansloc}, {ans}, {{"marker", "*"}, {"markersize", "20"}, {"color", "gold"}, {"label", "Optimal solution location"}, {"ls", "--"}});
    plt::plot({ansloc, 0}, {ans, ans}, {{"color", "gold"}, {"ls", "--"}});
    plt::plot({ansloc, ansloc}, {ans, 0}, {{"color", "gold"}, {"ls", "--"}});
}

/**
 * @brief Read the input and plot the maximum total group cost at different locations
 *
 * @return void
 */
void mtgc_overview() {
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

    for (int i = 1; i <= m; i++)
        plot(i);
    plot(0);
    plot_optimal();

    plt::legend({{"loc", "lower right"}});
    plt::title("The maximum total group cost at different locations", {{"loc", "center"}});

    double maxh = 0;
    for (double i = 0; i <= 1; i += 0.00001)
        maxh = std::max(maxh, max_h(i));
    plt::xlim(0, 1);
    plt::ylim((double)0, maxh * 1.1);
    plt::save("fig.svg");
}

int main() {
    mtgc_overview();

    return 0;
}