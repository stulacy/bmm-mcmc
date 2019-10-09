#include "utils.h"

// Sample alpha from Gamma(a, b) using method described by Escobar and West
// in Section 6
// https://pdfs.semanticscholar.org/df25/adb36860c1ad9edaac04b8855a2f19e79c5b.pdf
double update_alpha(double alpha_old, double a, double b, int N, int K) {
    double b_eps, pi1, pi2, pi, alpha_new;
    b_eps = b - log(R::rbeta(alpha_old+1, N));
    pi1 = a + K - 1;
    pi2 = N * b_eps;
    pi = pi1 / (pi1 + pi2);
    alpha_new = pi * R::rgamma(a+K, 1/(b_eps)) + (1-pi) * R::rgamma(a+K-1, 1/(b_eps));
    return alpha_new;
}

void print_clusters(std::vector < std::vector < int > > clusters) {
    for (unsigned int i=0; i < clusters.size(); i++) {
        Rcout << "Cluster: " << i << "\n";
        for (auto it : clusters[i]) {
            Rcout << "Individual: " << it << "\n";
        }
    }
}

