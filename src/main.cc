#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <ctime>
#include <random>

struct Gsmf /* Gibbs Sampling Motif Finder */
{
    std::vector<std::vector<int>> nseqs; /* sequence set */
    std::vector<int> startPos; /* start position set */
    std::vector<std::vector<double>> prof, weight; /* profile matrix and position weight matrix */
    std::string const alphs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"; /* allowed alphabets */
    int mlen; /* motif length */
    int sampler(std::vector<double> const &d)
    {
        std::random_device rd;
        std::mt19937 generator(rd());
        std::discrete_distribution<int> distribution(d.begin(), d.end());
        return distribution(generator);
    }
    void setWeight(std::vector<std::vector<double>> const &p)
    {
        for (size_t i = 0; i < nseqs.size(); ++ i)
            for (size_t j = 0; j < nseqs[i].size() - mlen + 1; ++ j)
                weight[i][j] = p[i][j]; /* avoid writing "weight = p" in case that p has larger table */
    }
    Gsmf(std::vector<std::string> const &s, int const m): mlen(m)
    {
        nseqs.resize(s.size());
        weight.resize(s.size());
        startPos.resize(s.size());
        for (size_t i = 0; i < s.size(); ++ i) {
            nseqs[i].resize(s[i].size());
            std::transform(s[i].begin(), s[i].end(), nseqs[i].begin(), [&](char c){return alphs.find(c);});
            weight[i] = std::vector<double>(s[i].size() - mlen + 1, 1.);
        }
        prof = std::vector<std::vector<double>>(m, std::vector<double>(alphs.size(), 1. / alphs.size()));
    }
    void train(unsigned int const maxIt)
    {
        srand(static_cast<unsigned int>(time(NULL)));
        std::transform(weight.begin(), weight.end(), startPos.begin(), [&](std::vector<double> &d){return sampler(d);});
        for (size_t i = 0, count = 0; count < nseqs.size() && i < maxIt; ++ i) {
            int hold = rand() % nseqs.size(); /* choose one sequence without which profile table is constructed */
            std::fill_n(prof.begin(), mlen, std::vector<double>(alphs.size(), 1. / nseqs.size() /* pseudocount */));
            for (size_t j = 0; j < nseqs.size(); ++ j)
                if (hold != (int)j)
                    for (int k = 0; k < mlen; ++ k)
                        prof[k][nseqs[j][startPos[j] + k]] += 1. / nseqs.size(); /* each row sums upto 1 (with pseudocount) */
            std::vector<double> weightNew(nseqs[hold].size() - mlen + 1, 1.);
            for (size_t j = 0; j < nseqs[hold].size() - mlen + 1; ++ j) {
                for (int k = 0; k < mlen; ++ k)
                    weightNew[j] *= prof[k][nseqs[hold][j + k]];
                weightNew[j] *= weight[hold][j]; /* profile-matrix-based posterior probability times position-weight */
            }
            double old = startPos[hold];
            count = old == (startPos[hold] = sampler(weightNew))? count + 1: 0; /* simple criteria of convergence */
        }
    }
};

int main(int const, char const * const argv[])
{
    std::ifstream fasta(argv[1]);
    std::vector<std::string> seqs;
    for (std::string s; std::getline(fasta, s);)
        if ('>' == s.at(0))
            seqs.push_back("");
        else
            seqs.back() += s;
    int const mlen = std::atoi(argv[2]);
    Gsmf g(seqs, mlen);
    g.train(1000000);
    for (size_t i = 0; i < seqs.size(); ++ i)
        std::cerr << g.startPos[i] << "\t"<< seqs[i].substr(g.startPos[i], mlen) << std::endl;
    return 0;
}
