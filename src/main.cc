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
    std::vector<std::vector<double>> prof; /* profile matrix */
    std::vector<double> null0; /* total character frequency */
    std::string const alphs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"; /* allowed alphabets */
    int mlen; /* motif length */
    int sampler(std::vector<double> const &d) /* sample from given histogram */
    {
        std::random_device rd;
        std::mt19937 generator(rd());
        std::discrete_distribution<int> distribution(d.begin(), d.end());
        return distribution(generator);
    }
    Gsmf(std::vector<std::string> const &s, int const m): mlen(m)
    {
        nseqs.resize(s.size());
        startPos.resize(s.size(), 0.);
        null0.resize(alphs.size(), 1. /* pseudocount */);
        for (size_t i = 0; i < s.size(); ++ i) {
            nseqs[i].resize(s[i].size());
            std::transform(s[i].begin(), s[i].end(), nseqs[i].begin(), [&](char c){return alphs.find(c);});
            std::for_each(nseqs[i].begin(), nseqs[i].end(), [&](int c){++ null0[c];});
        }
        prof = std::vector<std::vector<double>>(m + 1, std::vector<double>(alphs.size(), 1. / alphs.size()));
    }
    void normalize(std::vector<double> v)
    {
        double sum = 0.;
        for (size_t i = 0; i < v.size(); ++ i) sum += v[i];
        for (size_t i = 0; i < v.size(); ++ i) v[i] /= sum;
    }
    void train(unsigned int const maxIt)
    {
        srand(static_cast<unsigned int>(time(NULL)));
        for (size_t i = 0; i < nseqs.size(); ++ i) startPos[i] = rand() % (nseqs[i].size() - mlen + 1);
        for (size_t i = 0, count = 0; count < 10 + nseqs.size() && i < maxIt; ++ i) {
            int hold = rand() % nseqs.size(); /* choose one sequence without which profile table is constructed */
            std::fill_n(prof.begin(), mlen + 1, std::vector<double>(alphs.size(), 1. /* pseudocount */));
            prof[0] = null0;
            for (size_t j = 0; j < nseqs.size(); ++ j) /* update profile and null */
                if (hold != (int)j)
                    for (int k = 0; k < mlen; ++ k) {
                        ++ prof[k + 1][nseqs[j][startPos[j] + k]];
                        -- prof[0][nseqs[j][startPos[j] + k]];
                    }
            for (size_t j = 0; j < nseqs[hold].size(); ++ j) -- prof[0][nseqs[hold][j]]; /* update null */
            for (size_t j = 0; j < prof.size(); ++ j) normalize(prof[j]);
            std::vector<double> weightNew(nseqs[hold].size() - mlen + 1, 1.);
            for (size_t j = 0; j < nseqs[hold].size() - mlen + 1; ++ j)
                for (int k = 0; k < mlen; ++ k)
                    weightNew[j] *= prof[k + 1][nseqs[hold][j + k]] / prof[0][nseqs[hold][j + k]]; /* weight of startPos */
            double old = startPos[hold];
            count = old == (startPos[hold] = sampler(weightNew))? count + 1: 0; /* simple criteria of convergence */
        }
    }
};

int main(int const, char const * const argv[])
{
    std::ifstream fasta(argv[1]);
    std::vector<std::string> seqs;
    for (std::string s; getline(fasta, s);)
        if ('>' == s.at(0))
            seqs.push_back("");
        else
            seqs.back() += s;
    int const mlen = atoi(argv[2]);
    Gsmf g(seqs, mlen);
    g.train(1000000);
    for (size_t i = 0; i < seqs.size(); ++ i)
        std::cerr << g.startPos[i] << "\t"<< seqs[i].substr(g.startPos[i], mlen) << std::endl;
    return 0;
}
