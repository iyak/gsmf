#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <ctime>
#include <random>
using namespace std;

struct Gsmf /* Gibbs Sampling Motif Finder */
{
    vector<vector<int>> nseqs; /* sequence set */
    vector<int> startPos; /* start position set */
    vector<vector<double>> prof; /* profile matrix */
    string const alphs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"; /* allowed alphabets */
    int mlen; /* motif length */
    int sampler(vector<double> const &d)
    {
        random_device rd;
        mt19937 generator(rd());
        discrete_distribution<int> distribution(d.begin(), d.end());
        return distribution(generator);
    }
    Gsmf(vector<string> const &s, int const m): mlen(m)
    {
        nseqs.resize(s.size());
        startPos.resize(s.size());
        for (size_t i = 0; i < s.size(); ++ i) {
            nseqs[i].resize(s[i].size());
            transform(s[i].begin(), s[i].end(), nseqs[i].begin(), [&](char c){return alphs.find(c);});
        }
        prof = vector<vector<double>>(m, vector<double>(alphs.size(), 1. / alphs.size()));
    }
    void train(unsigned int const maxIt)
    {
        srand(static_cast<unsigned int>(time(NULL)));
        transform(nseqs.begin(), nseqs.end(), startPos.begin(), [&](vector<int> &s){return rand() % (s.size() - mlen + 1);});
        for (size_t i = 0, count = 0; count < nseqs.size() && i < maxIt; ++ i) {
            int hold = rand() % nseqs.size(); /* choose one sequence without which profile table is constructed */
            fill_n(prof.begin(), mlen, vector<double>(alphs.size(), 1. / nseqs.size() /* pseudocount */));
            for (size_t j = 0; j < nseqs.size(); ++ j)
                if (hold != (int)j)
                    for (int k = 0; k < mlen; ++ k)
                        prof[k][nseqs[j][startPos[j] + k]] += 1. / nseqs.size(); /* each row sums upto 1 (with pseudocount) */
            vector<double> weightNew(nseqs[hold].size() - mlen + 1, 1.);
            for (size_t j = 0; j < nseqs[hold].size() - mlen + 1; ++ j)
                for (int k = 0; k < mlen; ++ k)
                    weightNew[j] *= prof[k][nseqs[hold][j + k]];
            double old = startPos[hold];
            count = old == (startPos[hold] = sampler(weightNew))? count + 1: 0; /* simple criteria of convergence */
        }
    }
};

int main(int const, char const * const argv[])
{
    ifstream fasta(argv[1]);
    vector<string> seqs;
    for (string s; getline(fasta, s);)
        if ('>' == s.at(0))
            seqs.push_back("");
        else
            seqs.back() += s;
    int const mlen = atoi(argv[2]);
    Gsmf g(seqs, mlen);
    g.train(1000000);
    for (size_t i = 0; i < seqs.size(); ++ i)
        cerr << g.startPos[i] << "\t"<< seqs[i].substr(g.startPos[i], mlen) << endl;
    return 0;
}
