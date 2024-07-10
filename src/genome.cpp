#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <string_view>

#include <Rcpp.h>
using namespace Rcpp;

using std::ifstream;
using std::string;
using std::vector;
using std::unordered_map;
using std::string_view;

class Genome {
public:
    Genome(string filename);
    vector<string> seqs() {
        return _seqs;
    };

    // input is 1-based
    char base_at(string const &chrom, int pos) {
        return _genome.at(chrom).at(pos - 1);
    };

    // input is 1-based
    string_view bases_at(string const &chrom, int start, int end) {
        // convert to 0-based
        int sv_start = start - 1;
        // length of the substring is calculated as end - start + 1
        int sv_len = end - sv_start;

        return string_view(_genome.at(chrom)).substr(sv_start, sv_len);
    };

    // given a genomic position and a motif with index of the base corresponding to the genomic position, check if the motif is present
    bool check_if_motif(string const &chrom, int pos, string const &motif, int motif_ind) {
        int motif_start = pos - motif_ind;
        int motif_end = motif_start + motif.size() - 1;
        string_view chrom_view = this->bases_at(chrom, motif_start, motif_end);
        return _genome.at(chrom).at(pos - 1) == motif.at(motif_ind);
    }
private:
    unordered_map<string, string> _genome;
    vector<string> _seqs;
};

// read from fasta file
Genome::Genome(string filename) {
    ifstream file(filename);

    if (!file.is_open()) {
        Rcpp::stop("Could not open file");
    }

    string line;
    string seq;
    string chrom;

    while (std::getline(file, line)) {
        if (line[0] == '>') {
            if (!seq.empty()) {
                _genome[chrom] = seq;
                _seqs.push_back(chrom);
            }

            chrom = line.substr(1);
            seq = "";
        } else {
            seq += line;
        }
    }

    if (!seq.empty()) {
        _genome[chrom] = seq;
        _seqs.push_back(chrom);
    }
}
