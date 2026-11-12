
#ifndef KACHE_HASH_EXAMPLES_FASTX_HPP
#define KACHE_HASH_EXAMPLES_FASTX_HPP



// A small reader of uncompressed FASTA and FASTQ files for the examples, which
// reads a file whole, as records of names and sequences; and a splitter of
// sequences into pieces for threads to share. kache-hash's counter library has
// a parallel reader, `Fastx_Reader`, which reads files and pipes in chunks.

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>


namespace examples
{

// A record of a FASTA or FASTQ file.
struct Record
{
    std::string name;   // Its header line, without the leading `>` or `@`.
    std::string seq;    // Its sequence: a FASTA record's lines, joined.
};


// Returns the records of the uncompressed FASTA or FASTQ file at `path`. The
// first character sets the format: '>' for FASTA, '@' for FASTQ, whose records
// have four lines each.
inline std::vector<Record> read_fastx(const std::string& path)
{
    std::ifstream in(path);
    if(!in)
        throw std::runtime_error("cannot open " + path);

    const auto next = [&](std::string& line)
    {
        if(!std::getline(in, line))
            return false;
        if(!line.empty() && line.back() == '\r')
            line.pop_back();
        return true;
    };

    std::vector<Record> recs;
    std::string line;
    if(in.peek() == '@')
        while(next(line))
        {
            recs.push_back({line.substr(1), ""});
            if(!next(recs.back().seq) || !next(line) || !next(line))    // The sequence, `+`, and the qualities.
                throw std::runtime_error("a truncated FASTQ record in " + path);
        }
    else
        while(next(line))
            if(!line.empty() && line[0] == '>')
                recs.push_back({line.substr(1), ""});
            else if(!recs.empty())
                recs.back().seq += line;

    return recs;
}


// Returns the pieces of the sequences of the records `recs` of at most `piece`
// k-mer positions each, overlapping by `k - 1` characters, so that the pieces'
// k-mer positions are the sequences': threads can share even one genome.
inline std::vector<std::string_view> pieces(const std::vector<Record>& recs, const std::size_t k, const std::size_t piece = std::size_t(1) << 14)
{
    std::vector<std::string_view> p;
    for(const auto& r : recs)
        for(std::size_t i = 0; i + k <= r.seq.size(); i += piece)
            p.emplace_back(r.seq.data() + i, std::min(piece + k - 1, r.seq.size() - i));

    return p;
}


// Returns the number of k-mer positions of the sequences of the records
// `recs`: a bound on their distinct k-mers, to size a table for.
inline std::size_t positions(const std::vector<Record>& recs, const std::size_t k)
{
    std::size_t n = 0;
    for(const auto& r : recs)
        n += r.seq.size() >= k ? r.seq.size() - k + 1 : 0;

    return n;
}

}



#endif
