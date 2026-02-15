
#ifndef KACHE_HASH_REF_PARSER_HPP
#define KACHE_HASH_REF_PARSER_HPP



#include "zlib.h"
#include "kseq/kseq.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <queue>
#include <vector>
#include <iostream>


namespace kache_hash
{


// Declare the type of file handler and the read() function.
// Required for FASTA/FASTQ file reading using the kseq library.
KSEQ_INIT(gzFile, gzread)


// FASTA/FASTQ file parser.
class Parser
{
    typedef _KSEQ_DATA kseq_t;

private:

    std::queue<std::string> ref_paths;  // Collection of the reference file paths.
    gzFile file_ptr = nullptr;  // Pointer to the reference file being parsed.
    kseq_t* parser = nullptr;   // The kseq parser for the reference file being parsed.

    std::string curr_ref_path;  // Path to the reference currently being parsed.
    uint64_t ref_count = 0; // Number of the reference currently being parsed.
    uint64_t seq_id_; // Number of the current sequence (in the current reference).


    // Constructs a parser for the reference input collection `refs`.
    Parser(const std::vector<std::string>& refs);

    // Opens the reference at path `reference_path`.
    void open_reference(const std::string& reference_path);

    // Opens the next reference to be parsed from the collection `ref_file_paths`.
    bool open_next_reference();


public:

    // Constructs a parser for the file at path `file_path`.
    Parser(const std::string& file_path);

    // Returns the path to the reference currently being parsed.
    const std::string& curr_ref() const { return curr_ref_path; }

    // If sequences are remaining to be read, reads the next one
    // into memory and returns `true`. Returns `false` otherwise.
    bool read_next_seq();

    // Returns a pointer to the current sequence in the buffer.
    const char* seq() const { return parser->seq.s; }

    // Returns the length of the current sequence in the buffer.
    size_t seq_len() const { return parser->seq.l; }

    // Returns the current size of the buffer.
    size_t buff_sz() const { return parser->seq.m; }

    // Returns the id (number) of the current reference being parsed.
    uint64_t ref_id() const { return ref_count; }

    // Returns the id (number) of the current sequence in the buffer.
    uint64_t seq_id() const { return seq_id_; }

    // Returns the name (as parsed) of the current sequence in the buffer.
    const char* seq_name() const { return parser->name.s; }

    // Closes the internal kseq parser for the current reference.
    void close();
};


inline Parser::Parser(const std::string& file_path)
{
    ref_paths.push(file_path);

    // Open the first reference for subsequent parsing.
    open_next_reference();
}


inline Parser::Parser(const std::vector<std::string>& refs):
    ref_paths(std::deque<std::string>(refs.begin(), refs.end()))
{}


inline void Parser::open_reference(const std::string& reference_path)
{
    file_ptr = gzopen(reference_path.c_str(), "r");  // Open the file handler.
    if(file_ptr == nullptr)
    {
        std::cerr << "Error opening reference file " << reference_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    parser = kseq_init(file_ptr);   // Initialize the kseq parser.

    curr_ref_path = reference_path;
    ref_count++;
    seq_id_ = 0;
}


inline bool Parser::open_next_reference()
{
    if(ref_paths.empty())
        return false;


    open_reference(ref_paths.front());
    ref_paths.pop();

    return true;
}


inline bool Parser::read_next_seq()
{
    // Sequences still remain at the current reference being parsed.
    if(parser != nullptr && kseq_read(parser) >= 0)
    {
        seq_id_++;
        return true;
    }

    // The current reference has been parsed completely. Close its handles.
    close();

    // Start parsing the next reference, if exists.
    if(open_next_reference())
        return read_next_seq();

    return false;
}


inline void Parser::close()
{
    if(file_ptr != nullptr)
    {
        kseq_destroy(parser);   // Close the kseq parser.
        gzclose(file_ptr);  // Close the file handler.

        parser = nullptr;
        file_ptr = nullptr;
    }
}

}



#endif
