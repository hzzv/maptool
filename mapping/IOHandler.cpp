#include <istream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <ios>
#include <ostream>
#include <stdexcept>
#include <cstring>

#include <iomanip>

#include "../../ocaml-bgzf/bgzf.h"

#include "include/Sequence.h"
#include "include/Query.h"
#include "include/IOHandler.h"

using std::istream;
using std::ifstream;
using std::map;
using std::string;
using std::vector;
using std::pair;
using std::make_pair;


IOHandler::IOHandler(char header_fname[], char bin_fname[], bool compressed,
                     char maf_fname[]):
    compressed_(compressed)
{
    strcpy(header_fname_, header_fname);
    strcpy(bin_fname_, bin_fname);
    strcpy(maf_fname_, maf_fname);
    map_ = false;
    preprocess_ = false;
    if (maf_fname_[0] == '\0') map_ = true;
    else preprocess_ = true;
    map_opened_ = false;
}

IOHandler::~IOHandler()
{
    if (map_opened_)
    {
        if (compressed_) bgzf_close(bgzf_);
        else ibin_.close();
    }
    for (auto it = cache_.begin(); it != cache_.end(); ++it)
    {
        delete it->second.first;
    }
}

// Opens the BGZF or BIN file with preprocessed alignments or an empty file
void IOHandler::open_to_map()
{
    try
    {
        if (map_)
        {
            if (compressed_) bgzf_ = bgzf_open(bin_fname_, "r");
            else ibin_.open(bin_fname_, std::ios::in | std::ios::binary);
        }
        if (preprocess_)
        {
            if (compressed_) bgzf_ = bgzf_open(bin_fname_, "w");
            else obin_.open(bin_fname_, std::ios::out | std::ios::binary);
        }
    }
    catch (std::exception &e)
    {
        std::cerr << "Error opening " << bin_fname_ << "." << std::endl;
        exit(1);
    }
    map_opened_ = true;
}

// Read 'size' bytes and returns number represented by those bytes
// 'size' must be <= 8
uint64_t IOHandler::bytes_to_number(istream &s, const int size)
{
    uint64_t number = 0;
    char data[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    s.read(data, size);
    for (int i = 0; i < size; ++i)
    {
        number <<= 8;
        number += (unsigned uint8_t)data[i];
    }
    return number;
}

// Read information from header into given structures
void IOHandler::read_header(map <string, bioid_t> &genome_map,
    vector <map <string, pair <bioid_t, seqpos_t> > > &chr_maps,
    map <bioid_t, vector <IndexItem*> > &index)
{
    ifstream s(header_fname_, std::ios::in | std::ios::binary);
    s.unsetf(std::ios::skipws);
    s.exceptions(istream::failbit | istream::badbit);
    
    // Read into 'genome_map': name and id
    biocount_t genome_count = bytes_to_number(s, OLD_BIOCOUNT_SIZE1);
    for (int i = 0; i < genome_count; ++i)
    {
        bioid_t name_len = bytes_to_number(s, OLD_BIOID_SIZE1);
        string name;
        name.resize(name_len);
        s.read(&name[0], name_len);        
        bioid_t genome_id = bytes_to_number(s, OLD_BIOID_SIZE1);
        genome_map[name] = genome_id;
    }
    
    // Read into 'chr_maps': for each genome and chromosome read
    // chromosome name and chromosome id
    for (unsigned int i = 0; i < genome_map.size(); ++i)
    {
        bioid_t chr_count = bytes_to_number(s, OLD_BIOID_SIZE2);
        map <string, pair <bioid_t, seqpos_t> > chr_map;
        for (int j = 0; j < chr_count; ++j)
        {
            bioid_t name_len = bytes_to_number(s, OLD_BIOID_SIZE1);
            string name;
            name.resize(name_len);
            s.read(&name[0], name_len);            
            bioid_t chr_id = bytes_to_number(s, OLD_BIOID_SIZE2);
            seqpos_t chr_len = bytes_to_number(s, OLD_SEQPOS_SIZE);
            chr_map[name] = make_pair(chr_id, chr_len);
        }
        chr_maps.push_back(chr_map);
    }
    
    // Read into 'index': for each chromosome in reference and each
    // reference block read position on chromosome, count of bases
    // and pointer to bgzf file
    for (unsigned int i = 0; i < chr_maps[0].size(); ++i)
    {
        bioid_t chr_id = bytes_to_number(s, OLD_BIOID_SIZE2);
        bioid_t ref_count = bytes_to_number(s, OLD_BIOCOUNT_SIZE2);
        vector<IndexItem*> chr_index;
        for (int j = 0; j < ref_count; ++j)
        {
            bool strand = bytes_to_number(s, STRAND_SIZE);
            seqpos_t chr_pos = bytes_to_number(s, OLD_SEQPOS_SIZE);
            seqpos_t bases_count = bytes_to_number(s, OLD_SEQPOS_SIZE);
            uint64_t pointer = bytes_to_number(s, FILE_OFFSET_SIZE);
            chr_index.push_back(new IndexItem(strand, chr_pos, bases_count,
                pointer));
        }
        index[chr_id] = chr_index;
    }
    
    s.close();
}

// Read from BGZF/BIN file 'size' bytes to a number
uint64_t IOHandler::read_bin_number(const int size)
{
    if (!map_) return 0;
    uint64_t number = 0;
    char data[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    if (compressed_)
    {
        if (bgzf_read(bgzf_, data, size) != size)
        {
            throw std::runtime_error("Unreadable BGZF file" +
                                     string(bin_fname_));
        }
    }
    else
    {
        ibin_.read(data, size);
        if (ibin_.gcount() != size)
        {
            throw std::runtime_error("Unreadable binary file" +
                                     string(bin_fname_));
        }
    }
    for (int i = 0; i < size; ++i)
    {
        number <<= 8;
        number += (unsigned uint8_t)data[i];
    }
    return number;
}

int IOHandler::max(int a, int b)
{
    if (a > b) return a;
    return b;
}

// Read from BGZF/BIN file binary sequence of length 'length'
vector<bool>* IOHandler::read_bin_sequence(seqpos_t length,
                                           vector<int>* rankselect,
                                           bool selecting)
//TODO: this needs to be changed if the format of data in BGZF file will change
{
    if (!map_) return new vector<bool>;
    seqpos_t real_length = length/8;
    if (length % 8 != 0) real_length += 1;
    // Read sequence
    char* data = new char[real_length];
    if (compressed_)
    {
        if (bgzf_read(bgzf_, data, real_length) != real_length)
        {
            delete[] data;
            throw std::runtime_error("Unreadable BGZF file" +
                                     string(bin_fname_));
        }
    }
    else
    {
        ibin_.read(data, real_length);
        if (ibin_.gcount() != real_length)
        {
            delete[] data;
            throw std::runtime_error("Unreadable binary file" +
                                     string(bin_fname_));
        }
    }
    // Convert read data to vector of booleans, optionally fill rank/select
    vector<bool>* ret = new vector<bool>;
    ret->reserve(length);
    if (rankselect != NULL)
    {
        if (selecting) rankselect->push_back(-1);
    }
    int one_bits = 0;
    for (unsigned int i = 0; i < real_length; ++i)
    {
        if ((rankselect != NULL) && (!selecting) && ((i*8) % RANK_BITS == 0))
            rankselect->push_back(one_bits);
        for (unsigned int j = max(0, 8 - length + i*8); j < 8; ++j)
        {
            ret->push_back((bool)((data[i] >> (7-j)) & 1));
            one_bits += ((data[i] >> (7-j)) & 1);
            if ((rankselect != NULL) && selecting && (one_bits == SELECT_BITS))
            {
                one_bits = 0;
                rankselect->push_back(ret->size()-1);
            }
        }
    }
    delete[] data;
    return ret;
}

// Read from BGZF/BIN file references on given indices
vector<Reference*>* IOHandler::read_references(vector<IndexItem*> &index_items,
                                               bioid_t ref_chr_id,
                                               int indices[])
{
    vector<Reference*>* references = new vector<Reference*>;
    if (!map_opened_) return references;
    seqpos_t length, chr_pos, seq_pos, seq_len, bases_count;
    biocount_t inf_number, inf_block_num;
    bioid_t inf_id, chr_id;
    bool strand;
    for (int i = indices[0]; i <= indices[1]; ++i)
    {
        Reference* cached = get_from_cache(index_items[i]->get_pointer());
        if (cached != NULL)
        {
            references->push_back(cached);
            continue;
        }
        // Seek in BGZF/BIN file
        if (compressed_)
        {
            if (bgzf_seek(bgzf_, index_items[i]->get_pointer(), SEEK_SET) == -1)
            {
                throw std::runtime_error("Unreadable binary file" +
                                        string(bin_fname_));
            }
        }
        else
        {
            ibin_.seekg((int)index_items[i]->get_pointer());
            if (ibin_.fail())
            {
                throw std::runtime_error("Unreadable binary file" +
                                        string(bin_fname_));
            }
        }
        // Read reference information
        length = read_bin_number(OLD_SEQPOS_SIZE);
        vector<int>* select = new vector<int>;
        vector<bool>* sequence = read_bin_sequence(length, select, true);
        references->push_back(new Reference(sequence, select, ref_chr_id,
                              index_items[i]->get_chr_pos(),
                              index_items[i]->get_strand(),
                              index_items[i]->get_bases_count()));
        // Read reference's informant information
        inf_number = read_bin_number(OLD_BIOCOUNT_SIZE1);
        vector< pair<bioid_t, biocount_t> > infs;
        for (int j = 0; j < inf_number; ++j)
        {
            inf_id = read_bin_number(OLD_BIOID_SIZE1);
            inf_block_num = read_bin_number(OLD_BIOCOUNT_SIZE2);
            infs.push_back(make_pair(inf_id, inf_block_num));
        }
        for (auto it = infs.begin(); it != infs.end(); ++it)
        {
            for (int k = 0; k < it->second; ++k)
            {
                chr_id = read_bin_number(OLD_BIOID_SIZE2);
                strand = read_bin_number(STRAND_SIZE);
                chr_pos = read_bin_number(OLD_SEQPOS_SIZE);
                seq_pos = read_bin_number(OLD_SEQPOS_SIZE) - 1;
                seq_len = read_bin_number(OLD_SEQPOS_SIZE);
                bases_count = read_bin_number(OLD_SEQPOS_SIZE);
                // Uncomment the commented lines to compute rank for informant
                // (do not forget about commented lines in Sequence.cpp)
//                 vector<int>* rank = new vector<int>;
                vector<bool>* sequence = read_bin_sequence(seq_len);
//                     read_bin_sequence(seq_len, rank, false);
                (*references)[references->size()-1]->add_informant(
                    it->first,
                    new Informant(sequence, /*rank,*/ chr_id, chr_pos, strand,
                                  bases_count, seq_pos,
                                  (*references)[references->size()-1]));
            }
        }
    }
    for (int i = indices[0]; i <= indices[1]; ++i)
    {
        add_to_cache(index_items[i]->get_pointer(),
                     (*references)[i-indices[0]]);
    }
    return references;
}

void IOHandler::add_to_cache(uint64_t pointer, Reference* reference)
{
    if ((cache_.count(pointer) > 0) && (cache_[pointer].second == 0)) return;
    for (auto it = cache_.begin(); it != cache_.end(); ++it)
    {
        ++it->second.second;
    }
    if (cache_.count(pointer) > 0) cache_[pointer].second = 0;
    else
    {
        if (cache_.size() == max_cache_size_)
        {
            int max = 0;
            uint64_t to_erase = 0;
            for (auto it = cache_.begin(); it != cache_.end(); ++it)
            {
                if (it->second.second > max)
                {
                    max = it->second.second;
                    to_erase = it->first;
                }
            }
            delete cache_[to_erase].first;
            cache_.erase(to_erase);
        }
        cache_.insert(make_pair(pointer, make_pair(reference, 0)));
    }
}

Reference* IOHandler::get_from_cache(uint64_t pointer)
{
    if (cache_.count(pointer) > 0)
    {
        return cache_[pointer].first;
    }
    else return NULL;
}