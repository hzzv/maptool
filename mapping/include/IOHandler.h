#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <istream>
#include <fstream>
#include <ostream>
#include <queue>

#include "../../../ocaml-bgzf/bgzf.h"

#include "Sequence.h"
#include "Query.h"

class IOHandler
{
    public:
        IOHandler(char header_fname[], char bin_fname[],
                  bool compressed, char maf_fname[]);
        ~IOHandler();
        void read_header(std::map<std::string, bioid_t> &genome_map,
                         std::vector< std::map<std::string,
                         std::pair <bioid_t, seqpos_t> > > &chr_maps,
                         std::map <bioid_t, std::vector<IndexItem*> > &index);
        void open_to_map();
        std::vector<Reference*>* read_references(std::vector<IndexItem*>
                                                 &index_items,
                                                 bioid_t ref_chr_id,
                                                 int indices[]);
        
    private:
        char header_fname_[1000];
        char bin_fname_[1000];
        bool compressed_;
        char maf_fname_[1000];
        bool map_, preprocess_, map_opened_;
        BGZF* bgzf_;
        std::ifstream ibin_;
        std::ofstream obin_;
        std::map <uint64_t, std::pair<Reference*, int> > cache_;
        static const int max_cache_size_ = 10;
        
        uint64_t bytes_to_number(std::istream &s, const int size);
        uint64_t read_bin_number(const int size);
        int max(int a, int b);
        std::vector<bool>* read_bin_sequence(seqpos_t length,
                                             std::vector<int>* rankselect =NULL,
                                             bool selecting = true);
        void add_to_cache(uint64_t pointer, Reference* reference);
        Reference* get_from_cache(uint64_t pointer);

};

#endif /* IOHANDLER_H */