#ifndef QUERY_H
#define QUERY_H

#include <string>
#include <vector>

#include "Sequence.h"


#include <iostream>

class IndexItem
{
    public:
        IndexItem(bool strand, seqpos_t chr_pos, seqpos_t bases_count,
            uint64_t pointer)
        : strand_(strand), chr_pos_(chr_pos), bases_count_(bases_count),
        pointer_(pointer) {};
        
        bool get_strand();
        seqpos_t get_chr_pos();
        seqpos_t get_bases_count();
        uint64_t get_pointer();
        
    private:
        bool strand_;
        seqpos_t chr_pos_;
        seqpos_t bases_count_;
        uint64_t pointer_;
};

class BedQuery
{
    public:
        explicit BedQuery(std::string &bedline);
        BedQuery(const BedQuery &bq, std::string &chromosome, bool strand,
                 seqpos_t chr_size, seqpos_t start);
        BedQuery();
        ~BedQuery();
        std::string get_chr();
        seqpos_t get_start();
        seqpos_t get_end();
        std::string get_name();
        bool get_strand();
        seqpos_t get_thick_start();
        seqpos_t get_thick_end();
        unsigned get_exon_count();
        std::vector<seqpos_t>* get_exon_starts();
        std::vector<seqpos_t>* get_exon_ends();
        std::string get_bedline();
        bool merge_query(BedQuery *query, bool query_strand);
        bool merge_thick(BedQuery *query);
        bool merge_exons(std::vector<BedQuery*> &queries);
        void to_closed();
        void to_half_closed();
        
    private:
        std::string chromosome_, name_, rgb_;
        seqpos_t start_, end_, thick_start_, thick_end_, chr_size_;
        int score_;
        unsigned exon_count_;
        bool strand_, closed_, original_strand_;
        std::vector<seqpos_t>* exon_starts_;
        std::vector<seqpos_t>* exon_ends_;
        
        void set_numbers(std::string str, std::vector<seqpos_t>* numbers);
};

#endif /* QUERY_H */