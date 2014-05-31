#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <vector>
#include <map>
#include <cstdint>

#include <iostream>

typedef uint16_t bioid_t;
typedef int64_t seqpos_t;
typedef uint16_t biocount_t;

const int BIOID_SIZE = 2, SEQPOS_SIZE = 8, BIOCOUNT_SIZE = 2, STRAND_SIZE = 1,
    FILE_OFFSET_SIZE = 8, NAME_SIZE = 100, SELECT_BITS = 32, RANK_BITS = 32;
const int OLD_BIOID_SIZE1 = 1, OLD_BIOID_SIZE2 = 2, OLD_BIOCOUNT_SIZE1 = 1,
    OLD_SEQPOS_SIZE = 4, OLD_BIOCOUNT_SIZE2 = 4;


class Sequence
{
    public:
        Sequence(std::vector<bool>* sequence, bioid_t chr_id, seqpos_t chr_pos,
                 bool strand, seqpos_t bases_count)
        : sequence_(sequence), chr_id_(chr_id), chr_pos_(chr_pos),
        strand_(strand), bases_count_(bases_count)
        {};
        Sequence(std::vector<bool>* sequence, std::vector<int>* rankselect,
                 bioid_t chr_id, seqpos_t chr_pos, bool strand,
                 seqpos_t bases_count)
        : sequence_(sequence), rankselect_(rankselect), chr_id_(chr_id),
        chr_pos_(chr_pos), strand_(strand), bases_count_(bases_count)
        {
            has_rankselect_ = true;
        };
        
        virtual ~Sequence();
        
        void add_sequence(std::vector<bool>* sequence);
        std::vector<bool>* get_sequence();
        seqpos_t get_chr_pos();
        seqpos_t get_bases_count();
        seqpos_t get_chr_id();
        bool get_strand();
        seqpos_t length();
        virtual void print_info();
        virtual void print_seq();
        int select(int number);
        int rank(seqpos_t seq_pos);
        seqpos_t min(seqpos_t x, seqpos_t y);
        seqpos_t max(seqpos_t x, seqpos_t y);
        
        //TODO: implement or delete this
        char* to_bytes();
        
    private:
        std::vector<bool>* sequence_;
        std::vector<int>* rankselect_;
        bioid_t chr_id_;
        seqpos_t chr_pos_;
        bool strand_, has_rankselect_ = false;
        seqpos_t bases_count_;
};

class Reference;

class Informant: public Sequence
{
    public:
        Informant(std::vector<bool>* sequence, bioid_t chr_id, seqpos_t chr_pos,
                  bool strand, seqpos_t bases_count, seqpos_t seq_pos)
        : Sequence(sequence, chr_id, chr_pos, strand, bases_count)
        {
            seq_pos_ = seq_pos;
        }
        Informant(std::vector<bool>* sequence, bioid_t chr_id, seqpos_t chr_pos,
                  bool strand, seqpos_t bases_count, seqpos_t seq_pos,
                  Reference* aligned_to)
        : Sequence(sequence, chr_id, chr_pos, strand, bases_count)
        {
            seq_pos_ = seq_pos;
            aligned_to_ = aligned_to;
        }
        Informant(std::vector<bool>* sequence, std::vector<int>* rank,
                  bioid_t chr_id, seqpos_t chr_pos, bool strand,
                  seqpos_t bases_count, seqpos_t seq_pos, Reference* aligned_to)
        : Sequence(sequence, rank, chr_id, chr_pos, strand, bases_count)
        {
            seq_pos_ = seq_pos;
            aligned_to_ = aligned_to;
        }
        
        void print_info();
        seqpos_t get_seq_pos();
        Reference* get_ref();
        bool find_aligned_one(int way, int &jinf, int &jref);
        
        //TODO: implement or delete this
        char* to_bytes();
    private:
        seqpos_t seq_pos_;
        Reference* aligned_to_;
};

class Reference: public Sequence
{
    public:
        Reference(std::vector<bool>* sequence, bioid_t chr_id, seqpos_t chr_pos,
                 bool strand, seqpos_t bases_count)
        : Sequence(sequence, chr_id, chr_pos, strand, bases_count)
        {};
        
        Reference(std::vector<bool>* sequence, std::vector<int>* select,
                  bioid_t chr_id, seqpos_t chr_pos, bool strand,
                  seqpos_t bases_count)
        : Sequence(sequence, select, chr_id, chr_pos, strand, bases_count)
        {};
        
        ~Reference();
        
        std::vector<Informant*>* get_informant_vector(bioid_t inf_id);
        void add_informant(bioid_t inf_id, Informant* informant);
        void print_info();
        bool find_informant(/*std::vector<Informant*>::iterator &inf_it,*/
                            seqpos_t &inf_index,
                            bioid_t inf_id, seqpos_t seq_pos, int way);
        bool find_aligned_one(std::vector<Informant*>::iterator &inf_it,
                              bioid_t inf_id, seqpos_t seq_pos, int way,
                              seqpos_t &inf_seq_pos);
        
        //TODO: implement or delete this
        char* to_bytes();
        
    private:
        std::map<bioid_t, std::vector<Informant*> > informants_;
};

#endif /* SEQUENCE_H */