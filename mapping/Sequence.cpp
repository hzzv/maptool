#include <vector>
#include <map>
#include <cstdint>

#include <iostream>
#include <typeinfo>

#include "include/Sequence.h"

using std::vector;
using std::map;

Sequence::~Sequence()
{
    delete sequence_;
    if (has_rankselect_) delete rankselect_;
}

void Sequence::add_sequence(std::vector<bool>* sequence)
{
    for (auto it = sequence->begin(); it != sequence->end(); ++it)
    {
        sequence_->push_back(*it);
    }
}

std::vector<bool>* Sequence::get_sequence()
{
    return sequence_;
}

seqpos_t Sequence::get_chr_pos()
{
    return chr_pos_;
}

seqpos_t Sequence::get_bases_count()
{
    return bases_count_;
}

seqpos_t Sequence::get_chr_id()
{
    return chr_id_;
}

bool Sequence::get_strand()
{
    return strand_;
}

seqpos_t Sequence::length()
{
    return sequence_->size();
}

void Sequence::print_info()
{
    std::cout << chr_id_ << "\t" << chr_pos_ << "\t" << strand_ << "\t" <<\
        bases_count_ << std::endl;
}

void Sequence::print_seq()
{
    for (auto it = sequence_->begin(); it != sequence_->end(); ++it)
        std::cout << *it;
    std::cout << std::endl;
}

// Find index of 'number'-th '1' in this sequence
int Sequence::select(int number)
{
    int seq_pos = (*rankselect_)[number/SELECT_BITS];
    number = number % SELECT_BITS;
    while ((seq_pos+1 < (int)sequence_->size()) &&
           ((number > 0) || (!(*sequence_)[seq_pos+1])))
    {
        number -= (*sequence_)[++seq_pos];
    }
    return seq_pos+1;
}


int Sequence::rank(seqpos_t seq_pos)
{
    // Uncomment the commented lines to put rank in use (do not forget about
    // commented lines in IOHandler.cpp)
//     int ret = this->get_chr_pos() + (*rankselect_)[seq_pos/RANK_BITS];
//     int from = seq_pos - (seq_pos % RANK_BITS);
//     for (int i = 0; i < (seq_pos % RANK_BITS); ++i)
//     {
//         if ((*(this->get_sequence()))[from+i]) ++ret;
//     }
//     return ret;
    seqpos_t ret = this->get_chr_pos();
    for (int i = 0; i < seq_pos; ++i)
    {
        if ((*(this->get_sequence()))[i]) ++ret;
    }
    return ret;
}

seqpos_t Sequence::min(seqpos_t x, seqpos_t y)
{
    if (x < y) return x;
    return y;
}

seqpos_t Sequence::max(seqpos_t x, seqpos_t y)
{
    if (x > y) return x;
    return y;
}

void Informant::print_info()
{
    std::cout << this->get_bases_count() << " "  << seq_pos_ << std::endl;
}

seqpos_t Informant::get_seq_pos()
{
    return seq_pos_;
}

Reference* Informant::get_ref()
{
    return aligned_to_;
}

// Find '1' in inf. sequence corresponding to given '1' in ref. if possible
bool Informant::find_aligned_one(int way, int &jinf, int &jref)
{
    while (!((*get_sequence())[jinf]) ||
           !((*(aligned_to_->get_sequence()))[jref]))
    {
        jinf += way;
        jref += way;
        if (jinf >= length()) return false;
        if (jinf < 0) return false;
    }
    return true;
}

Reference::~Reference()
{
    for (map<bioid_t, vector<Informant*> >::iterator it = informants_.begin();
         it != informants_.end(); ++it)
    {
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            delete (*it2);
        }
    }
}

std::vector<Informant*>* Reference::get_informant_vector(bioid_t inf_id)
{
    return &(informants_[inf_id]);
}

void Reference::add_informant(bioid_t inf_id, Informant* informant)
{
    if (informants_.count(inf_id) == 0)
    {
        informants_[inf_id] = vector<Informant*>();
    }
    informants_[inf_id].push_back(informant);
}

// Finds the informant which is aligned to 'seq_pos'
bool Reference::find_informant(seqpos_t &inf_index,
                               bioid_t inf_id, seqpos_t seq_pos, int way)
{
    unsigned lo = 0, hi = informants_[inf_id].size(), mid;
    if ((hi == 0) ||
        ((seq_pos < informants_[inf_id][0]->get_seq_pos()) && (way == -1)) ||
        ((seq_pos >= informants_[inf_id][hi-1]->get_seq_pos() +
         informants_[inf_id][hi-1]->length()) && (way == 1)))
    {
        return false;
    }
    else if ((seq_pos < informants_[inf_id][0]->get_seq_pos() +
              informants_[inf_id][0]->length()) && (way == 1))
    {
        inf_index = 0;
        return true;
    }
    else if ((seq_pos >= informants_[inf_id][hi-1]->get_seq_pos()) &&
             (way == -1))
    {
        inf_index = informants_[inf_id].size() - 1;
        return true;
    }
    // Now is guaranteed that seq_pos is contained in this interval:
    // [first informant, last informant]
    // Binary search seq_pos
    while (lo < hi)
    {
        mid = (lo+hi)/2;
        if ((seq_pos - informants_[inf_id][mid]->get_seq_pos() <
             informants_[inf_id][mid]->length()) &&
            (seq_pos - informants_[inf_id][mid]->get_seq_pos() >= 0))
        {
            hi = lo;
        }
        if ((mid + 1 < informants_[inf_id].size()) &&
            (seq_pos >= informants_[inf_id][mid]->get_seq_pos() +
             informants_[inf_id][mid]->length()) &&
            (seq_pos < informants_[inf_id][mid+1]->get_seq_pos()))
        {
            if (way == 1)
            {
                ++mid;
                hi = lo;
            }
            else hi = lo;
        }
        else if ((seq_pos >= informants_[inf_id][mid]->get_seq_pos() +
                  informants_[inf_id][mid]->length()) &&
                 (seq_pos - informants_[inf_id][mid]->get_seq_pos() >= 0))
        {
            lo = mid;
        }
        else if ((seq_pos - informants_[inf_id][mid]->get_seq_pos() <
                 informants_[inf_id][mid]->length()) &&
                 (seq_pos - informants_[inf_id][mid]->get_seq_pos() < 0))
        {
            hi = mid;
        }
        else hi = lo;
    }
    inf_index = mid;
    return true;
}

void Reference::print_info()
{
    Sequence::print_info();
    std::cout << "Informants:" << std::endl;
    for (auto it = informants_.begin(); it != informants_.end(); ++it)
    {
        std::cout << "Inf id: " << it->first << " count " << it->second.size() << std::endl;
        for (vector<Informant*>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            (*it2)->Informant::print_info();
        }
    }
}