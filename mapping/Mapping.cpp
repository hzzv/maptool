#include <map>
#include <string>
#include <vector>
#include <utility>

#include "../../ocaml-bgzf/bgzf.h"

#include "include/Query.h"
#include "include/Sequence.h"
#include "include/Mapping.h"


#include <iostream>

using std::string;
using std::map;
using std::vector;
using std::pair;
using std::make_pair;


Mapping::Mapping(IOHandler* ioh, string &informant, int inf_maxgap,
                 int ref_maxgap, bool inner, bool alwaysmap,
                 map < string, bioid_t > *genome_map,
                 vector < map < string, pair <bioid_t, seqpos_t> > > *chr_maps,
                 map <bioid_t, vector<IndexItem*> > *index):
    ioh_(ioh), informant_(informant), inf_maxgap_(inf_maxgap),
    ref_maxgap_(ref_maxgap), inner_(inner), alwaysmap_(alwaysmap),
    genome_map_(genome_map), chr_maps_(chr_maps), index_(index)
{
    inf_id_ = (*genome_map)[informant];
    id_to_len_ = new map <bioid_t, pair <string, seqpos_t> >();
    for (map <string, pair <bioid_t, seqpos_t> >::iterator it =\
        (*chr_maps)[inf_id_].begin(); it != (*chr_maps)[inf_id_].end(); ++it)
    {
        (*id_to_len_)[(it->second).first] = make_pair(it->first,
                                                      (it->second).second);
    }
    query_ = NULL;
}

Mapping::~Mapping()
{
    delete this->id_to_len_;
}

string Mapping::get_error_message(string error_name)
{
    string to_add = "";
    for (int i = 0; i < known_error_count_; ++i)
    {
        if (error_name.compare(known_error_names_[i]) == 0)
        {
            if ((error_name.compare("inf_gap") == 0) ||
                (error_name.compare("ref_gap") == 0))
                to_add = std::to_string(found_gap_);
            return known_error_messages_[i] + to_add;
        }
    }
    return "";
}

void Mapping::print_errors()
{
    std::cerr << query_->get_name() << " ";
    for (auto it = errors_.begin(); it != errors_.end(); ++it)
    {
        std::cerr << *it << " " << get_error_message(*it) << std::endl;
    }
    errors_.clear();
}

// Set one error and throw exception
void Mapping::error(string error_name)
{
    errors_.push_back(error_name);
    throw MappingError();
}

void Mapping::set_query(BedQuery* query)
{
    query_ = query;
}

seqpos_t Mapping::min(seqpos_t x, seqpos_t y)
{
    if (x < y) return x;
    return y;
}

seqpos_t Mapping::max(seqpos_t x, seqpos_t y)
{
    if (x > y) return x;
    return y;
}

// Return references containing given positions
vector<Reference*>* Mapping::get_references(seqpos_t start, seqpos_t end)
{
    seqpos_t lo = 0, hi = (*index_)[ref_chr_id_].size(), mid;
    int indices[2] = {-1, -1};
    seqpos_t queries[2] = {start, end};
    for (int i = 0; i < 2; ++i)
    {
        while (hi > lo)
        {
            mid = (lo+hi)/2;
            if ((*index_)[ref_chr_id_][mid]->get_chr_pos() <= queries[i] &&\
                (*index_)[ref_chr_id_][mid]->get_chr_pos() +\
                (*index_)[ref_chr_id_][mid]->get_bases_count() > queries[i])
            {
                indices[i] = mid;
                hi = lo;
            }
            else if ((*index_)[ref_chr_id_][mid]->get_chr_pos() > queries[i])
                hi = mid;
            else if (lo != mid) lo = mid;
            else hi = lo;
        }
        if (indices[i] == -1) error("no_mapping");
        lo = indices[i];
        hi = (*index_)[ref_chr_id_].size();
    }
    return ioh_->read_references((*index_)[ref_chr_id_], ref_chr_id_, indices);
}

// Fill 'informants' by all informants belonging to references in 'references'
void Mapping::fill_informant_vector(vector<Reference*> &references,
                                    vector<Informant*> &informants)
{
    for (auto ref_it = references.begin(); ref_it != references.end(); ++ref_it)
    {
        for (auto inf_it = (*ref_it)->get_informant_vector(inf_id_)->begin();
            inf_it != (*ref_it)->get_informant_vector(inf_id_)->end(); ++inf_it)
        {
            informants.push_back(*inf_it);
        }
    }
}

// Get number of informants belonging to references from 'from' to 'to'
seqpos_t Mapping::get_inf_count(vector<Reference*>::iterator &from,
                                const vector<Reference*>::iterator &to)
{
    seqpos_t overall_inf_index = 0;
    while (from != to)
    {
        overall_inf_index += (*from)->get_informant_vector(inf_id_)->size();
        ++from;
    }
    return overall_inf_index;
}

// Map one position from reference to informant
BedQuery* Mapping::map_position(vector<Reference*> &references,
                                vector<Informant*> &informants,
                                seqpos_t position, int way,
                                vector<Reference*>::iterator &ref_it,
                                vector<Informant*>::iterator &inf_it_ret)
{
    position -= (*ref_it)->get_chr_pos();
    // Find index of position-th '1' in references
    seqpos_t seq_pos = (*ref_it)->select(position);
    // Find instance of Informant in which is position corresponding to seq_pos
    seqpos_t inf_index;
    int gap = 0;
    while (!((*ref_it)->find_informant(inf_index, inf_id_, seq_pos, way)))
    {
        if (way == 1)
        {
            gap += (*ref_it)->length() - seq_pos;
            ++ref_it;
            if (ref_it == references.end()) error("pos_to_gap");
            seq_pos = 0;
        }
        else
        {
            if (ref_it == references.begin()) error("pos_to_gap");
            gap += seq_pos;
            --ref_it;
            seq_pos = (*ref_it)->length() - 1;
        }
        if (gap > ref_maxgap_) error("ref_gap");
    }
    auto ref_begin = references.begin();
    inf_index += get_inf_count(ref_begin, ref_it);
    vector<Informant*>::iterator inf_it = informants.begin() + inf_index;
    // Find '1' in inf. sequence corresponding to given '1' in ref.
    int jinf;
    if (way == 1) jinf = max(0, seq_pos - (*inf_it)->get_seq_pos());
    else jinf = min((*inf_it)->length() - 1,
                    seq_pos - (*inf_it)->get_seq_pos());
    int jref = (*inf_it)->get_seq_pos() + jinf;
    Reference* ref = (*inf_it)->get_ref();
    while (!((*inf_it)->find_aligned_one(way, jinf, jref)))
    {
        if (jinf >= (*inf_it)->length())
        {
            if (inf_it != informants.end()) ++inf_it;
            else error("pos_to_gap");
            if (ref->get_chr_pos() != ((*inf_it)->get_ref())->get_chr_pos())
            {
                ref = (*inf_it)->get_ref();
                jref = 0;
            }
            jinf = 0;
        }
        if (jinf < 0)
        {
            if (inf_it == informants.begin()) error("pos_to_gap");
            --inf_it;
            if (ref->get_chr_pos() != ((*inf_it)->get_ref())->get_chr_pos())
            {
                ref = (*inf_it)->get_ref();
                jref = ref->length() - 1;
            }
            jinf = (*inf_it)->length() - 1;
        }
    }
    // Find order of the '1' on jinf-th positon
    seqpos_t inf_pos = (*inf_it)->rank(jinf);
    inf_it_ret = inf_it;
    BedQuery* ret = new BedQuery(*query_,
                                 (*id_to_len_)[(*inf_it)->get_chr_id()].first,
                                 (*inf_it)->get_strand(),
                                 (*id_to_len_)[(*inf_it)->get_chr_id()].second,
                                 inf_pos);
    return ret;
}

// Check if given informants are correctly preceeding each other
bool Mapping::check_informants(vector<Informant*>::iterator &inf_it1,
                               vector<Informant*>::iterator &inf_it2)
{
    if (inf_it1 > inf_it2)
    {
        errors_.push_back("inf_preceed");
        return false;
    }
    seqpos_t last_inf_end, last_ref_end;
    bool last_strand;
    bioid_t last_chr_id;
    while (inf_it1 != inf_it2)
    {
        last_inf_end = (*inf_it1)->get_chr_pos() +\
            (*inf_it1)->get_bases_count();
        last_ref_end = (*inf_it1)->get_seq_pos() +\
            (*inf_it1)->get_bases_count();
        last_strand = (*inf_it1)->get_strand();
        last_chr_id = (*inf_it1)->get_chr_id();
        ++inf_it1;
        
        if ((last_inf_end > (*inf_it1)->get_chr_pos()) ||
            (last_strand != (*inf_it1)->get_strand()) ||
            (last_chr_id != (*inf_it1)->get_chr_id()) ||
            ((inf_maxgap_ > -1) &&
             ((*inf_it1)->get_chr_pos() - last_inf_end > inf_maxgap_)) ||
            ((ref_maxgap_ > -1) &&
             ((*inf_it1)->get_seq_pos() - last_ref_end > ref_maxgap_)))
        {
            if (last_inf_end > (*inf_it1)->get_chr_pos())
                errors_.push_back("inf_preceed");
            if (last_strand != (*inf_it1)->get_strand())
                errors_.push_back("inf_strand");
            if (last_chr_id != (*inf_it1)->get_chr_id())
                errors_.push_back("inf_contig");
            if ((inf_maxgap_ > -1) &&
                ((*inf_it1)->get_chr_pos() - last_inf_end > inf_maxgap_))
            {
                errors_.push_back("inf_gap");
                found_gap_ = (*inf_it1)->get_chr_pos() - last_inf_end;
            }
            if ((ref_maxgap_ > -1) &&
                ((*inf_it1)->get_seq_pos() - last_ref_end > ref_maxgap_))
            {
                errors_.push_back("ref_gap");
                found_gap_ = (*inf_it1)->get_seq_pos() - last_ref_end;
            }
            return false;
        }
    }
    return true;
}

// Get a mapping of the given interval
BedQuery* Mapping::get_mapping(seqpos_t start, seqpos_t end,
                               vector <Reference*> &references,
                               vector<Informant*> &informants,
                               vector<Reference*>::iterator &ref_it1,
                               vector<Reference*>::iterator &ref_it2)
{    
    if (start > end) error("invalid_query");
    // Get mapping of two positions
    vector<Informant*>::iterator inf_it1;
    vector<Informant*>::iterator inf_it2;
    int way = 1;
    if (!inner_) way = -1;
    BedQuery *answer1 = map_position(references, informants, start, way,
                                     ref_it1, inf_it1);
    BedQuery *answer2 = map_position(references, informants, end, (-1) * way,
                                     ref_it2, inf_it2);
    
    // Check if the positions make up an interval
    if (!(answer1->merge_query(answer2, query_->get_strand())))
        error("no_mapping");
    if (!check_informants(inf_it1, inf_it2))
    {
        throw MappingError();
    }
    
    delete answer2;
    return answer1;
}

// Sets given iterators such that corresponding Reference-s contain start, end
void Mapping::set_ref_iterators(seqpos_t start, seqpos_t end,
                                vector<Reference*> &references,
                                vector<Reference*>::iterator &ref_it1,
                                vector<Reference*>::iterator &ref_it2)
{
    while (((*ref_it1)->get_chr_pos() +
            (*ref_it1)->get_bases_count() <= start) ||
           ((*ref_it2)->get_chr_pos() > end))
    {
        if ((*ref_it1)->get_chr_pos() + (*ref_it1)->get_bases_count() <= start)
            ++ref_it1;
        if ((*ref_it2)->get_chr_pos() > start)
            --ref_it2;
        if ((ref_it1 == references.end()) || ((ref_it2 == references.begin()) &&
             ((*ref_it2)->get_chr_pos() > start)) || (ref_it1 > ref_it2))
        {
            error("no_thick_mapping");
        }
    }  
}

// Get mapping of a given BED line - interval, thick interval and exons
BedQuery* Mapping::get_answer()
{
    if ((query_ == NULL) || (query_->get_start() > query_->get_end()))
        error("invalid_query");
    // Get references
    ref_chr_id_ = (*chr_maps_)[0][query_->get_chr()].first;
    vector <Reference*> *references = get_references(query_->get_start(),
                                                     query_->get_end());
    if (references->size() == 0) error("no_mapping");
    vector<Reference*>::iterator ref_it1 = references->begin();
    vector<Reference*>::iterator ref_it2 = --(references->end());
    vector<Informant*> informants;
    fill_informant_vector(*references, informants);
    // Interval
    int imaxgap = inf_maxgap_;
    if (query_->get_exon_count() > 0) inf_maxgap_ = -1;
    BedQuery* answer = get_mapping(query_->get_start(), query_->get_end(),
                                   *references, informants, ref_it1, ref_it2);
    
    // Thick interval
    if (query_->get_thick_start() != -1)
    {
        BedQuery* thick_answer;
        if ((query_->get_thick_start() == query_->get_start()) &&
            (query_->get_thick_end() == query_->get_end()))
        {
            thick_answer = answer;
        }
        else
        {
            set_ref_iterators(query_->get_thick_start(),
                              query_->get_thick_end(), *references,
                              ref_it1, ref_it2);
            thick_answer = get_mapping(query_->get_thick_start(),
                                       query_->get_thick_end(),
                                       *references, informants, ref_it1,
                                       ref_it2);
        }
        answer->merge_thick(thick_answer);
    }
    
    // Exons
    if (query_->get_exon_count() > 0)
    {
        vector<BedQuery*> exons;
        inf_maxgap_ = imaxgap;
        for (unsigned i = 0; i < query_->get_exon_count(); ++i)
        {
            ref_it1 = references->begin();
            ref_it2 = --(references->end());
            set_ref_iterators(query_->get_start() +
                              (*(query_->get_exon_starts()))[i],
                              query_->get_start() +
                              (*(query_->get_exon_ends()))[i],
                              *references, ref_it1, ref_it2);
            exons.push_back(get_mapping(query_->get_start() +
                                        (*(query_->get_exon_starts()))[i],
                                        query_->get_start() +
                                        (*(query_->get_exon_ends()))[i],
                                        *references, informants, ref_it1,
                                        ref_it2));
        }
        if (!(answer->merge_exons(exons)) && !alwaysmap_)
            error("no_exon_mapping");
    }
    return answer;
}