#ifndef MAPPING_H
#define MAPPING_H

#include <map>
#include <string>
#include <vector>
#include <stdexcept>

#include "../../../ocaml-bgzf/bgzf.h"

#include "Query.h"
#include "Sequence.h"
#include "IOHandler.h"


class MappingError: public std::runtime_error
{
    public:
        MappingError() throw(): std::runtime_error("MappingError") {};
        MappingError(const char m[]) throw(): std::runtime_error(m) {}; 
        MappingError(const MappingError &other) throw():
            std::runtime_error(other) { }
};

class Mapping
{
    public:
        Mapping(IOHandler* ioh, std::string &informant, int inf_maxgap,
                int ref_maxgap, bool inner, bool alwaysmap,
                std::map<std::string, bioid_t> *genome_map,
                std::vector< std::map<std::string,
                std::pair <bioid_t, seqpos_t> > > *chr_maps,
                std::map <bioid_t, std::vector<IndexItem*> > *index);
        ~Mapping();
        
        void add_error(std::string error);
        void set_query(BedQuery* qry);
        BedQuery* get_answer();
        
        void print_errors();
    
    private:
        IOHandler* ioh_;
        std::string informant_;
        bioid_t inf_id_;
        bioid_t ref_chr_id_;
        int inf_maxgap_, ref_maxgap_;
        bool inner_, alwaysmap_;
        BedQuery* query_;
        std::map<std::string, bioid_t> *genome_map_;
        std::vector< std::map<std::string,
                    std::pair <bioid_t, seqpos_t> > > *chr_maps_;
        std::map <bioid_t, std::vector<IndexItem*> > *index_;
        std::map <bioid_t, std::pair <std::string, seqpos_t> > *id_to_len_;
        std::vector < std::string > errors_;
        static const int known_error_count_ = 10;
        std::string known_error_names_[known_error_count_] = {"no_mapping",
            "pos_to_gap", "inf_preceed", "inf_strand", "inf_contig", "inf_gap",
            "invalid_query", "no_exon_mapping", "no_thick_mapping", "ref_gap"};
        std::string known_error_messages_[known_error_count_] = {
            "There is no mapping of the interval (maybe try --outer?)",
            "Position maps to gap",
            "In informant: one sequence does not preceed the next one",
            "In informant: sequences are from different strands",
            "In informant: sequences are from different contigs",
            "In informant: there is a gap of width ",
            "The query is invalid",
            "There is no mapping of the exons "
                "(could be overriden by -alwaysmap)",
            "There is no mapping of the thick region "
                "(could be overriden by -alwaysmap)",
            "In reference: there is a gap of width "};
        int found_gap_ = 0;
        
        std::vector<Reference*>* get_references(seqpos_t start, seqpos_t end);
        BedQuery* map_position(std::vector <Reference*> &references,
                               std::vector<Informant*> &informants,
                               seqpos_t position, int way,
                               std::vector<Reference*>::iterator &ref_it,
                               std::vector<Informant*>::iterator &inf_it_ret);
        seqpos_t min(seqpos_t x, seqpos_t y);
        seqpos_t max(seqpos_t x, seqpos_t y);
        bool check_informants(std::vector<Informant*>::iterator &inf_it1,
                              std::vector<Informant*>::iterator &inf_it2);
        BedQuery* get_mapping(seqpos_t start, seqpos_t end,
                              std::vector<Reference*> &references,
                              std::vector<Informant*> &informants,
                              std::vector<Reference*>::iterator &ref_it1,
                              std::vector<Reference*>::iterator &ref_it2);
        std::string get_error_message(std::string error_name);
        void error(std::string error_name);
        void set_ref_iterators(seqpos_t start, seqpos_t end,
                               std::vector<Reference*> &references,
                               std::vector<Reference*>::iterator &ref_it1,
                               std::vector<Reference*>::iterator &ref_it2);
        void fill_informant_vector(std::vector<Reference*> &references,
                                   std::vector<Informant*> &informants);
        seqpos_t get_inf_count(std::vector<Reference*>::iterator &from,
                               const std::vector<Reference*>::iterator &to);
};

#endif /* MAPPING_H */