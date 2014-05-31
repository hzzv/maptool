#include <string>
#include <sstream>
#include <vector>

#include <iostream>

#include "include/Query.h"

using std::string;
using std::to_string;
using std::istringstream;
using std::vector;

bool IndexItem::get_strand()
{
    return strand_;
}

seqpos_t IndexItem::get_chr_pos()
{
    return chr_pos_;
}

seqpos_t IndexItem::get_bases_count()
{
    return bases_count_;
}

uint64_t IndexItem::get_pointer()
{
    return pointer_;
}

// Set numbers given as string separated by ',' in given vector
void BedQuery::set_numbers(string str, vector<seqpos_t>* numbers)
{
    seqpos_t num = -1;
    for (unsigned int i = 0; i < str.size(); ++i)
    {
        if (str[i] == ' ') continue;
        if (str[i] == ',')
        {
            numbers->push_back(num);
            num = -1;
        }
        if (str[i] != ',')
        {
            if (num == -1) num = 0;
            num *= 10;
            num += str[i]-'0';
        }
    }
    if (num != -1) numbers->push_back(num);
}

// Parse given BED line to a BedQuery
BedQuery::BedQuery(string &bedline)
{
    istringstream bedstream(bedline);
    string strand, exon_sizes_str, exon_starts_str;
    end_ = -1;
    name_ = "default_name";
    score_ = -1;
    thick_start_ = -1;
    thick_end_ = -1;
    rgb_ = "0,0,0";
    exon_count_ = 0;
    bedstream >> chromosome_ >> start_ >> end_ >> name_ >> score_ >> strand >>\
        thick_start_ >> thick_end_ >> rgb_ >> exon_count_ >> exon_sizes_str >>\
        exon_starts_str;
    if (strand.size() == 0) strand = "+";
    strand_ = (strand[0] == '+');
    original_strand_ = strand_;
    if (exon_sizes_str.size() == 0) exon_count_ = 0;
    vector<seqpos_t>* exon_sizes_vec = new vector<seqpos_t>;
    exon_starts_ = new vector<seqpos_t>;
    exon_ends_ = new vector<seqpos_t>;
    set_numbers(exon_sizes_str, exon_sizes_vec);
    set_numbers(exon_starts_str, exon_starts_);
    if (exon_sizes_vec->size() < exon_count_)
        exon_count_ = exon_sizes_vec->size();
    if (exon_starts_->size() < exon_count_)
        exon_count_ = exon_starts_->size();
    for (unsigned i = 0; i < exon_count_; ++i)
    {
        exon_ends_->push_back((*exon_starts_)[i] + (*exon_sizes_vec)[i]);
    }
    
    chr_size_ = -1;
    closed_ = false;
    
    delete exon_sizes_vec;
}

// Create an empty BedQuery
BedQuery::BedQuery()
{
    start_ = -1;
    end_ = -1;
    thick_start_ = -1;
    thick_end_ = -1;
    chr_size_ = -1;
    score_ = -1;
    strand_ = true;
    original_strand_ = strand_;
    rgb_ = "";
    name_ = "";
    chromosome_ = "";
    exon_count_ = 0;
    exon_starts_ = new vector<seqpos_t>;
    exon_ends_ = new vector<seqpos_t>;
    closed_ = false;
}

// Create BedQuery from given mapping of a position
BedQuery::BedQuery(const BedQuery &bq, string &chromosome, bool strand,
                   seqpos_t chr_size, seqpos_t start)
{
    chromosome_ = string(chromosome);
    start_ = start;
    end_ = -1;
    name_ = string(bq.name_);
    score_ = bq.score_;
    strand_ = strand;
    original_strand_ = strand_;
    thick_start_ = -1;
    thick_end_ = -1;
    rgb_ = bq.rgb_;
    exon_count_ = 0;
    exon_starts_ = new vector<seqpos_t>;
    exon_ends_ = new vector<seqpos_t>;
    
    chr_size_ = chr_size;
    closed_ = true;
}

BedQuery::~BedQuery()
{
    delete this->exon_starts_;
    delete this->exon_ends_;
}

string BedQuery::get_chr()
{
    return chromosome_;
}

seqpos_t BedQuery::get_start()
{
    return start_;
}

seqpos_t BedQuery::get_end()
{
    return end_;
}

string BedQuery::get_name()
{
    return name_;
}

bool BedQuery::get_strand()
{
    return strand_;
}

seqpos_t BedQuery::get_thick_start()
{
    return thick_start_;
}

seqpos_t BedQuery::get_thick_end()
{
    return thick_end_;
}

unsigned BedQuery::get_exon_count()
{
    return exon_count_;
}

vector<seqpos_t>* BedQuery::get_exon_starts()
{
    return exon_starts_;
}

vector<seqpos_t>* BedQuery::get_exon_ends()
{
    return exon_ends_;
}

// Create BED line
string BedQuery::get_bedline()
{
    string strand = "+";
    if (!strand_) strand = "-";
    string bedline = chromosome_ + "\t" + to_string(start_) + "\t" +
        to_string(end_);
    if (name_.compare("default_name") != 0 or score_ != -1)
        bedline += "\t" + name_;
    if (score_ != -1) bedline += "\t" + to_string(score_);
    if (!strand_ or (thick_start_ != -1 and thick_end_ != -1))
        bedline += "\t" + strand;
    if (thick_start_ != -1 and thick_end_ != -1)
        bedline += "\t" + to_string(thick_start_) + "\t" +
            to_string(thick_end_);
    if (rgb_.compare("0,0,0") != 0 or exon_count_ > 0)
        bedline += "\t" + rgb_;
    if (exon_count_ > 0)
    {
        bedline += "\t" + to_string(exon_count_) + "\t";
        int c = 0;
        if (closed_) c = 1;
        for (unsigned i = 0; i < exon_count_; ++i)
        {
            if (i != 0) bedline += ",";
            bedline += to_string((*exon_ends_)[i] - (*exon_starts_)[i] + c);
        }
        bedline += "\t";
        for (unsigned i = 0; i < exon_count_; ++i)
        {
            if (i != 0) bedline += ",";
            bedline += to_string((*exon_starts_)[i]);
        }
    }
    return bedline;
}

// Merge given 'query' with 'this'
bool BedQuery::merge_query(BedQuery *query, bool query_strand)
{
    if ((strand_ != query->strand_) ||
        (chromosome_.compare(query->chromosome_) != 0) ||
        (start_ > query->get_start()) ||
        (closed_ != query->closed_))
        return false;
    seqpos_t start = start_, end = query->get_start();
    if (!strand_)
    {
        start = chr_size_ - query->get_start() - 1;
        end = chr_size_ - start_ - 1;
    }
    else if (!query_strand and strand_) strand_ = false;
    if (!query_strand and !original_strand_) strand_ = true;
    start_ = start;
    end_ = end;
    return true;
}

// Add given 'query' as thick interval to 'this'
bool BedQuery::merge_thick(BedQuery *query)
{
    if ((strand_ != query->strand_) ||
        (chromosome_.compare(query->chromosome_) != 0) ||
        (closed_ != query->closed_)) return false;
    thick_start_ = query->start_;
    thick_end_ = query->end_;
    return true;
}

// Add given 'queries' as exons to 'this'
bool BedQuery::merge_exons(std::vector<BedQuery*> &queries)
{
    if (queries.size() == 0)
    {
        exon_count_ = 0;
        return true;
    }
    int way = 1, m = 0, n = queries.size();
    if (!original_strand_)
    {
        way = -1;
        m = queries.size()-1;
        n = -1;
    }
    for (int i = m + way; i != n; i += way)
    {
        if ((queries[i-way]->strand_ != queries[i]->strand_) ||
            (queries[i-way]->chromosome_.compare(queries[i]->chromosome_) != 0) ||
            (queries[i-way]->end_ >= queries[i]->start_) ||
            (queries[i-way]->closed_ != queries[i]->closed_))
            return false;
    }
    if ((strand_ != queries[0]->strand_) ||
        (chromosome_.compare(queries[0]->chromosome_) != 0) ||
        (closed_ != queries[0]->closed_))
        return false;
    for (int i = m; i != n; i += way)
    {
        exon_starts_->push_back(queries[i]->start_ - start_);
        exon_ends_->push_back(queries[i]->end_ - start_);
    }
    exon_count_ = queries.size();
    return true;
}

// Change representation from BED-like to MAF-like
void BedQuery::to_closed()
{
    if (closed_) return;
    --end_;
    --thick_end_;
    for (unsigned i = 0; i < exon_count_; ++i) --(*exon_ends_)[i];
    closed_ = true;
}

// Change representation from MAF-like to BED-like
void BedQuery::to_half_closed()
{
    if (!closed_) return;
    ++end_;
    ++thick_end_;
    for (unsigned i = 0; i < exon_count_; ++i) ++(*exon_ends_)[i];
    closed_ = false;
}