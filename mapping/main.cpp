
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include <cstring>
#include <cstdlib>

using std::map;
using std::vector;
using std::string;
using std::pair;
using std::getline;

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include "include/Sequence.h"
#include "include/Query.h"
#include "include/IOHandler.h"
#include "include/Mapping.h"

const int WRONG_ARGNUM = 1, USAGE_ALL = 0, USAGE_PREP = 1, USAGE_BED = 2,
    USAGE_INFO = 3, FILE_INACCESSIBLE = 2, WRONG_ARGS = 3;

bool check_file_existence(char filename[])
{
    if (FILE *file = fopen(filename, "r"))
    {
        fclose(file);
        return true;
    }
    else
    {
        return false;
    }   
}

bool print_error(int error, int usage=USAGE_ALL, void* message = NULL)
{
    if ((error == WRONG_ARGNUM) || (error == WRONG_ARGS))
    {
        if (error == WRONG_ARGNUM)
        {
            std::cerr << "Wrong number of arguments. Usage:" << endl;
        }
        else std::cerr << "Usage:" << endl;
        //TODO: Implement preprocessing.
//         if (usage == USAGE_PREP || usage == USAGE_ALL)
//         {
//             std::cerr << "./maptool preprocess <alignment.maf> <header.bin> "
//                 "<compressed.bgzf>" << endl;
//         }
        if (usage == USAGE_BED || usage == USAGE_ALL)
        {
            std::cerr << "./maptool bed <header.bin> <compressed.bgzf> "
                "<informant> [--maxgap N] [--outer] [--alwaysmap]" << endl;
        }
        if (usage == USAGE_INFO || usage == USAGE_ALL)
        {
            std::cerr << "./maptool info <header.bin>" << endl;
        }
    }
    else if (error == FILE_INACCESSIBLE)
    {
        std::cerr << "File " << message << " is not accessible." << std::endl;
    }
    return false;
}

bool parse_options(char* opt[], int optnum, char command[], char file1[],
    char file2[], char file3[], char informant[], int &maxgap, bool &inner,
    bool &alwaysmap, bool &compressed)
{
    strcpy(file3, "");
    if (optnum <= 1) return false;
    if (strcmp(opt[1], "bed") == 0)
    {
        if (optnum < 5 || optnum > 10)
            return print_error(WRONG_ARGNUM, USAGE_BED);
        if (!check_file_existence(opt[2]))
            return print_error(FILE_INACCESSIBLE, 0, opt[2]);
        if (!check_file_existence(opt[3]))
            return print_error(FILE_INACCESSIBLE, 0, opt[3]);
        bool ok[] = {false, false, false, false};
        for (int i = 5; i < optnum; ++i) {
            if ((strcmp(opt[i], "--maxgap") == 0) && (optnum > i+1))
            {
                ok[i-5] = true;
                ok[i-4] = true;
                maxgap = atoi(opt[i+1]);
            }
            if (strcmp(opt[i], "--outer") == 0)
            {
                ok[i-5] = true;
                inner = false;
            }
            if (strcmp(opt[i], "--uncompressed") == 0)
            {
                ok[i-5] = true;
                compressed = false;
            }
            if (strcmp(opt[i], "--alwaysmap") == 0)
            {
                ok[i-5] = true;
                alwaysmap = true;
            }
        }
        for (int i = 5; i < optnum; ++i)
        {
            if (!ok[i-5]) return false;
        }
        strcpy(command, "bed");
        strcpy(file1, opt[2]);
        strcpy(file2, opt[3]);
        strcpy(informant, opt[4]);
    }
    else if (strcmp(opt[1], "info") == 0)
    {
        if (optnum != 3) return print_error(WRONG_ARGNUM, USAGE_INFO);
        if (!check_file_existence(opt[2]))
            return print_error(FILE_INACCESSIBLE, 0, opt[2]);
        strcpy(command, "info");
        strcpy(file1, opt[2]);
    }
    return true;
}

void delete_index(map <bioid_t, vector <IndexItem*> > &index)
{
    for (auto mit = index.begin(); mit != index.end(); ++mit)
    {
        for (auto vit = mit->second.begin(); vit != mit->second.end(); ++vit)
        {
            delete (*vit);
        }
    }
}

int main(int argc, char* argv[]) {
    char command[10];
    char file1[1000] = "", file2[1000] = "", file3[1000] = "";
    char informantc[NAME_SIZE] = "";
    int maxgap = 10;
    bool inner = true, alwaysmap = false, compressed = true;
    
    if (!parse_options(argv, argc, command, file1, file2, file3, informantc,
        maxgap, inner, alwaysmap, compressed))
    {
        print_error(WRONG_ARGS);
        exit(0);
    }
    string informant(informantc);
    
    map <string, bioid_t> genome_map;
    vector <map <string, pair <bioid_t, seqpos_t> > > chr_maps;
    map <bioid_t, vector <IndexItem*> > index;
    
    IOHandler ioh(file1, file2, compressed, file3);
    ioh.read_header(genome_map, chr_maps, index);
    
    if (strcmp(command, "info") == 0)
    {
        cout << "Informants (total " << genome_map.size()-1 << "):";
        for (map<string, bioid_t>::iterator it = genome_map.begin();
            it != genome_map.end(); ++it)
        {
            if (it->second != 0) cout << " " << it->first;
        }
        cout << endl;
        cout << "Reference chromosomes (total " << chr_maps[0].size() << "):";
        for (map <string, pair <bioid_t, seqpos_t> > ::iterator it =\
            chr_maps[0].begin(); it != chr_maps[0].end(); ++it)
        {
            cout << " " << it->first;
        }
        cout << endl;
    }
    else if (strcmp(command, "bed") == 0)
    {
        try
        {
            string informant(informantc);
            ioh.open_to_map();
            Mapping to_map(&ioh, informant, maxgap, maxgap, inner, alwaysmap,
                        &genome_map, &chr_maps, &index);
            string bedline;
            while (true)
            {
                // For each BED-line on input:
                getline(cin, bedline);
                if (cin.eof()) break;
                if (bedline.compare("") == 0) continue;
                BedQuery* bedquery = new BedQuery(bedline);
                // Transform it to closed interval
                bedquery->to_closed();
                // Try to map the interval
                to_map.set_query(bedquery);
                BedQuery* bq;
                try
                {
                    bq = to_map.get_answer();
                    bq->to_half_closed();
                    cerr << bq->get_name() << "\tmapped" << endl;
                    // If the mapping was successful, print it
                    cout << bq->get_bedline() << endl;
                }
                catch (MappingError e)
                {
                    to_map.print_errors();
                }
                to_map.delete_old();
                delete bedquery;
            }
        }
        catch (std::runtime_error e)
        {
            delete_index(index);
            throw e;
        }
        delete_index(index);
    }
}