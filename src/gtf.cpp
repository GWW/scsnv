/*
Copyright (c) 2018-2020 Gavin W. Wilson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include "gtf.hpp"
#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>
#include "misc.hpp"
#include "tokenizer.hpp"
#include "read_buffer.hpp"
#include <limits>
#include <array>
#include <cstring>

using namespace std;
using namespace gwsc;

void gwsc::parse_GTF(const std::string & file, GeneModel & m) {
    unordered_map<string, Gene * >       gmap;
    unordered_map<string, Transcript * > tmap;
    //split_tokens toks, ttoks;
    string line;
    FileWrapper in(file);

    /*
    if(!in){
        cerr << "Error opening " << file << " for parsing\n";
    }
    */

    unsigned int line_no = 0;
    std::array< pair< string, string >, 6 > tags { {
        { "gene_name", string() },
        { "gene_id", string() },
        { "gene_biotype", string() },
        { "transcript_name", string() },
        { "transcript_id", string() },
        { "transcript_biotype", string() }
    } };

    string gid, tid, gbio, tbio, gname, tname;
    std::string last_chrom;
    auto chrom_ptr = m.chroms.end();

    Gene       * gptr = nullptr;
    Transcript * tptr = nullptr;

    //std::vector<char *> toks, ttoks;

    //Tokenizer::tokens toks, ttoks;
    Tokenizer::tokens ttoks;
    ParserTokens toks;
    while(in.tokenize_line(toks) >= 0){
        if(!toks.empty() && toks[0][0] == '#') {
            line_no++;
            continue;
        }
	if(toks.size() < 9 || toks.size() > 10){
	    std::cout << "Error malformed GTF file: " << file << " at line: " << line_no << " number of toks = " << toks.size() << "\n";
	    exit(1);
	}

        bool exon = false, stop_codon = false, start_codon = false;
        if(toks[2] ==  "exon"){
            exon = true;
        }else if(toks[2] == "start_codon"){
            start_codon = true;
        }else if(toks[2] == "stop_codon"){
            stop_codon = true;
        }
        /*
        if(strcmp(toks[2], "exon") == 0){
            exon = true;
        }else if(strcmp(toks[2], "start_codon") == 0){
            start_codon = true;
        }else if(strcmp(toks[2], "stop_codon") == 0){
            stop_codon = true;
        }
        */

        if(!exon && !start_codon && !stop_codon){
            continue;
        }

        for(auto & t : tags){
            t.second.clear();
        }

        Tokenizer::get(toks[8], ';', ttoks);

        for(auto & tk : ttoks){
            while(*tk == ' ') tk++;
            for(auto & t : tags){
                if(strncmp(t.first.c_str(), tk, t.first.size()) == 0){
                    char *st = strchr(tk, '"');
                    while(*(++st) != '"'){
                        t.second.append(1, *st);
                    }
                    //cout << "Found " << t.first << " val = " << t.second << "\n";
                }
            }
        }

        if(tags[1].second.empty() || tags[4].second.empty()){
	    std::cout << "Error malformed GTF file: " << file << " at line: " << line_no << " number of toks = " << toks.size() << "\n";
	    exit(1);
        }

        if(toks[0] != last_chrom){
            chrom_ptr = m.chroms.insert(make_pair(toks[0], GeneModel::genes())).first;
            last_chrom = toks[0];
        }

        if(!gptr || gptr->id != tags[1].second){
            auto res = gmap.insert(make_pair(tags[1].second, nullptr));
            if(res.second){
                chrom_ptr->second.push_back(Gene());
                gptr = &chrom_ptr->second.back();
                gptr->id = tags[1].second;
                gptr->name = tags[0].second.empty() ? tags[1].second : tags[0].second;
                gptr->strand = toks[6][0];
                gptr->ref = toks[0];
                gptr->biotype = tags[2].second;
            }
        }

        if(!tptr || tptr->id != tags[4].second){
            auto res = tmap.insert(make_pair(tags[4].second, nullptr));
            if(res.second){
                gptr->children.push_back(Transcript());
                tptr = &gptr->children.back();
                tptr->id = tags[4].second;
                tptr->name = tags[3].second.empty() ? tags[4].second : tags[3].second;
                tptr->strand = toks[6][0];
                tptr->ref = toks[0];
                tptr->biotype = tags[5].second;
            }
        }

        if(exon){
            tptr->children.push_back(Exon(std::stoi(toks[3]) - 1, std::stoi(toks[4]) - 1, toks[6][0]));
        }else if(start_codon){
            if(toks[6][0] == '+'){
                tptr->cds_start = std::stoi(toks[3]) - 1;
            }else{
                tptr->cds_end = std::stoi(toks[4]) - 1;
            }
        }else if(stop_codon){
            if(toks[6][0] == '+'){
                tptr->cds_end = std::stoi(toks[4]) - 1;
            }else{
                tptr->cds_start = std::stoi(toks[3]) - 1;
            }
        }
    }
    sort_model(m);
}
