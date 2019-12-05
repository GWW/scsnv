/*
Copyright (c) 2018-2019 Gavin W. Wilson

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

#include "build.hpp"
#include <bwa/bwa.h>
#include "../util/task_log.hpp"
#include "../util/gzstream.hpp"
#include "aux.hpp"
#include <fstream>
#include <locale>

using namespace gwsc;

size_t write_transcript(gzFile zout, gzofstream & zout_tx, size_t tid, size_t gid, size_t tx_id, 
        const Sequence & seq, const Transcript & t, std::vector<Block> & exons, unsigned int chunk=60)
{
    Sequence s;
    for(size_t i = 0; i < t.children.size(); i++){
        auto const & e = t.children[i];
        exons.push_back(e);
        s.append(seq.begin() + e.lft, seq.begin() + e.rgt + 1);
    }

    if(t.strand == '-'){
        s.reverse_cmpl();
    }
    //cout << "  " << t.id << " length = " << s.size() << "\n";
    std::string n = ">" + t.id + "\n";
    gzwrite(zout, &n[0], n.length());
    for(unsigned int i = 0; i < s.size(); i += chunk){
        unsigned int w = std::min(i + chunk, static_cast<unsigned int>(s.size()));
        if(w == 0) break;
        n.clear();
        n.assign(s.begin() + i, s.begin() + w);
        n += '\n';
        gzwrite(zout, &n[0], n.length());
    }

    zout_tx << t.ref << "\t" << tid << "\t" << gid << "\t" << tx_id << "\t" << t.id << "\t" << t.name << "\t"
        << t.lft << "\t" << t.rgt << "\t" << t.strand << "\t" << t.children.size();

    zout_tx << "\t" << t.children.front().lft;
    for(size_t i = 1; i < t.children.size(); i++) zout_tx << "," << t.children[i].lft;
    zout_tx << "\t" << t.children.front().rgt;
    for(size_t i = 1; i < t.children.size(); i++) zout_tx << "," << t.children[i].rgt;

    if(t.biotype == "protein_coding"){
        zout_tx << "\t" << (t.cds_start == -1 ? t.lft : t.cds_start)
                << "\t" << (t.cds_end == -1 ? t.rgt : t.cds_end)
                << "\n";
    }else{
        zout_tx << "\t-1\t-1\n";

    }
    return s.size();
}

void write_gene_introns(gzofstream & zout, unsigned int tid, unsigned int gid, const Gene & g, std::vector<Block> & exons){
    std::sort(exons.begin(), exons.end());
    if(exons.size() > 1){
        auto res = exons.begin();
        auto it = res;
        while (++it != exons.end()) {
            //IF they don't overlap add after res
            if(!res->overlaps(*it) && ++res != it){
                *res = std::move(*it);
            }else{
                res->lft = std::min(res->lft, it->lft);
                res->rgt = std::max(res->rgt, it->rgt);
            }
        }
        exons.erase(++res, exons.end());
    }

    zout << g.ref << "\t" << tid << "\t" << gid << "\t" << g.id << "\t" << g.name << "\t" << g.lft << "\t" << g.rgt << "\t" << g.strand;
    if(exons.size() < 2){
        zout << "\t0\t\t\n";
    }else{
        std::vector<Block> introns;
        for(size_t i = 1; i < exons.size(); i++){
            introns.push_back(Block(exons[i - 1].rgt + 1, exons[i].lft - 1));
        }
        zout << "\t" << introns.size();
        zout << "\t" << introns.front().lft;
        for(size_t i = 1; i < introns.size(); i++) zout << "," << introns[i].lft;
        zout << "\t" << introns.front().rgt;
        for(size_t i = 1; i < introns.size(); i++) zout << "," << introns[i].rgt;
        zout << "\n";
    }
}


TXIndexBuild::TXIndexBuild(const std::string & prefix, bool build_bwa) 
    : prefix_(prefix), build_bwa_(build_bwa)
{
}

void TXIndexBuild::build(const std::string & ref, const std::string & gtf, unsigned int tx_length, bool retained_introns){
    //TODO: Also write a list of chromosome sizes
    std::cout.imbue(std::locale(""));
    tout << "Loading the GTF file" << std::endl;
    GeneModel    gm;
    parse_GTF(gtf, gm);
    tout << "Building Transcript Index" << std::endl;
    gm.make_map();
    FastaReader far(ref);
    size_t tshort = 0, wrbases = 0, wbases = 0, retained = 0, txid = 0, gid = 0, tid = 0;
    Fasta fa;
    std::string out = prefix_ + "_transcripts.fa.gz";
    gzFile zout = gzopen(out.c_str(), "wb");

    gzofstream zout_tx(prefix_ + "_transcripts.txt.gz");
    zout_tx << "ref\ttid\tgid\ttx_id\ttranscript_id\ttranscript_name\tlft\trgt\tstrand\texons\telfts\tergts\tcoding_start\tcoding_end\n";

    gzofstream zout_genes(prefix_ + "_genes.txt.gz");
    zout_genes << "ref\ttid\tgid\tgene_id\tgene_name\tlft\trgt\tstrand\tintrons\tilfts\tirgts\n";

    std::ofstream lout(prefix_ + "_lenghts.txt");
    std::vector<Block> exons;
    while(far.read(fa)){
        size_t total_transcripts = 0, bases = 0, total_retained = 0, total_short = 0;
        size_t total_genes = 0, rbases = 0;
        auto it = gm.chroms.find(fa.name);
        lout << fa.name << "\t" << fa.seq.size() << "\t" << fa.comment << "\n";
        if(it == gm.chroms.end()) {
            tid++;
            continue;
        }
        for(auto & g : it->second){
            bool used = false;
            exons.clear();
            for(const auto & t : g.children){
                if(t.tlen() < tx_length){
                    total_short++;
                }else if(t.biotype != "retained_intron" && !t.children.empty()){
                    bases += write_transcript(zout, zout_tx, tid, gid, txid, fa.seq, t, exons);
                    txid++;
                    total_transcripts++;
                    used = true;
                }else{
                    if(retained_introns){
                        rbases += write_transcript(zout, zout_tx, tid, gid, txid, fa.seq, t, exons);
                        txid++;
                        used = true;
                    }
                    total_retained++;
                }
            }
            if(used){
                write_gene_introns(zout_genes, tid, gid, g, exons);
                gid++;
                total_genes++;
            }
        }
        tout << "Processed " << fa.name << " bases = " << bases << ", genes = " << total_genes << " transcripts = " << total_transcripts;
        if(retained_introns){
            std::cout << ", retained intron bases = " << rbases << " retained intron transcripts = " << total_retained;
        }else{
            std::cout << ", skipped retained intron transcripts = " << total_retained;
        }
        std::cout << std::endl;

        wbases += bases;
        wrbases += rbases;
        retained += total_retained;
        tshort += total_short;
        tid++;
    }

    gzclose(zout);
    zout_tx.close();
    zout_genes.close();

    tout << "Done bases = " << wbases << ", genes = " << gid << " transcripts = " << txid;
    if(retained_introns){
        std::cout << ", retained intron bases = " << wrbases << " retained intron transcripts = " << retained;
    }else{
        std::cout << ", skipped retained intron transcripts = " << retained;
    }
    std::cout << std::endl;

    out = prefix_ + "_transcripts.fa.gz";
    std::string iout = prefix_ + "_bwa";
    if(build_bwa_){
        tout << "Building the BWA-mem index\n";
        bwa_idx_build(out.c_str(), iout.c_str(), BWTALGO_AUTO, 1000000000);
    }else{
        tout << "Skipping the BWA-mem index step\n";
    }
    tout << "Done!\n";
}
