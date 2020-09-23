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
#include "map_worker.hpp"
#include "sequence.hpp"
#include <set>
using namespace gwsc;

void MapWorker::operator()() {
    unsigned int N = in_(rps_, reads_, counts);
    while(N > 0){
        for(size_t i = 0; i < N; i++){
            map_(reads_[i]);
            if(data.countable){
                aligns.push_back(data.summary);
            }

            if(write_bam_){
                bout->write(data, reads_[i], buff, rcount, *tx_idx);
            }
            /*
            if(data.res == AlignGroup::CDNA){
                std::sort(data.alns.begin(), data.alns.begin() + data.acount);
                for(size_t i = 0; i < data.acount; i++){
                    auto & a = data.alns[i];
                    std::cout << "  " << i << " " << a.print() << "\n";
                }
                std::cout << "\n";
            }
            */
            //std::cout << "  " << reads_[i].name << " res = " << AlignGroup::alignres2str(data.res) << " tag = " << reads_[i].tag << "\n";
        }
        N = in_(rps_, reads_, counts);
        std::fill(counts.begin(), counts.end(), 0);
    }
}

/*
int find_poly(const std::string & tag, char target){
    char lb = ' ';
    int pcount = 0, pstart = 0;
    int count = 0, start = 0;
    int mcount = 0, mstart = 0;
    for(size_t i = 0; i < tag.size(); i++){
        char c = tag[i];
        if(c == target){
            if(lb == target){
                count++;
            }else{
                start = i;
                count = 1;
            }
        }else{
            if(count > 0){
                if(pcount > 0 && (start - (pstart + pcount)) == 1 && (count + pcount + 1) > mcount){
                    mstart = pstart;
                    mcount = count + pcount + 1;
                }
                if(count > mcount){
                    mstart = pstart;
                    mcount = count;
                }
                pstart = start;
                pcount = count;
            }
            start = -1;
            count = 0;
        }
        lb = c;
    }
    if(count > 0){
        if(pcount > 0 && (start - (pstart + pcount)) == 1 && (count + pcount + 1) > mcount){
            mstart = pstart;
            mcount = count + pcount + 1;
        }
        if(count > mcount){
            mstart = pstart;
            mcount = count;
        }
    }

    if(mcount >= 25){
        return mstart;
    }
    return tag.size() - 1;
}
*/

void MapWorker::map_(Read & read){
    data.reset();
    AlignSummary::bint barcode_index = 0;
    uint32_t umi_encoded = 0;
    int res = bc_(read.barcode, barcode_index);
    //std::cout << read.name << "\n";
    if(res == 1){
        counts[AlignGroup::BARCODE_FAIL]++;
        data.res = AlignGroup::BARCODE_FAIL;
        return;
    }

    //std::cout << "  Barcode Passed = " << read.barcode << " res = " << res << " Barcode index = " << barcode_index << "\n";

    if(res == 2) {
        barcode_corrected++;
        bc_rates[barcode_index][AlignGroup::BARCODE_FAIL]++;
    }else{
        barcode_correct++;
    }

    //Check the UMI
    if(!seq2int<gwsc::ADNA4, uint32_t>(read.umi, umi_encoded) || 
            read.umi.find_first_not_of(read.umi[0]) == std::string::npos){
        bc_rates[barcode_index][AlignGroup::UMI_FAIL]++;
        counts[AlignGroup::UMI_FAIL]++;
        data.res = AlignGroup::UMI_FAIL;
        return;
    }

    //std::cout << "  UMI Passed\n";

    int end = read.tag.size() - 1;
    if(smode_ == StrandMode::TAG_REV){
        /*
        if(internal_){
            end = find_poly(read.tag, 'T');
            if(end < ((int)read.tag.size() - 1)) itrimmed++;
        }
        */
        while(end >= 0 && read.tag[end] == 'T'){
            end--;
        }
    }else{
        /*
        if(internal_){
            end = find_poly(read.tag, 'A');
            if(end < ((int)read.tag.size() - 1)) itrimmed++;
        }
         */
        while(end >= 0 && read.tag[end] == 'A'){
            end--;
        }
    }
    unsigned int N = (read.tag.size() - end - 1);
    if(N < 6){
        end = read.tag.size() - 1;
    }else if(N > (read.tag.size() / 2) || (max_dust_ > 0 && dust_.calculate(read.tag, end + 1) < max_dust_)){
        data.res = AlignGroup::TAG_FAIL;
        bc_rates[barcode_index][AlignGroup::TAG_FAIL]++;
        counts[AlignGroup::TAG_FAIL]++;
        return;
    }

    //std::cout << read.tag << "\n";

    
    read.tend = end + 1;
    //std::cout << "Tend = " << read.tend << " out of " << read.tag.size() << "\n";


    data.summary.barcode = barcode_index;
    data.summary.umi = umi_encoded;
    data.summary.gene_id = std::numeric_limits<uint32_t>::max();
    tx_align->align(data, read.tag, read.tend);
    genome_align->align(data, read.tag, read.tend);
    //idx_.align(read.tag, data);
    // times 3 in case we do some end trimming
    int max_score = std::max(data.transcript_score, data.genome_score);
    //std::cout << "max score = " << max_score << "\n";
    if(max_score < static_cast<int>(tx_align->ascore.min_score)){
        data.res = AlignGroup::UNMAPPED;
        counts[data.res]++;
        bc_rates[barcode_index][data.res]++;
        //std::cout << "-----------------------\n\n";
        return;
    }

    //std::cout << "  Max score passed\n";

    max_score =  std::max(tx_align->ascore.min_score, max_score - (tx_align->ascore.max_diff * 3));
    int nmax_score = 0;
    //std::cout << "Max Score = " << max_score << " tend = " << read.tend << "\n";
    for(size_t i = 0; i < data.transcript_ar.n; i++){
        if(data.transcript_ar.a[i].score < max_score || data.transcript_ar.a[i].is_alt) continue;
        tx_align->get(data, i, read.tag, read.tend);
        //genome_align->verify(data.alns[data.acount - 1], read.tag);
        nmax_score = std::max(nmax_score, data.alns[data.acount - 1].score);
    }
    for(size_t i = 0; i < data.genome_ar.n; i++){
        if(data.genome_ar.a[i].score < max_score || data.genome_ar.a[i].is_alt) continue;
        genome_align->get(data, i, read.tag, read.tend);
        //genome_align->verify(data.alns[data.acount - 1], read.tag);
        nmax_score = std::max(nmax_score, data.alns[data.acount - 1].score);
    }
    std::sort(data.alns.begin(), data.alns.begin() + data.acount);
    max_score = nmax_score - tx_align->ascore.max_diff;
    std::array<unsigned int, AlignType::ELEM_COUNT> cnts{};
    uint32_t lgid = std::numeric_limits<uint32_t>::max();
    uint32_t igid = std::numeric_limits<uint32_t>::max();

    for(size_t i = 0; i < data.acount; i++){
        auto & a = data.alns[i];
        if(a.score >= max_score){
            a.passed = true;
            if(a.atype == AlignType::CDNA){
                if(a.gid != lgid){
                    cnts[AlignType::CDNA]++;
                    lgid = a.gid;
                }
            }else{
                bool count = true;
                if(a.atype == AlignType::INTRON){
                    igid = a.gid;
                }else if(a.atype == AlignType::AMBIGUOUS){
                    // Don't need to count this one
                    for(auto g : a.gids){
                        if(g == lgid){
                            count = false;
                            a.passed = false;
                            break;
                        }
                    }
                }
                cnts[a.atype] += count;
            }
        }
    }


    unsigned int gcount = cnts[AlignType::CDNA] + cnts[AlignType::INTRON] + cnts[AlignType::INTERGENIC];
    data.total_passed = 1;
    if(cnts[AlignType::INTRON] == 1 && cnts[AlignType::CDNA] == 1 && lgid == igid){
        //Ambiguous alignment where a read maps to an intron and cdna equally well
        data.res = AlignGroup::AMBIGUOUS;
    }else if(gcount > 1 || (cnts[AlignType::INTERGENIC] > 0 && cnts[AlignType::ANTISENSE]) > 0){
        //map to more than one gene or intergenic region
        //or map to an antisense and intergenic region
        data.res = AlignGroup::MULTIMAPPED;
        data.total_passed = gcount + cnts[AlignType::INTERGENIC] + cnts[AlignType::ANTISENSE];
    }else if(gcount == 0 && cnts[AlignType::AMBIGUOUS] > 0){
        data.res = AlignGroup::AMBIGUOUS;
    }else if(gcount == 0 && (cnts[AlignType::ANTISENSE] > 0)){
        data.res = AlignGroup::ANTISENSE;
    }else if(gcount == 1){
        if(cnts[AlignType::INTRON]){
            data.summary.intronic = true;
            data.summary.gene_id = igid;
            data.countable = true;
            data.res = AlignGroup::INTRONIC;
        }else if(cnts[AlignType::CDNA]){
            data.countable = true;
            data.summary.gene_id = lgid;
            data.res = AlignGroup::CDNA;
            /*
                std::cout << read.tag << "\n";
                for(size_t i = 0; i < data.acount; i++){
                    auto & a = data.alns[i];
                    if(a.score >= max_score){
                        std::cout << "  " << i << " " << a.print() << "\n";
                    }
                }
                std::cout << "\n";
            */
        }else{
            data.res = AlignGroup::INTERGENIC;
        }

        data.summary.pos = (static_cast<uint64_t>(data.alns[0].lft) << 32) | data.alns[0].rgt;

        /*
        if(data.alns[0].rev){
            data.summary.pos = data.alns[0].rgt;
        }else{
            data.summary.pos = data.alns[0].lft;
        }
        */
    }else{
        data.res = AlignGroup::AMBIGUOUS;
        /*
        for(size_t i = 0; i < data.acount; i++){
            auto & a = data.alns[i];
            if(a.score >= max_score){
                std::cout << "  " << i << " " << a.print() << "\n";
            }
        }

        for(size_t i = 0; i < AlignType::ELEM_COUNT; i++){
            std::cout << "  " << aligntype2str(static_cast<AlignType>(i)) << " " << cnts[i] << "\n";

        }
        std::cout << "-----------------------\n\n";
        */
    }

    counts[data.res]++;
    bc_rates[barcode_index][data.res]++;

    return;
}
