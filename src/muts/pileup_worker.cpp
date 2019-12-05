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

#include "pileup_worker.hpp"
#include "collapse_aux.hpp"

using namespace gwsc;


void PileupWorker::operator()() {
    BamBuffer::rpair range;
    pcount = 0;
    coverage.clear();
    while(buffer_->get_next(range)){
        if(range.first == range.second) break;
        gids_.clear();
        for(auto it = range.first; it != range.second; it++){
            BamDetail & d = *(*it);
            d.corrected = false;
            d.xt = bam_aux2A(bam_aux_get(d.b, "RE"));
            d.dup = (d.b->core.flag & BAM_FDUP) > 0;
            tmp_ = bam_aux2Z(bam_aux_get(d.b, "CB"));
            if(cellranger_){
                std::size_t pos = tmp_.find("-"); 
                tmp_ = tmp_.substr(0, pos);
            }
            d.barcode = bhash_(tmp_);
            tmp_ = bam_aux2Z(bam_aux_get(d.b, "UB"));
            seq2int<gwsc::ADNA4, uint32_t>(tmp_, d.umi);
            d.make_hash();
            d.lft = d.b->core.pos;
            d.rgt = bam_endpos(d.b) - 1;

            if(bam_is_rev(d.b)){
                d.pos = d.rgt;
            }else{
                d.pos = d.lft;
            }
            if(gids_.empty() || d.gid != gids_.back()) gids_.push_back(d.gid);
            reads++;
        }
        //std::cout << "Pileup " << (range.second - range.first) << " reads \n";
        unsigned int base_pos_delta = build_genes_();


        pup_.build_reads(range);

        while(pup_.next(out_)){
            /*
            std::sort(out_.begin(), out_.end(), [](const PileupOut & r1, const PileupOut & r2) {
                //return std::make_tuple(r1.base, r1.rev, r1.d->barcode) < std::make_tuple(r2.base, r2.rev, r2.d->barcode);
                return std::make_tuple(r1.d->barcode, r1.base, r1.rev) < std::make_tuple(r2.d->barcode, r2.base, r2.rev);
            });
            */
            pcount++;
            bcounter_.clear();
            if(pcount >= positions.size()){
                positions.push_back(PositionCount());
            }

            auto & p = positions[pcount - 1];
            p.reset();
            p.tid = pup_.tid;
            p.pos = pup_.pos;
            p.ref = genome_[p.tid].seq[p.pos];
            switch(p.ref){
                case 'A': p.refi = 0; break;
                case 'C': p.refi = 1; break;
                case 'G': p.refi = 2; break;
                case 'T': p.refi = 3; break;
                default:  p.refi = 4; break;
            }
            unsigned int plus_cov = 0, minus_cov = 0;
            bool has_plus = false, has_minus = false;
            for(auto & r : out_){
                if(r.base < 4 && r.qual >= min_qual_){
                    has_minus |= r.rev;
                    has_plus |= !r.rev;
                    auto & b = p.bases[r.base];
                    if(r.rev){
                        b.m_count++;
                        //b.m_gaps += (r.dbases + r.ibases);
                        //b.m_NM += r.NM;
                        //b.m_lens += r.abases;
                        minus_cov++;
                    }else{
                        b.p_count++;
                        //b.p_gaps += (r.dbases + r.ibases);
                        //b.p_NM += r.NM;
                        //b.p_lens += r.abases;
                        plus_cov++;
                    }
                    p.coverage++;
                    auto it = bcounter_.insert(std::make_pair(r.d->barcode, std::array<uint16_t, 8>{}));
                    it.first->second[r.base + 4 * r.rev]++;
                }else{
                    p.ambig++;
                }
            }
            count_barcodes_(p);
            if(p.barcodes < min_barcodes_ || p.ref == 'N'){
                pcount--;
            }else{
                //std::cout << p.tid << " " << p.pos << " coverage = " << p.coverage << " has = " << has_plus << "/" << has_minus << " pcount = " << pcount << "\n";
                bool passed = false;

                unsigned int rc = 0, ac = 0;
                for(size_t i = 0; i < 4; i++){
                    unsigned int c = p.bases[i].p_count + p.bases[i].m_count;
                    unsigned int b = p.bases[i].t_barcodes;

                    char bs = "ACGT"[i];
                    if(bs == p.ref){
                        rc = c;
                    }else if(b >= min_alternative_ && c > ac){
                        ac = c;
                    }
                }
                if((1.0 * ac / (ac + rc)) >= min_af_){
                    passed = true;
                }

                uint16_t bt = 0;
                char pbase = '-', mbase = '-';
                unsigned int pidx = 6, midx = 6;
                if(has_plus){
                    pidx = 5;
                    pbase = 'I';
                }
                if(has_minus){
                    midx = 5;
                    mbase = 'I';
                }

                if(p.pos >= base_pos_delta && p.pos < (base_pos_delta + btypes_.size())){
                    bt = btypes_[p.pos - base_pos_delta];
                    if(has_minus){
                        if((bt & SBaseType::MEXON)){
                            midx = 0;
                            mbase = 'C';
                        }else if((bt & SBaseType::MNC)){
                            midx = 1;
                            mbase = 'E';
                        }else if(bt & SBaseType::M3UTR){
                            midx = 2;
                            mbase = '3';
                        }else if(bt & SBaseType::M5UTR){
                            midx = 3;
                            mbase = '5';
                        }else if(bt & SBaseType::MINTRON){
                            midx = 4;
                            mbase = 'N';
                        }
                    }

                    if(has_plus){
                        if((bt & SBaseType::PEXON)){
                            pidx = 0;
                            pbase = 'C';
                        }else if((bt & SBaseType::PNC)){
                            pidx = 1;
                            pbase = 'E';
                        }else if(bt & SBaseType::P3UTR){
                            pidx = 2;
                            pbase = '3';
                        }else if(bt & SBaseType::P5UTR){
                            pidx = 3;
                            pbase = '5';
                        }else if(bt & SBaseType::PINTRON){
                            pidx = 4;
                            pbase = 'N';
                        }
                    }
                }
                p.pbase = pbase;
                /*
                std::cout << p.tid << " " << p.pos << " coverage = " << p.coverage << " bases = " << pbase << "/" << mbase << " has = " << has_plus << "/" << has_minus << " ref = " << p.ref << "\n";
                std::cout << "  barcodes = " << p.barcodes << " bcounts = " << p.bases[0].t_barcodes;
                for(size_t j = 1; j < 4; j++){
                    std::cout << ", " << p.bases[j].t_barcodes;
                }
                std::cout << "\n";
                */
                p.mbase = mbase;
                bases++;
                coverage.push_back(PositionCoverage(p.tid, p.pos, plus_cov, minus_cov));
                if(has_plus)  plus_bases++;
                if(has_minus) minus_bases++;
                if(midx < 6) total_rates.mrates[midx + 6 * (3 - p.refi)]++;
                if(pidx < 6) total_rates.prates[pidx + 6 *  p.refi]++;

                /*
                for(auto & b : p.bcounts){
                    std::cout << "  barcode = " << b.barcode << " plus = " << b.pbases[0];
                    for(size_t j = 1; j < 4; j++){
                        std::cout << ", " << b.pbases[j] ;
                    }
                    std::cout << " minus = " << b.mbases[0];

                    for(size_t j = 1; j < 4; j++){
                        std::cout << ", " << b.mbases[j] ;
                    }
                    std::cout << "\n";
                }
                */

                for(auto & b : bcounter_){
                    if(barcode_rates.size() < b.first){
                        barcode_rates.resize(b.first + 1);
                    }
                    if(midx < 6 && std::accumulate(b.second.begin() + 4, b.second.end(),0)){
                        barcode_rates[b.first].mrates[midx + 6 * (3 - p.refi)]++;
                    }
                    if(pidx < 6 && std::accumulate(b.second.begin(),b.second.begin() + 4,0)){
                        barcode_rates[b.first].prates[pidx + 6 * p.refi]++;
                    }
                }

                if(!passed){
                    // Need to remove the barcode counts collected from this position
                    pcount--;
                }else{
                    pbases++;
                }
            }
        }
    }
}

void PileupWorker::count_barcodes_(PositionCount & p){
    p.barcodes = bcounter_.size();
    for(auto & b : bcounter_){
        p.bcounts.push_back(BarcodeCount(b.first));
        for(size_t i = 0; i < 4; i++){
            p.bases[i].t_barcodes += (b.second[i] + b.second[i + 4]) > 0;
            p.bcounts.back().pbases[i] += b.second[i];
            p.bcounts.back().mbases[i] += b.second[i + 4];
        }
    }
}

unsigned int PileupWorker::build_genes_(){

    std::sort(gids_.begin(), gids_.end());
    gids_.erase(std::unique(gids_.begin(), gids_.end()), gids_.end());
    btypes_.clear();

    unsigned int gl = std::numeric_limits<unsigned int>::max(), gr = 0;
    for(auto gid : gids_){
        const auto & g = txidx_.gene(gid);
        gl = std::min(gl, g.lft);
        gr = std::max(gr, g.rgt);
    }
    btypes_.resize(gr - gl + 1);
    for(auto gid : gids_){
        const auto & g = txidx_.gene(gid);
        for(auto & intron : g.introns){
            for(size_t i = intron.lft; i <= intron.rgt; i++){
                btypes_[i - gl] |= g.strand == '+' ? SBaseType::PINTRON : SBaseType::MINTRON;
            }
        }

        for(size_t txid = g.tstart; txid < g.tend; txid++){
            const auto & t = txidx_.transcript(txid);
            for(size_t ei = 0; ei < t.exons.size(); ei++){
                auto const & e = t.exons[ei];
                for(size_t i = e.lft; i <= e.rgt; i++){
                    if(t.coding_start > -1 && static_cast<int>(i) < t.coding_start){
                        btypes_[i - gl] |= t.strand == '+' ? SBaseType::P5UTR : SBaseType::M3UTR;
                    }else if(t.coding_end > -1 && static_cast<int>(i) > t.coding_end){
                        btypes_[i - gl] |= t.strand == '+' ? SBaseType::P3UTR : SBaseType::M5UTR;
                    }else if(t.coding_start > -1){
                        btypes_[i - gl] |= t.strand == '+' ? SBaseType::PEXON : SBaseType::MEXON;
                    }else{
                        btypes_[i - gl] |= t.strand == '+' ? SBaseType::PNC : SBaseType::MNC;
                    }
                }
            }
        }

    }
    return gl;
}
