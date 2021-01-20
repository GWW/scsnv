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

#include "iprocessor.hpp"
#include "h5misc.hpp"

using namespace gwsc;


unsigned int base2int(char base){
    switch(base){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    };
    return 4;
}

void IndexProcessor::operator()(){
    auto rng = getrange_(5);
    while(rng.first < rng.second){
        for(size_t i = rng.first; i < rng.second; i++){
            AccSNVGroup & grp = (*snvs_)[i];
            AccSNV & snv = *grp.start;
            tid_ = bam_name2id(bh_, snv.chrom.c_str());
            curr_  = 0;
            read_ = 0;
            hts_itr_t *iter = sam_itr_queryi(bi_, tid_, snv.pos, grp.end_pos + 1);
            if(recs_.empty()) for(size_t i = 0; i < 100; i++) recs_.push_back(new BamDetail());
            int ret;
            while((ret=sam_itr_next(bf_, iter, recs_[curr_]->b)) > 0){
                read_++;
                if(!(*processor_)(*recs_[curr_], 0)) continue;
                if(filter_barcodes_){
                    btmp_ = bam_aux2Z(bam_aux_get(recs_[curr_]->b, "CB"));
                    size_t dash = btmp_.find('-');
                    if(dash != std::string::npos){
                        btmp_ = btmp_.erase(dash);
                    }
                    auto it = bchash_.find(btmp_);
                    if(it == bchash_.end()){
                        continue;
                    }
                    recs_[curr_]->barcode = it->second;
                }else{
                    recs_[curr_]->barcode = 0;
                    char xs = "+-"[bam_is_rev(recs_[curr_]->b)];
                    bam_aux_append(recs_[curr_]->b, "XS", 'A', 1, reinterpret_cast<const uint8_t*>(&xs));
                    recs_[curr_]->gid = std::numeric_limits<uint32_t>::max();
                }

                curr_++;
                if(curr_ >= recs_.size()) for(size_t i = 0; i < 100; i++) recs_.push_back(new BamDetail());
            }
            hts_itr_destroy(iter);

            if(curr_ > 0) {
                std::pair<std::vector<BamDetail*>::iterator, std::vector<BamDetail*>::iterator> rpair = {recs_.begin(), recs_.begin() + curr_};
                pileup_.reset(tid_, snv.pos);
                pileup_.process_range(rpair);
                auto xstart = grp.start;
                auto xend = grp.end;
                for(auto it = pileup_.positions.begin(); it != pileup_.positions.end(); it++){
                    while(xstart != xend && xstart->pos < it->pos){
                        xstart++;
                    }
                    if(xstart == xend) break;
                    if(xstart->pos == it->pos){
                        AccSNV & tsnv = *xstart;
                        uint32_t refi = base2int(tsnv.ref);
                        uint32_t alti = base2int(tsnv.alt);
                        auto & e = *it;
                        tsnv.pref = e.bases[refi].p_count;
                        tsnv.palt = e.bases[alti].p_count;
                        tsnv.mref = e.bases[refi].m_count;
                        tsnv.malt = e.bases[alti].m_count;

                        if(filter_barcodes_){
                            for(auto & b : e.bcounts){
                                calls.push_back(BarcodeCall(b.barcode, xstart->index, b.pbases[refi], b.mbases[refi], b.pbases[alti], b.mbases[refi]));
                            }
                            tsnv.r_barcodes = e.bases[refi].t_barcodes;
                            tsnv.a_barcodes = e.bases[alti].t_barcodes;

                            tsnv.pr_barcodes = e.bases[refi].p_barcodes;
                            tsnv.pa_barcodes = e.bases[alti].p_barcodes;
                            tsnv.par_barcodes = e.bases[alti].par_barcodes;

                            tsnv.mr_barcodes = e.bases[refi].m_barcodes;
                            tsnv.ma_barcodes = e.bases[alti].m_barcodes;
                            tsnv.mar_barcodes = e.bases[alti].mar_barcodes;

                            tsnv.tar_barcodes = e.bases[alti].tar_barcodes;
                            tsnv.total_barcodes = e.barcodes;
                        }
                    }
                }
            }
        }

        rng = getrange_(5);
    }
}

void DataManager::write(const std::string & out, const std::string & header, const std::vector<std::string> & barcode_ids) {
    std::string sout = out + ".txt.gz";
    gzofstream zout(sout);
    zout << header << "\tsnv_idx_h5";
    bool fbarcodes = processors.front()->fbarcodes();
    if(fbarcodes){
        zout 
            << "\t" << name << "_plus_ref"
            << "\t" << name << "_plus_alt"
            << "\t" << name << "_minus_ref"
            << "\t" << name << "_minus_alt"
            << "\t" << name << "_ref_barcodes"
            << "\t" << name << "_alt_barcodes"
            << "\t" << name << "_plus_ref_barcodes"
            << "\t" << name << "_plus_alt_barcodes"
            << "\t" << name << "_plus_altref_barcodes"
            << "\t" << name << "_minus_ref_barcodes"
            << "\t" << name << "_minus_alt_barcodes"
            << "\t" << name << "_minus_altref_barcodes"
            << "\t" << name << "_total_altref_barcodes"
            << "\t" << name << "_total_barcodes"
            << "\n";

    }else{
        zout 
            << "\t" << name << "_ref"
            << "\t" << name << "_alt"
            << "\n";
    }

    for(auto & s : snvs){
        zout << s.line << "\t" << s.index;
        if(fbarcodes){
            zout
                << "\t" << s.pref
                << "\t" << s.palt
                << "\t" << s.mref
                << "\t" << s.malt
                << "\t" << s.r_barcodes
                << "\t" << s.a_barcodes

                << "\t" << s.pr_barcodes
                << "\t" << s.pa_barcodes
                << "\t" << s.par_barcodes

                << "\t" << s.mr_barcodes
                << "\t" << s.ma_barcodes
                << "\t" << s.mar_barcodes

                << "\t" << s.tar_barcodes
                << "\t" << s.total_barcodes 
                << "\n";
        }else{
            zout
                << "\t" << (s.pref + s.mref)
                << "\t" << (s.palt + s.malt)
                << "\n";
        }
    }
    tout << name << " Done!!\n";

    if(fbarcodes){

        using namespace H5;
        std::string hout = out + "_barcodes.h5";
        H5File file(hout, H5F_ACC_TRUNC);
        //H5::Group group(file.createGroup("/barcode_rates"));
        std::vector<const char *> ctmp;
        for(auto & b : barcode_ids) ctmp.push_back(b.c_str());
        write_h5_string("barcodes", ctmp, file);

        std::vector<uint32_t> itmp;
        for(auto & c : calls) itmp.push_back(c.sindex);
        write_h5_numeric("snv_id_h5", itmp, file, PredType::NATIVE_UINT32);

        itmp.clear();
        for(auto & c : calls) itmp.push_back(c.barcode);
        write_h5_numeric("barcode_ids", itmp, file, PredType::NATIVE_UINT32);

        itmp.clear();
        for(auto & c : calls) itmp.push_back(c.pref);
        write_h5_numeric("plus_ref", itmp, file, PredType::NATIVE_UINT32);

        itmp.clear();
        for(auto & c : calls) itmp.push_back(c.mref);
        write_h5_numeric("minus_ref", itmp, file, PredType::NATIVE_UINT32);

        itmp.clear();
        for(auto & c : calls) itmp.push_back(c.palt);
        write_h5_numeric("plus_alt", itmp, file, PredType::NATIVE_UINT32);

        itmp.clear();
        for(auto & c : calls) itmp.push_back(c.malt);
        write_h5_numeric("minus_alt", itmp, file, PredType::NATIVE_UINT32);
    }

}

void DataManager::run() {
    using namespace std::placeholders;
    IndexProcessor::getrange f = std::bind(&DataManager::get_range, &(*this), _1);
    for(size_t i = 0; i < processors.size(); i++){
        processors[i]->start(f, groups);
    }
    for(size_t i = 0; i < processors.size(); i++){
        processors[i]->join();
    }

    for(auto & p : processors){
        calls.insert(calls.end(), p->calls.begin(), p->calls.end());
    }
    std::sort(calls.begin(), calls.end());
}
