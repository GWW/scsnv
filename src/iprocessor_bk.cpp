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

using namespace gwsc;

void IndexProcessor::run_position(const std::string & chrom, unsigned int pos){
    {
        std::lock_guard<std::mutex> lk(mtx_wait_);
        chrom_ = chrom;
        pos_ = pos;
        ready_ = true;
        processed_ = false;
    }
    cv_.notify_one();
}

void IndexProcessor::wait(){
    std::unique_lock<std::mutex> lk(mtx_wait_);
    while(!processed_) cv_.wait(lk, [this]{return this->processed_;});
}

void IndexProcessor::done(){
    std::unique_lock<std::mutex> lk(mtx_wait_);
    done_ = true;
    ready_ = true;
    processed_ = false;
    lk.unlock();
    cv_.notify_one();
}

//void IndexProcessor::pileup_position(){
//}

void IndexProcessor::operator()(){
    while(!done_){
        std::unique_lock<std::mutex> lk(mtx_wait_);
        while(!ready_) cv_.wait(lk, [this]{ return ready_; });
        if(done_) break;

        tid_ = bam_name2id(bh_, chrom_.c_str());
        curr_  = 0;
        read_ = 0;
        valid_ = false;
        if(tid_ == -1) return;
        hts_itr_t *iter = sam_itr_queryi(bi_, tid_, pos_, pos_ + 1);
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
            pileup_.reset(tid_, pos_);
            pileup_.process_range(rpair);
            for(auto it = pileup_.positions.begin(); it != pileup_.positions.end(); it++){
                if(it->pos == pos_){
                    valid_ = true;
                    entry_ = it;
                    break;
                }
            }

        }
        processed_ = true;
        ready_ = false;
        lk.unlock();
        cv_.notify_one();
    }
}

void IndexProcessor::header(gzofstream & zout) const {
    if(filter_barcodes_){
        zout 
            << "\t" << name_ << "_plus_ref"
            << "\t" << name_ << "_plus_alt"
            << "\t" << name_ << "_minus_ref"
            << "\t" << name_ << "_minus_alt"
            << "\t" << name_ << "_ref_barcodes"
            << "\t" << name_ << "_alt_barcodes"
            << "\t" << name_ << "_plus_ref_barcodes"
            << "\t" << name_ << "_plus_alt_barcodes"
            << "\t" << name_ << "_plus_altref_barcodes"
            << "\t" << name_ << "_minus_ref_barcodes"
            << "\t" << name_ << "_minus_alt_barcodes"
            << "\t" << name_ << "_minus_altref_barcodes"
            << "\t" << name_ << "_total_altref_barcodes"
            << "\t" << name_ << "_total_barcodes";
    }else{
        zout 
            << "\t" << name_ << "_ref"
            << "\t" << name_ << "_alf";
    }

}

void IndexProcessor::write(gzofstream & zout, unsigned int refi, unsigned int alti) const {
    if(filter_barcodes_){
        if(valid()){
            auto e = *entry();
            zout 
                << "\t" << e.bases[refi].p_count
                << "\t" << e.bases[alti].p_count
                << "\t" << e.bases[refi].m_count
                << "\t" << e.bases[alti].m_count
                << "\t" << e.bases[refi].t_barcodes
                << "\t" << e.bases[alti].t_barcodes

                << "\t" << e.bases[refi].p_barcodes
                << "\t" << e.bases[alti].p_barcodes
                << "\t" << e.bases[alti].par_barcodes

                << "\t" << e.bases[refi].m_barcodes
                << "\t" << e.bases[alti].m_barcodes
                << "\t" << e.bases[alti].mar_barcodes

                << "\t" << e.bases[alti].tar_barcodes
                << "\t" << e.barcodes
                ;
        }else{
            zout << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
        }
    }else{
        if(valid()){
            auto e = *entry();
            zout 
                << "\t" << (e.bases[refi].m_count + e.bases[refi].p_count)
                << "\t" << (e.bases[alti].m_count + e.bases[alti].p_count);
        }else{
            zout << "\t0\t0";
        }
    }
}
