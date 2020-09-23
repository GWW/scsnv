#pragma once
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

#include <vector>
#include <array>
#include "index.hpp"

namespace gwsc {

struct SpliceLocations {

    using sit = std::vector<unsigned int>::const_iterator;
    SpliceLocations() {

    }
    SpliceLocations(const TXIndex & index) : index(&index) {
    }

    int find_before_(sit & start, sit end, int pos){
        int d = sdist + 1;
        while(start != end && static_cast<int>(*start) < (pos - sdist)){
            start++;
        }
        sit it = start;
        while(it != end && static_cast<int>(*it) <= pos){
            int delta = pos - *it;
            d = std::min(delta, d);
            it++;
        }
        return d;
    }

    int find_after_(sit & start, sit end, int pos){
        int d = sdist + 1;

        while(start != end && static_cast<int>(*start) < pos){
            start++;
        }

        sit it = start;
        while(it != end && static_cast<int>(*it) <= (pos + sdist)){
            int delta = *it - pos;
            d = std::min(delta, d);
            it++;
        }
        return d;
    }

    std::array<int, 4> closest(int tid, unsigned int pos){
        if(tid < 0) return {sdist + 1, sdist + 1, sdist + 1, sdist + 1};
        if(cref == nullptr || cref->tid != static_cast<unsigned int>(tid)) {
            cref = &index->ref(tid);
            pl_start = cref->plus_lsplices.begin();
            pr_start = cref->plus_rsplices.begin();
            ml_start = cref->minus_lsplices.begin();
            mr_start = cref->minus_rsplices.begin();
        }

        int pl_dist = sdist + 1;
        int pr_dist = sdist + 1;
        int ml_dist = sdist + 1;
        int mr_dist = sdist + 1;
        pl_dist = find_before_(pl_start, cref->plus_lsplices.end(), pos);
        ml_dist = find_before_(ml_start, cref->minus_lsplices.end(), pos);
        pr_dist = find_after_(pr_start, cref->plus_rsplices.end(), pos);
        mr_dist = find_after_(mr_start, cref->minus_rsplices.end(), pos);
        //pa_dist = dist_(pl_start, cref->p_acceptors.end(), pos);
        //pd_dist = dist_(pd_start, cref->p_donors.end(), pos, false);
        //ma_dist = dist_(ma_start, cref->m_acceptors.end(), pos, false);
        //md_dist = dist_(md_start, cref->m_donors.end(), pos, true);
        //std::cout << "  Dists = " << pa_dist << ", " << pd_dist << ", " << ma_dist << ", " << md_dist << "\n";
        return {pl_dist, pr_dist, ml_dist, mr_dist};
    }

    const TXIndex      * index = nullptr;
    const TXIndex::Ref * cref = nullptr;

    sit                  pl_start;
    sit                  pr_start;
    sit                  ml_start;
    sit                  mr_start;
    int                  sdist = 10;
};

class TargetFinder {
    public:
        struct Target{
            Target() {

            }

            Target(unsigned int target, char base)
                : target(target), base(base)
            {

            }

            bool operator<(const Target & rhs) const {
                return target < rhs.target;
            }

            unsigned int target;
            char         base = 'N';
        };

        using tvect = std::vector<Target>;
        using trefs = std::vector<tvect>;

        TargetFinder(){

        }

        void set_refs(const trefs & refs){
            refs_ = &refs;
        }

        bool loaded() const {
            return refs_ != nullptr;
        }

        void rewind(int32_t tid, uint32_t pos){
            ltid_ = tid;
            lindex_ = 0;
            if(tid >= static_cast<int>(refs_->size())){
                return;
            }
            auto it = std::lower_bound(refs_->at(tid).begin(), refs_->at(tid).end(), pos, [](const Target & t, uint32_t pos) -> bool { return t.target < pos; } );
            lindex_ = it - refs_->at(tid).begin();
        }

        const Target * check(int32_t tid, uint32_t pos) {
            if(refs_ == nullptr || tid < 0 || tid >= static_cast<int>(refs_->size())){
                return nullptr;
            }
            if(ltid_ != tid){
                ltid_ = tid;
                lindex_ = 0;
            }

            auto & tv = refs_->at(ltid_);
            while(lindex_ < tv.size() && tv[lindex_].target < pos){
                lindex_++;
            }
            if(lindex_ < tv.size() && tv[lindex_].target == pos){
                return &tv[lindex_];
            }
            return nullptr;
        }

    private:
        const trefs        * refs_ = nullptr;
        size_t               lindex_ = 0;
        int32_t              ltid_ = -1;
};


struct PositionCoverage {
    PositionCoverage(int tid = -1, unsigned int pos = 0, unsigned int pcoverage = 0, unsigned int mcoverage = 0, 
            unsigned int tbarcodes = 0, unsigned int pbarcodes = 0, unsigned int mbarcodes = 0) 
        :tid(tid), pos(pos), pcoverage(pcoverage), mcoverage(mcoverage), tbarcodes(tbarcodes), pbarcodes(pbarcodes), mbarcodes(mbarcodes)
    {
    }

    bool operator<(const PositionCoverage & p) const {
        return tid < p.tid || (tid == p.tid && pos < p.pos);
    }
    int32_t tid = -1;
    uint32_t pos = 0;
    uint32_t pcoverage = 0;
    uint32_t mcoverage = 0;
    uint32_t tbarcodes = 0;
    uint32_t pbarcodes = 0;
    uint32_t mbarcodes = 0;
};

struct BarcodeCount{
    BarcodeCount(uint32_t barcode)
        : barcode(barcode)
    {
        for(size_t i = 0; i < 4; i++) pbases[i] = mbases[i] = 0;
    }

    void inc_base(uint8_t b, uint32_t count, bool rev) {
        if(rev){
            mbases[b] = std::min((uint32_t)std::numeric_limits<uint16_t>::max(), (uint32_t)mbases[b] + count);
        }else{
            pbases[b] = std::min((uint32_t)std::numeric_limits<uint16_t>::max(), (uint32_t)pbases[b] + count);
        }
    }

    uint32_t barcode;
    uint16_t pbases[4];
    uint16_t mbases[4];
};

struct BarcodeRate {
    static const size_t BTOTAL = 24;
    BarcodeRate()
    {
        for(size_t i = 0; i < BTOTAL; i++) prates[i] = mrates[i] = 0;
    }

    BarcodeRate & operator+=(const BarcodeRate & rhs){
        for(size_t i = 0; i < BTOTAL; i++) {
            prates[i] += rhs.prates[i];
            mrates[i] += rhs.mrates[i];
        }
        return *this;
    }
    uint32_t prates[BTOTAL];
    uint32_t mrates[BTOTAL];
};

struct PositionCount{
    void reset() {
        bases = {};
        tid = -1;
        pos = 0;
        coverage = 0;
        barcodes = 0;
        ambig = 0;
        ref = 'N';
        refi = 4;
        pbase = '-';
        mbase = '-';
        max_nr = 'N';
        bcounts.clear();
    }

    bool operator<(const PositionCount & p) const {
        return tid < p.tid || (tid == p.tid && pos < p.pos);
    }

    struct BaseCount{
        uint32_t         p_count = 0;
        uint32_t         m_count = 0;
        uint32_t         p_edge_dist = 0;
        uint32_t         m_edge_dist = 0;
        uint32_t         t_barcodes = 0;
        uint32_t         m_barcodes = 0;
        uint32_t         p_barcodes = 0;
        uint32_t         par_barcodes = 0;
        uint32_t         mar_barcodes = 0;
        uint32_t         tar_barcodes = 0;
    };

    std::array<BaseCount, 4>  bases;
    std::array<int, 4>        sdists;
    std::vector<BarcodeCount> bcounts;
    int                       tid = -1;
    unsigned int              pos = 0;
    unsigned int              coverage = 0;
    unsigned int              barcodes = 0;
    unsigned int              pbarcodes = 0;
    unsigned int              mbarcodes = 0;
    unsigned int              ambig = 0;
    char                      ref = 'N';
    char                      pbase = '-';
    char                      mbase = '-';
    char                      max_nr = 'N';
    uint8_t                   refi = 4;
};

}
