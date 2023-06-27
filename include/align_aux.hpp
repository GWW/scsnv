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

#include <string>
#include <cstdint>
#include <vector>
#include <ostream>
#include <sstream>
#include <limits>
#include "interval_tree.hpp"
#include "bwamem.h"
#include "htslib/htslib/sam.h"

namespace gwsc {

enum StrandMode {
    TAG_UNKNOWN = 0,
    TAG_FWD = 1, // Tag maps to the forward strand // For 10X 3' R2 maps to the forward strand of the transcript
    TAG_REV = 2  // Tag maps to the reverse strand // For 10X 5' R2 maps to the reverse strand of the transcript
};

static const char * CIGAR_CHARS = (char*)"MIDNSHP=XB";

enum class Cigar : std::uint8_t {
    MATCH,
    INS,
    DEL,
    REF_SKIP,
    SOFT_CLIP,
    HARD_CLIP,
    PAD,
    EQUAL,
    DIFF,
    BACK
};

struct CigarElement {
    CigarElement(uint32_t cigar){
        unpack(cigar);
    }

    CigarElement(uint32_t len, Cigar op) : len(len), op(op) {

    }

    CigarElement() : len(0), op(Cigar::PAD) {

    }


    bool operator==(const CigarElement & c) const {
        return c.len == len && op == c.op;
    }

    bool operator!=(const CigarElement & c) const {
        return c.len == len && op == c.op;
    }

    uint32_t packed() const {
        return (len << 4) + static_cast<uint8_t>(op);
    }

    void unpack(uint32_t v) {
        len = v >> 4;
        op = static_cast<Cigar>(v & 0xF);
    }

    friend std::ostream & operator<< (std::ostream &os, const CigarElement &e) {
        os << (e.len) << CIGAR_CHARS[static_cast<std::size_t>(e.op)];
        return os;
    }

    uint32_t len;
    Cigar    op;
};

class CigarString {
    public:
        typedef std::vector<CigarElement> ElementVect;
        typedef ElementVect::const_iterator const_iterator;
        typedef ElementVect::iterator       iterator;

        typedef ElementVect::const_reverse_iterator const_reverse_iterator;
        typedef ElementVect::reverse_iterator       reverse_iterator;

        CigarString() { }

        void clear() {
            cigar_.clear();
        }

        bool operator==(const CigarString & c) const {
            return cigar_ == c.cigar_;
        }

	bool operator!=(const CigarString & c) const {
	    return !(*this == c);
	}

        size_t size() const {
            return cigar_.size();
        }

	bool empty() const {
	    return cigar_.empty();
	}

        void resize(size_t s, const CigarElement &e) {
            cigar_.resize(s, e);
        }

        void append(uint32_t cigar) {
            cigar_.push_back(CigarElement(cigar));
        }

        void append(const CigarElement &c) {
            cigar_.push_back(c);
        }

        void append(uint32_t len, Cigar op) {
            cigar_.push_back(CigarElement(len, op));
        }

        void push_back(uint32_t cigar) {
            cigar_.push_back(CigarElement(cigar));
        }

        void push_back(const CigarElement &c) {
            cigar_.push_back(c);
        }

        void push_back(uint32_t len, Cigar op) {
            cigar_.push_back(CigarElement(len, op));
        }

        void pop_back() {
            cigar_.pop_back();
        }

        iterator insert(iterator pos, uint32_t cigar) {
            return cigar_.insert(pos,CigarElement(cigar));
        }

        iterator insert(iterator pos, uint32_t len, Cigar op) {
            return cigar_.insert(pos, CigarElement(len, op));
        }

        iterator insert(iterator pos, const CigarElement &c) {
            return cigar_.insert(pos,CigarElement(c));
        }

        CigarElement & front() {
            return cigar_.front();
        }

        const CigarElement & front() const {
            return cigar_.front();
        }

        CigarElement & back() {
            return cigar_.back();
        }

        const CigarElement & back() const {
            return cigar_.back();
        }

        iterator erase(iterator pos) {
            return cigar_.erase(pos);
        }

        iterator erase(iterator start, iterator end) {
            return cigar_.erase(start,end);
        }

        iterator begin() {
            return cigar_.begin();
        }

        const_iterator begin() const {
            return cigar_.begin();
        }

        iterator  end() {
            return cigar_.end();
        }

        const_iterator end() const {
            return cigar_.end();
        }

        reverse_iterator  rbegin() {
            return cigar_.rbegin();
        }

        const_reverse_iterator rbegin() const {
            return cigar_.rbegin();
        }

        reverse_iterator  rend() {
            return cigar_.rend();
        }

        const_reverse_iterator rend() const {
            return cigar_.rend();
        }

        friend std::ostream& operator<< (std::ostream &os, const CigarString &c) {
            for(const_iterator it = c.begin(); it != c.end(); ++it)
                os << (*it);
            return os;
        }

        unsigned int rbases() const {
            unsigned int b = 0;
            for(auto c : cigar_){
                if(c.op == Cigar::MATCH || c.op == Cigar::DEL || c.op == Cigar::REF_SKIP) b += c.len;
            }
            return b;
        }

        unsigned int qbases() const {
            unsigned int b = 0;
            auto it = begin();
            if(it->op == Cigar::SOFT_CLIP){
                b = it->len;
                it++;
            }

            while(it != end()){
                if(it->op == Cigar::MATCH || it->op == Cigar::INS) b += it->len;
            }
            return b;
        }

    private:
        ElementVect cigar_;
};


union BamFlag {
    struct {
        uint16_t paired:1;
        uint16_t proper_pair:1;
        uint16_t unmapped:1;
        uint16_t m_unmapped:1;
        uint16_t strand:1;
        uint16_t m_strand:1;
        uint16_t read1:1;
        uint16_t read2:1;
        uint16_t secondary:1;
        uint16_t qc_fail:1;
        uint16_t duplicate:1;
        uint16_t supplementary:1;
        uint16_t _pad:4;
    }f;
    uint16_t flag_val;
};


/*
static const char * FLAG_NAMES[] = {
    "PAIRED",
    "PROPER_PAIR",
    "UNMAP",
    "MUNMAP",
    "REVERSE",
    "MREVERSE",
    "READ1",
    "READ2",
    "SECONDARY",
    "QCFAIL",
    "DUP",
    "SUPPLEMENTARY"
};
*/

enum AlignType {
    UNKNOWN=0,
    CDNA = 1,
    INTRON = 2,
    INTERGENIC=3,
    AMBIGUOUS=4,
    ANTISENSE=5,
    ELEM_COUNT=6
};

inline const char * aligntype2str(AlignType a){
    switch(a){
        case AlignType::UNKNOWN:
            return "Unknown";
        case AlignType::CDNA:
            return "Exonic";
        case AlignType::INTRON:
            return "Intronic";
        case AlignType::INTERGENIC:
            return "Intergenic";
        case AlignType::AMBIGUOUS:
            return "Ambiguous";
        case AlignType::ANTISENSE:
            return "Antisense";
        case AlignType::ELEM_COUNT:
            return "";
    }
    return "";
}

inline char aligntype2code(AlignType a){
    switch(a){
        case AlignType::UNKNOWN:
            return '?';
        case AlignType::CDNA:
            return 'E';
        case AlignType::INTRON:
            return 'N';
        case AlignType::INTERGENIC:
            return 'I';
        case AlignType::AMBIGUOUS:
            return '?';
        case AlignType::ANTISENSE:
            return 'A';
        case AlignType::ELEM_COUNT:
            return '?';
    }
    return '?';
}


// Only countable alignments
struct AlignSummary{
    using bint = uint32_t;
    uint64_t pos;
    bint     barcode:31;
    bint     intronic:1;
    uint32_t gene_id;
    uint32_t umi;        // All the UMI's I have seen so far are less than 16 bases long

    AlignSummary(){

    }

    AlignSummary(AlignSummary::bint barcode, uint32_t gene_id, uint32_t umi, bool intronic, uint64_t pos) 
        : pos(pos), barcode(barcode), intronic(intronic), gene_id(gene_id), umi(umi)
    {

    }
    void reset() {
        barcode = 0;
        gene_id = 0;
        umi = 0;
        pos = 0;
        intronic = false;
    }

    bool operator<(const AlignSummary & rhs) const {
        return std::tie(barcode, gene_id, umi) < std::tie(rhs.barcode, rhs.gene_id, rhs.umi);
    }
};

struct AlignTagOut{
    AlignTagOut() 
        : barcode(0), intronic(false), gene_id(0), umi(0)
    {

    }

    AlignTagOut(AlignSummary::bint barcode, uint32_t gene_id, uint32_t umi, bool intronic) 
        :barcode(barcode), intronic(intronic), gene_id(gene_id), umi(umi)
    {

    }

    AlignSummary::bint barcode:31;
    AlignSummary::bint intronic:1;
    uint32_t           gene_id;
    uint32_t           umi;
};

struct AlignData {
    void reset(){
        tid = -1;
        total_bases = 0;
        dels = 0;
        NM = 0;
        mapq = 0;
        lft = rgt = 0;
        txid = std::numeric_limits<uint32_t>::max();
        gid = std::numeric_limits<uint32_t>::max();
        gids.clear();
        intergenic = 0;
        intronic = 0;
        exonic = 0;
        xs = '?';
        atype = AlignType::UNKNOWN;
        cigar.clear();
        rescore = false;
        passed = false;
        valid_transcript = false;
        rev = false;
    }

    std::string print(const std::string & rname = "") const {
        std::stringstream ss;
        ss << rname << ": " << lft << "-" << rgt << " rev = " << rev << " xs = " << xs 
           << " cigar: " << cigar << " type: " << aligntype2str(atype) << " tid = " << tid;
        if(atype != AlignType::UNKNOWN && atype != AlignType::INTERGENIC){
            ss << " gid = " << gid << " txid = " << txid;
        }
        ss  << " exonic = " << exonic << " intergenic = " << intergenic << " intronic = " << intronic << " total bases = " << total_bases 
            << " NM = " << NM << " gaps = " << gaps << " score = " << score; //<< " MD = " << MD;
        return ss.str();
    }

    unsigned int clipped(){
        return (cigar.front().op == Cigar::SOFT_CLIP ? cigar.front().len : 0) + 
               (cigar.back().op == Cigar::SOFT_CLIP ? cigar.back().len : 0);
    }

    bool operator<(const AlignData & a) const {
        return std::make_tuple(atype == AlignType::ANTISENSE, gid, -1 * static_cast<int>(score), txid) < 
            std::make_tuple(a.atype == AlignType::ANTISENSE, a.gid, -1 * static_cast<int>(a.score), a.txid);
    }

    CigarString  cigar;
    std::vector<uint32_t> gids;
    int32_t      tid;
    uint32_t     gid;
    uint32_t     txid;
    int32_t      score;
    uint32_t     lft;
    uint32_t     rgt;
    uint32_t     total_bases;
    uint32_t     dels; // Including deletions
    uint32_t     intergenic;
    uint32_t     intronic;
    uint32_t     exonic;
    uint32_t     gaps;
    uint32_t     NM;
    uint8_t      mapq;
    AlignType    atype;
    char         xs;          //Strand the read is derived from
    bool         rev;         //Alignment is reverse
    bool         valid_transcript;
    bool         rescore;
    bool         passed;
};

struct AlignScore {
    unsigned int match = 1;
    unsigned int mismatch = 4;
    unsigned int max_diff = 3; 
    unsigned int max_diff_search = 20;
    unsigned int min_score = 45; // will also limit the number of aligned bases
};

struct AlignBases {
    unsigned int exonic = 0;
    unsigned int intronic = 0;
    unsigned int intergenic = 0;
    bool         antisense = false;
};

struct AlignGroup {
    enum Result {
        BARCODE_FAIL,
        UMI_FAIL,
        TAG_FAIL,
        UNMAPPED,
        INTERGENIC,
        CDNA,
        INTRONIC,
        MULTIMAPPED,
        AMBIGUOUS,
        ANTISENSE,
        ELEM_COUNT
    };

    static const char * alignres2str(Result a, bool nice = false){
        switch(a){
            case Result::BARCODE_FAIL:
                return nice ? "Barcode QA Fail" : "barcode_qa_fail";
            case Result::UMI_FAIL:
                return nice ? "UMI QA Fail" : "umi_qa_fail";
            case Result::TAG_FAIL:
                return nice ? "TAG QA Fail" : "tag_qa_fail";
            case Result::UNMAPPED:
                return nice ? "Unmapped" : "unmapped";
            case Result::INTERGENIC:
                return nice ? "Intergenic" : "intergenic";
            case Result::CDNA:
                return nice ? "cDNA" : "cdna";
            case Result::INTRONIC:
                return nice ? "Intronic" : "intronic";
            case Result::MULTIMAPPED:
                return nice ? "Multimapped" : "multimapped";
            case Result::AMBIGUOUS:
                return nice ? "Ambiguous" : "ambiguous";
            case Result::ANTISENSE:
                return nice ? "Antisense" : "antisense";
            case Result::ELEM_COUNT:
                return "";
        }
        return "";
    }

    static char alignres2code(Result a){
        switch(a){
            case Result::BARCODE_FAIL:
                return 'B';
            case Result::UMI_FAIL:
                return 'U';
            case Result::TAG_FAIL:
                return 'T';
            case Result::UNMAPPED:
                return 'U';
            case Result::INTERGENIC:
                return 'I';
            case Result::CDNA:
                return 'E';
            case Result::INTRONIC:
                return 'N';
            case Result::MULTIMAPPED:
                return 'M';
            case Result::AMBIGUOUS:
                return '?';
            case Result::ANTISENSE:
                return 'A';
            case Result::ELEM_COUNT:
                return '?';
        }
        return '?';
    }

    using ResultCounts = std::array<unsigned int, AlignGroup::Result::ELEM_COUNT>;

    AlignGroup() {
        alns.resize(10);
        transcript_ar.a = NULL;
        transcript_ar.n = 0;
        genome_ar.a = NULL;
        genome_ar.n = 0;
    }

    ~AlignGroup(){
    }

    void reset(){
        if(genome_ar.a != NULL) free(genome_ar.a);
        genome_ar.a = NULL;
        genome_ar.n = 0;
        if(transcript_ar.a != NULL) free(transcript_ar.a);
        transcript_ar.a = NULL;
        transcript_ar.n = 0;
        summary.reset();
        res = UNMAPPED;
        acount = 0;
        transcript_score = 0;
        genome_score = 0;
        countable = false;
        total_genes = 0;
    }

    std::vector<Interval<unsigned int, unsigned int>> toverlaps;
    std::vector<AlignData>      alns;
    std::vector<AlignBases>     tbases;
    CigarString                 tcigar;   // for various operations to keep things thread safe and avoid reallocation
    std::string                 tmd;   // for various operations to keep things thread safe and avoid reallocation
    AlignSummary                summary;
    Result                      res;                       
    mem_alnreg_v                transcript_ar;
    mem_alnreg_v                genome_ar;
    size_t                      acount = 0;
    int                         transcript_score = 0;
    int                         genome_score = 0;
    unsigned int                total_genes = 0;
    unsigned int                total_passed = 0;
    bool                        countable = false;
};

struct UMIMap {
    UMIMap(AlignSummary::bint barcode, uint32_t gene_id, uint32_t umi_from, uint32_t umi_to) 
        : barcode(barcode), gene_id(gene_id), umi_from(umi_from), umi_to(umi_to)
    {

    }

    UMIMap() : barcode(0), gene_id(0), umi_from(0), umi_to(0){

    }

    bool operator<(const UMIMap & m) const {
        return
            std::tie(barcode, gene_id, umi_from) <
            std::tie(m.barcode, m.gene_id, m.umi_from);
    }

    bool operator>(const UMIMap & m) const {
        return
            std::tie(barcode, gene_id, umi_from) >
            std::tie(m.barcode, m.gene_id, m.umi_from);
    }

    bool operator==(const UMIMap & m) const {
        return
            std::tie(barcode, gene_id, umi_from) ==
            std::tie(m.barcode, m.gene_id, m.umi_from);
    }

    std::size_t hash() const {
        std::hash<uint32_t> hasher;
        std::size_t seed = 0;
        seed ^= hasher(barcode) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(gene_id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(umi_from) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }

    AlignSummary::bint barcode;
    uint32_t           gene_id;
    uint32_t           umi_from;
    uint32_t           umi_to;
};

struct UMIBad {
    UMIBad(AlignSummary::bint barcode, uint32_t gene_id, uint32_t umi) 
        : barcode(barcode), gene_id(gene_id), umi(umi)
    {

    }

    UMIBad() : barcode(0), gene_id(0), umi(0){

    }

    bool operator<(const UMIBad & m) const {
        return
            std::tie(barcode, gene_id, umi) <
            std::tie(m.barcode, m.gene_id, m.umi);
    }

    bool operator>(const UMIBad & m) const {
        return
            std::tie(barcode, gene_id, umi) >
            std::tie(m.barcode, m.gene_id, m.umi);
    }

    bool operator==(const UMIBad & m) const {
        return
            std::tie(barcode, gene_id, umi) ==
            std::tie(m.barcode, m.gene_id, m.umi);
    }

    bool operator<(const AlignSummary & m) const {
        return
            std::tie(barcode, gene_id, umi) <
            std::tie(m.barcode, m.gene_id, m.umi);
    }

    bool operator==(const AlignSummary & m) const {
        return
            std::tie(barcode, gene_id, umi) ==
            std::tie(m.barcode, m.gene_id, m.umi);
    }

    AlignSummary::bint barcode;
    uint32_t           gene_id;
    uint32_t           umi;
};

struct UMIMapHasher {
    std::size_t operator()(const UMIMap & m) const {
        return m.hash();
    }
};

struct GeneCount {
    GeneCount(AlignSummary::bint barcode, uint32_t gene_id, uint32_t molecules, uint32_t intronic) 
        : barcode(barcode), gene_id(gene_id), molecules(molecules), intronic(intronic)
    {

    }

    GeneCount() : barcode(0), gene_id(0), molecules(0), intronic(0){

    }
    AlignSummary::bint barcode;
    uint32_t           gene_id;
    uint32_t           molecules;
    uint32_t           intronic;
};

}
