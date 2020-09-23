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

#include "sbam_writer.hpp"
#include <iomanip>
#include <sstream>

using namespace gwsc;

void SortedBamWriter::write(const AlignGroup & g, const Read & r, read_buffer & buffer, unsigned int & rcount, const TXIndex & idx){
    bam1_t * s = NULL;
    if(rcount < buffer.size()){
        s = buffer[rcount++];
    }else{
        std::lock_guard<std::mutex> lock(mtx_write_);
        if((buffer.size() + file_reads_) >= out_sz_){
            // Must write the contents of the buffer
            force_write();
        }
        //std::cout << "Filling output slots " << file_reads_ << " - " << (file_reads_ + buffer.size()) << " out size = " << out_.size() << " buffer size = " << buffer.size() << " rcount = " << rcount << "\n";
        for(size_t i = 0; i < buffer.size(); i++, file_reads_++){
            std::swap(buffer[i], out_[file_reads_]);
        }
        rcount = 0;
        s = buffer[rcount++];
        if(s == NULL) {
            std::cout << "Couldn't get a new bam record from a flushed buffer\n";
            exit(1);
        }
    }

    if(g.res == AlignGroup::AMBIGUOUS || g.res == AlignGroup::UNMAPPED || 
            g.res == AlignGroup::BARCODE_FAIL || g.res == AlignGroup::BARCODE_FAIL ||
            g.res == AlignGroup::TAG_FAIL || g.res == AlignGroup::UMI_FAIL || 
            (g.res == AlignGroup::MULTIMAPPED && g.total_passed > 5))
    { 
        _align2unmapped(s, g, r);
    }else{
        _align2bam(s, g, r, idx);
    }
}

void SortedBamWriter::merge_buffer(read_buffer & buffer, unsigned int & rcount) {
    // For the last reads just put them all in the last file even if it goes over the required size
    for(size_t i = 0; i < rcount; i++, file_reads_++){
        if(file_reads_ >= out_.size()){
            out_.push_back(bam_init1());
        }
        std::swap(out_[file_reads_], buffer[i]);
    }
    rcount = 0;
}

void SortedBamWriter::force_write() {
    if(file_reads_ == 0) return;
    SortBamTidPos psort;
    psort.max_tid = bh.bam_hdr()->n_targets;

    std::sort(out_.begin(), out_.begin() + file_reads_, psort);
    std::stringstream ss;
    ss << prefix_  << "/scsnv_tmp_" << std::setw(4) << std::setfill('0') << file_number_ << ".bam";
    std::string outf = ss.str();

    samFile * bam_out = sam_open(outf.c_str(), "wb");
    if(pool_threads_ > 1) hts_set_thread_pool(bam_out, &pool_);
    if(sam_hdr_write(bam_out, bh.bam_hdr()) < 0) {
        std::cerr << "Error writing header\n";
        exit(1);
    }
    for(size_t i = 0; i < file_reads_; i++){
        if(sam_write1(bam_out, bh.bam_hdr(), out_[i]) < 0){
            std::cerr << "Error writing sam\n";
            exit(1);
        }
    }
    sam_close(bam_out);
    bam_out = nullptr;
    total_reads_ += file_reads_;
    file_reads_ = 0;
    file_number_++;
}

void SortedBamWriter::_align2unmapped(bam1_t * bam, const AlignGroup & g, const Read & r){
    BamFlag flag;
    flag.f.unmapped = true;
    bam->core.flag = flag.flag_val;
    bam->core.tid = -1;
    bam->core.mtid = -1;
    bam->core.mpos = -1;
    bam->core.isize = -1;
    bam->core.pos = 0;
    bam->core.bin = hts_reg2bin(0, 1, 14, 5);

    bam->core.n_cigar = 0;
    auto qlen = r.name.length() + 1;
    bam->core.l_extranul = (4 - (qlen & 3)) & 3;

    bam->core.l_qname = qlen + bam->core.l_extranul; // Because of the null terminator
    bam->core.l_qseq  = r.tag.length();

    uint32_t dl = bam->core.n_cigar*4 + bam->core.l_qname + bam->core.l_qseq + (bam->core.l_qseq + 1)/2;

    if(bam->m_data < dl) {
        bam->m_data = dl;
        kroundup32(bam->m_data);
        bam->data = (uint8_t*)realloc(bam->data, bam->m_data);
    }

    bam->l_data = dl;
    // Copy the qname
    uint8_t * p = bam->data;
    memcpy(p, r.name.c_str(), qlen);
    for(size_t i = 0; i < bam->core.l_extranul; i++){
        p[qlen + i] = 0;
    }
    p += bam->core.l_qname;

    for(size_t i = 0; i < r.tag.length(); i++){
        uint8_t * pp = p + (i / 2);
        uint8_t mask = 0xFU << (((~i) & 1) << 2);
        *pp = ((*pp) & ~mask) | seq_nt16_table[static_cast<size_t>(r.tag[i])] << (((~i) & 1) << 2);
    }

    p += (bam->core.l_qseq + 1)/2;
    for(size_t i = 0; i < r.q_tag.length(); i++){
        *p++ = r.q_tag[i] - 33;
    }

    if((p - bam->data) != dl){
        std::cout << "Data length wrong dl = " << dl << " compared to " << (p - bam->data) << "\n";
        exit(1);
    }

    //uint32_t ol = bam->m_data;
    bam_aux_append(bam, "UB", 'Z', r.umi.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(r.umi.c_str())));
    bam_aux_append(bam, "CB", 'Z', r.barcode.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(r.barcode.c_str())));

    char tmpc = AlignGroup::alignres2code(g.res);
    bam_aux_append(bam, "RA", 'A', 1, reinterpret_cast<uint8_t*>(&tmpc));
}

void SortedBamWriter::_align2bam(bam1_t * bam, const AlignGroup & g, const Read & r, const TXIndex & idx){
    auto it = g.alns.begin();
    auto aend = g.alns.begin() + g.acount;
    if(g.countable){
        // Make sure the alignment is to the gene that is optimal
        while(it != aend && !it->passed && it->gid != g.summary.gene_id){
            it++;
        }
        if(it == aend){
            std::cout << "Could not find alignment for matching gene id\n";
            exit(1);
        }
    }else{
        while(it != aend && !it->passed){
            it++;
        }
        if(it == aend){
            std::cout << "Write Read " << AlignGroup::alignres2str(g.res) << " acount = " << g.acount << "\n";
            for(auto dit = g.alns.begin(); dit != g.alns.end(); dit++){
                std::cout << "  P = " << dit->passed << " " << dit->print() << "\n";
            }
            std::cout << "Could not find alignment with a passing score?\n";
            exit(1);
        }
    }

    auto & a = *it;
    //std::cout << "  " << a.print() << "\n";
    BamFlag flag;
    //flag.f.read2 = true;
    flag.f.strand = a.rev;

    bam->core.flag = flag.flag_val;
    if(g.res == AlignGroup::MULTIMAPPED){
        bam->core.qual = 0; //a.mapq;  All "uniquely" mapped readsk
    }else{
        bam->core.qual = 255; //a.mapq;  All "uniquely" mapped reads
    }
    bam->core.n_cigar = a.cigar.size();
    bam->core.tid = a.tid;
    bam->core.mtid = -1;
    bam->core.mpos = -1;
    bam->core.isize = -1;
    bam->core.pos = a.lft;

    bam->core.bin = hts_reg2bin(a.lft, a.rgt + 1, 14, 5);

    auto qlen = r.name.length() + 1;
    bam->core.l_extranul = (4 - (qlen & 3)) & 3;

    bam->core.l_qname = qlen + bam->core.l_extranul; // Because of the null terminator
    bam->core.l_qseq  = r.tag.length();

    uint32_t dl = bam->core.n_cigar*4 + bam->core.l_qname + bam->core.l_qseq + (bam->core.l_qseq + 1)/2;


    if(bam->m_data < dl) {
        bam->m_data = dl;
        kroundup32(bam->m_data);
        bam->data = (uint8_t*)realloc(bam->data, bam->m_data);
    }

    bam->l_data = dl;
    // Copy the qname
    uint8_t * p = bam->data;
    memcpy(p, r.name.c_str(), qlen);
    for(size_t i = 0; i < bam->core.l_extranul; i++){
        p[qlen + i] = 0;
    }
    p += bam->core.l_qname;

    // copy the cigar
    for(auto c : a.cigar){
        uint32_t v = c.packed();
        memcpy(p, &v, 4);
        p+=4;
    }

    if(!flag.f.strand){
        for(size_t i = 0; i < r.tag.length(); i++){
            uint8_t * pp = p + (i / 2);
            uint8_t mask = 0xFU << (((~i) & 1) << 2);
            *pp = ((*pp) & ~mask) | seq_nt16_table[static_cast<size_t>(r.tag[i])] << (((~i) & 1) << 2);
        }

        p += (bam->core.l_qseq + 1)/2;
        for(size_t i = 0; i < r.q_tag.length(); i++){
            *p++ = r.q_tag[i] - 33;
        }
    }else{
        for(size_t i = 0; i < r.tag.length(); i++){
            uint8_t * pp = p + (i / 2);
            uint8_t mask = 0xFU << (((~i) & 1) << 2);
            char rev = reverse_cmpl_[static_cast<size_t>(r.tag[r.tag.length() - i - 1])];
            *pp = ((*pp) & ~mask) | seq_nt16_table[static_cast<size_t>(rev)] << (((~i) & 1) << 2);
        }
        p += (bam->core.l_qseq + 1)/2;
        for(size_t i = 0; i < r.q_tag.length(); i++){
            *p++ = r.q_tag[r.q_tag.size() - i - 1] - 33;
        }
    }


    if((p - bam->data) != dl){
        std::cout << "Data length wrong dl = " << dl << " compared to " << (p - bam->data) << "\n";
        exit(1);
    }

    //uint32_t ol = bam->m_data;
    uint32_t tmpu = a.NM;
    bam_aux_append(bam, "NM", 'i', 4, reinterpret_cast<uint8_t*>(&tmpu));
    char     tmpc = a.xs;
    bam_aux_append(bam, "XS", 'A', 1, reinterpret_cast<uint8_t*>(const_cast<char*>(&tmpc)));
    bam_aux_append(bam, "UB", 'Z', r.umi.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(r.umi.c_str())));
    bam_aux_append(bam, "CB", 'Z', r.barcode.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(r.barcode.c_str())));
    tmpu = a.score;
    bam_aux_append(bam, "AS", 'i', 4, reinterpret_cast<uint8_t*>(&tmpu));
    if(it->gid != std::numeric_limits<uint32_t>::max()){
        const std::string & gid = idx.gene(it->gid).gene_id;
        bam_aux_append(bam, "XG", 'Z', gid.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(gid.c_str())));
        tmpu = it->gid;
        bam_aux_append(bam, "XI", 'i', 4, reinterpret_cast<uint8_t*>(&tmpu));
    }

    tmpc = aligntype2code(a.atype);
    //std::cout << "  RE = " << tmpc << "\n";
    bam_aux_append(bam, "RE", 'A', 1, reinterpret_cast<uint8_t*>(&tmpc));

    if(g.res == AlignGroup::MULTIMAPPED){
        std::stringstream ss;
        bool first = true;
        uint32_t lgid = it->gid;
        /*
        uint32_t llgid = it->gid;
        for(auto aa = g.alns.begin(); aa != aend; aa++){
            if(aa->passed && (((aa->atype == CDNA || aa->atype == INTRON) && aa->gid != llgid) || aa->atype == INTERGENIC)){
                std::cout << "  MM: " << aa->print() << "\n";
                lgid = aa->gid;
            }
        }
        */
        while(++it != aend){
            if(it->passed && (((it->atype == CDNA || it->atype == INTRON) && it->gid != lgid) || it->atype == INTERGENIC)){
                if(!first) ss << ';';
                ss << bh.tname(it->tid) << ',' << it->xs << it->lft << ',' << it->cigar << ',';
                if(it->gid != std::numeric_limits<uint32_t>::max()){
                    const std::string & gid = idx.gene(it->gid).gene_id;
                    ss << gid;
                }else{
                    ss << '-';
                }
                ss << ',' << aligntype2code(a.atype);
                first = false;
                lgid = it->gid;
            }
        }
        std::string so = ss.str();
        bam_aux_append(bam, "XA", 'Z', so.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(so.c_str())));
        //std::cout << "  XA: " << so << "\n";
    }
    tmpc = AlignGroup::alignres2code(g.res);
    bam_aux_append(bam, "RA", 'A', 1, reinterpret_cast<uint8_t*>(&tmpc));
}

