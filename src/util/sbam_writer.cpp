#include "sbam_writer.hpp"
using namespace gwsc;
#include <iomanip>
#include <sstream>


void SortedBamWriter::write(const AlignGroup & g, const Read & r, read_buffer & buffer, unsigned int & rcount, const std::string & gid){
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

    _align2bam(s, g, r, gid);
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
    total_reads_ += file_reads_;
    file_reads_ = 0;
    file_number_++;
}

void SortedBamWriter::_align2bam(bam1_t * bam, const AlignGroup & g, const Read & r, const std::string & gid){
    auto it = g.alns.begin();
    while(it != g.alns.end() && it->gid != g.summary.gene_id){
        it++;
    }
    if(it == g.alns.end()){
        std::cout << "Could not find alignment for matching gene id\n";
        exit(1);
    }
    auto & a = *it;
    BamFlag flag;
    //flag.f.read2 = true;
    flag.f.strand = a.rev;

    bam->core.flag = flag.flag_val;
    bam->core.qual = 255; //a.mapq;  All "uniquely" mapped reads
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
    char     tmpc = a.xs;
    bam_aux_append(bam, "NM", 'i', 4, reinterpret_cast<uint8_t*>(&tmpu));
    bam_aux_append(bam, "XS", 'A', 1, reinterpret_cast<uint8_t*>(const_cast<char*>(&tmpc)));
    bam_aux_append(bam, "UB", 'Z', r.umi.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(r.umi.c_str())));
    bam_aux_append(bam, "CB", 'Z', r.barcode.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(r.barcode.c_str())));
    bam_aux_append(bam, "XG", 'Z', gid.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(gid.c_str())));
    tmpu = a.score;
    bam_aux_append(bam, "AS", 'i', 4, reinterpret_cast<uint8_t*>(&tmpu));
    tmpu = a.gid;
    bam_aux_append(bam, "XI", 'i', 4, reinterpret_cast<uint8_t*>(&tmpu));

    tmpc = '?';
    if(a.atype == ANTISENSE){
        tmpc = 'A';
    }else if(a.atype == AlignType::CDNA){
        tmpc = 'E';
    }else if(a.atype == AlignType::INTRON){
        tmpc = 'N';
    }else if(a.atype == AlignType::INTERGENIC){
        tmpc = 'I';
    }
    bam_aux_append(bam, "RE", 'A', 1, reinterpret_cast<uint8_t*>(&tmpc));
}

