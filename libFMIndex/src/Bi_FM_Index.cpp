#include "Bi_FM_Index.hpp"
Poco::Logger &Bi_FM_Index::log = Poco::Logger::get("Bi_FM_Index");

// bool Bi_FM_Index::backward_search(SARangePair &sa_pair, uint64_t c)
bool Bi_FM_Index::backward_search(SARangePair &sa_pair, uint8_t c){
    Range &bw_range = sa_pair.rangeSA;
    Range &fw_range = sa_pair.rangeSARev;
    sa_pair.length += 1;
    auto cc = backward_index.char2comp[c];

    auto const c_begin = backward_index.C[cc];
    auto const c_end = backward_index.C[cc + 1];

    if (sa_pair.width() == size()) {
        bw_range.begin = fw_range.begin = c_begin;
        bw_range.end = fw_range.end = c_end;
        // assert(bw_range == fw_range);
        return !sa_pair.empty();
    }

    if( c_begin >= c_end ){
        bw_range.begin = fw_range.begin = bw_range.end = fw_range.end = c_begin;
        return false;
    }    
    // backward_index.wavelet_tree.lex_count(bw_range.begin, )
    const auto [rank_l, smaller, greater] = backward_index.wavelet_tree.lex_count(bw_range.begin, bw_range.end, c);
    const auto rank_r = bw_range.end - bw_range.begin - smaller - greater + rank_l;
   
    fw_range.begin += smaller;
    fw_range.end -= greater;

    bw_range.begin = c_begin + rank_l;
    bw_range.end = c_begin + rank_r;
    
    return (!sa_pair.empty());

}

bool Bi_FM_Index::forward_search(SARangePair &sa_pair, uint8_t c){
    
    Range &bw_range = sa_pair.rangeSA;
    Range &fw_range = sa_pair.rangeSARev;
    sa_pair.length += 1;
    auto cc = backward_index.char2comp[c];

    auto const c_begin = backward_index.C[cc];
    auto const c_end = backward_index.C[cc + 1];

    if (sa_pair.width() == size()) {
        bw_range.begin = fw_range.begin = c_begin;
        bw_range.end = fw_range.end = c_end;
        // assert(bw_range == fw_range);
        return !fw_range.empty();
    }
    if( c_begin >= c_end ){
        // LOG(WARNING) << "Range not found";
        bw_range.begin = fw_range.begin = bw_range.end = fw_range.end = c_begin;
        return false;
    }    
    const auto [rank_l, smaller, greater] = forward_index.wavelet_tree.lex_count(fw_range.begin, fw_range.end, c);
    const auto rank_r = fw_range.end - fw_range.begin - smaller - greater + rank_l;

    bw_range.begin += smaller;
    bw_range.end -= greater;

    fw_range.begin = c_begin + rank_l;
    fw_range.end = c_begin + rank_r;
    
    return (!sa_pair.empty());
}

std::ostream& operator<<(std::ostream&os, const Bi_FM_Index &idx){
        os <<" i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << std::endl;
        sdsl::csXprintf(os, "%2I %2S %3s %3P %2p %3B   %:3T", idx.backward_index);
        std::cout << std::endl;
        os <<" i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << std::endl;
        sdsl::csXprintf(os, "%2I %2S %3s %3P %2p %3B   %:3T", idx.forward_index);
        std::cout << std::endl;

    return os;
}