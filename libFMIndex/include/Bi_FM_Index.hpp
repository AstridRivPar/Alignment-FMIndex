#pragma once

#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <sdsl/csa_wt.hpp>
#include <ranges>

#include <Poco/Logger.h>

#include "Helpers.hpp"

// using backward_index_type = sdsl::csa_wt<sdsl::wt_int<sdsl::rrr_vector<63>>,
//                                         16,
//                                         10'000'000,
//                                         sdsl::sa_order_sa_sampling<>,
//                                         sdsl::isa_sampling<>,
//                                         sdsl::int_alphabet<>>;

// using forward_index_type = sdsl::csa_wt<sdsl::wt_int<sdsl::rrr_vector<63>>,
//                                         10'000'000,
//                                         10'000'000,
//                                         sdsl::sa_order_sa_sampling<>,
//                                         sdsl::isa_sampling<>,
//                                         sdsl::int_alphabet<>>;


using backward_index_type = sdsl::csa_wt<sdsl::wt_int<sdsl::bit_vector>,
                                        16,            // Sampling rate of the suffix array
                                        10'000'000, // Sampling rate of the inverse suffix array
                                        sdsl::sa_order_sa_sampling<>,
                                        sdsl::isa_sampling<>,
                                        sdsl::int_alphabet<>>;

using forward_index_type = sdsl::csa_wt<sdsl::wt_int<sdsl::bit_vector>,
                                        10'000'000,
                                        10'000'000,
                                        sdsl::sa_order_sa_sampling<>,
                                        sdsl::isa_sampling<>,
                                        sdsl::int_alphabet<>>;

class Bi_FM_Index{
    static Poco::Logger &log;
    typedef uint8_t char_type;

    private:
        backward_index_type backward_index;
        forward_index_type forward_index;
        size_t FM_size;
        
    public:

        Bi_FM_Index() = default;
        // Bi_FM_Index(const std::vector<uint64_t> &vec){
        // sdsl::int_vector<> v(vec.begin(), vec.end()); //keep
        // sdsl::int_vector<> v_rev(vec.rbegin(), vec.rend());
        Bi_FM_Index(const std::string &forward_string){
            
            sdsl::construct_im(backward_index, forward_string,1);
            log.information("Backward Idx build successfully!", __FILE__,__LINE__);
            assert(backward_index.wavelet_tree.sigma < 256 && 
                "The size of the alphabet must be less than 255!");

            std::string backward_string = std::string(forward_string.rbegin(), forward_string.rend());
            sdsl::construct_im(forward_index, backward_string, 1);
            log.information("Forward Idx build successfully!", __FILE__,__LINE__);
            log.information("Index size: " + std::to_string(backward_index.size()), __FILE__,__LINE__);

            assert(backward_index.size()==forward_index.size()&&
                    "Error while building indexes!");

            FM_size = backward_index.size();

        }
        ~Bi_FM_Index() = default; //!< Defaulted.

        size_t size(){
            return FM_size;
        }
        
        bool backward_search(SARangePair &sa_pair, uint8_t next_c);
            

        bool forward_search(SARangePair &sa_pair, uint8_t next_c);
        // bool backward_search(SARangePair &sa_pair, uint64_t next_c);

        // bool forward_search(SARangePair &sa_pair, uint64_t next_c);
        // bool extend_char(SARangePair &sa_pair, uint64_t next_c);

        friend std::ostream& operator<<(std::ostream& os, const Bi_FM_Index &idx);

        const uint8_t get_alphabet_size(){
            return backward_index.wavelet_tree.sigma;
        }

        const uint8_t char2comp(char c){
            
            return backward_index.char2comp[c];

        }
        const char comp2char(uint8_t n){

            return backward_index.comp2char[n];
        }

        size_t count(const SARangePair & range) const noexcept
        {
            assert(FM_size > 0 && range.isCorrect());

            return range.width();
        }
        size_t offset(const SARangePair &range) const noexcept
        {
            assert(FM_size > range.width());
            return FM_size - range.width() - 1; // since the string is reversed during construction
        }
        std::vector<size_t> get_occs(const SARangePair &range){
            assert(FM_size > 0);
            std::vector<size_t> occ{};
            occ.reserve(range.width());
            // for (size_t i = 0; i < count(range); ++i)
            // {
            //     occ.emplace_back(0, offset(range) - backward_index[range.getRangeSA().getBegin() + i]);
            // }
            return occ;
        }

        size_t get_occ(const SARangePair &range){
            assert(FM_size > 0);
            return backward_index[range.getRangeSA().getBegin()];
        
        }

};

