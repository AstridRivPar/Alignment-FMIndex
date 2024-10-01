/* bidirectional search algorithms library
    Copyright (C) 2010 Thomas Schnattinger

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file index_bidirectional_waveletindex.hpp
    \brief Bidirectional Wavelet Index
	\author Thomas Schnattinger
*/
#ifndef INCLUDED_BIDIRECTIONAL_INDEX_BIDIRECTIONAL_WAVELETINDEX
#define INCLUDED_BIDIRECTIONAL_INDEX_BIDIRECTIONAL_WAVELETINDEX

#include <boost/test/utils/nullstream.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/algorithms.hpp>
#include "indexes.hpp"
#include "wavelet_tree.hpp"


using namespace sdsl;

/*! This is an implementation of the bidirectional wavelet index.
 * It is parameterized by a wavelet tree (supporting extract, 
 * occ and getBounds operations) and the number SampleDens specifying
 * how many samples of SA are explicitly stored
 */
template<class WaveletTree = wavelet_tree<unsigned char, bit_vector, rank_support_v<> >, uint32_t SampleDens = 5>
class index_bidirectional_waveletindex : public index {

	private:

		WaveletTree backward_index;
		WaveletTree forward_index;
		int_vector<> m_sa_sample;

	public:
	
		struct searchstate {
			size_t occ_begin;
			size_t occ_end;
		   	size_t occ_forward_begin;
			size_t occ_forward_end;
			size_t length;
			// std::string pMatch = "";
		
			searchstate() : occ_begin(0), occ_end(0), occ_forward_begin(0), occ_forward_end(0), length(0){}
			
			searchstate(size_t i1, size_t j1, size_t i2, size_t j2, size_t l) : occ_begin(i1), occ_end(i2), occ_forward_begin(j1), occ_forward_end(j2), length(l){}
			searchstate(const searchstate &s) : occ_begin(s.occ_begin), occ_end(s.occ_end), occ_forward_begin(s.occ_forward_begin), occ_forward_end(s.occ_forward_end), length(s.length){}
			
			// searchstate(size_t i1, size_t j1, size_t i2, size_t j2, size_t l, std::string pMatch) : occ_begin(i1), occ_end(i2), occ_forward_begin(j1), occ_forward_end(j2), length(l), pMatch(pMatch){}
			// searchstate(const searchstate &s) : occ_begin(s.occ_begin), occ_end(s.occ_end), occ_forward_begin(s.occ_forward_begin), occ_forward_end(s.occ_forward_end), length(s.length), pMatch(s.pMatch){}
			
					
		};

		//! Default constructor
		index_bidirectional_waveletindex() : index() {}

		/*! 
		 * Constructor for building the Index
		 * \param[in] str C-string of the text
		 */
		index_bidirectional_waveletindex(const unsigned char* str) : index() {

			size_t n = strlen((const char*)str);
			// std::cout<< n << std::endl;
			int_vector<> sa(n+1, 0, bit_magic::l1BP(n+1)+1);
			setText(str, n+1);

			unsigned char *bwt = new unsigned char[n+1];
			
			algorithm::calculate_sa(str, n+1, sa);   // calculate the suffix array sa of str

			{ /* Calculate Burrows-Wheeler-Transform */
				size_t i = 0;
				for(int_vector<>::const_iterator it = sa.begin(), end = sa.end(); it != end; ++it, ++i){
					bwt[i] = m_char2comp[str[(*it+n)%(n+1)]];
				}
			}
			
			backward_index = WaveletTree(bwt, n+1, m_sigma);
			/* Construct the SA-Samples */
			m_sa_sample.set_int_width( bit_magic::l1BP(sa.size())+1 );
			m_sa_sample.resize( (sa.size()+SampleDens-1)/SampleDens );
			size_t idx=0;
			size_t i=(sa.size()-1-SampleDens*(m_sa_sample.size()-1));
			for(int_vector<>::const_iterator it = sa.begin()+(ptrdiff_t)i; i < sa.size(); it += (ptrdiff_t)SampleDens, i += SampleDens, ++idx){
				m_sa_sample[idx] = *it;
			} 

			unsigned char* reverse = new unsigned char[n+1];
			for (size_t i=0; i<n; i++) reverse[i] = str[n-1-i];
			reverse[n] = '\0';
			//std::cout<< reverse <<std::endl;

			algorithm::calculate_sa(reverse, n+1, sa);   // calculate the suffix array sa of reverse string str

			{ /* Calculate Burrows-Wheeler-Transform */
				size_t i = 0;
				for(int_vector<>::const_iterator it = sa.begin(), end = sa.end(); it != end; ++it, ++i){
					bwt[i] = m_char2comp[reverse[(*it+n)%(n+1)]];
				}
			}

			forward_index = WaveletTree(bwt, n+1, m_sigma);
			
			delete [] bwt;
			delete [] reverse;

		}

		~index_bidirectional_waveletindex(){

		}
		/*!
		 * Initializes a searchstate to the searchstate of the empty word
		 * \param[out] s searchstate to be initialized
		 */
		void init_search_state(searchstate &s);

		/*!
		 * Search for character c in forward direction, starting with a given searchstate
		 * \param[in] c character to search for
		 * \param[in,out] s searchstate (will be overwritten)
		 * \return size of the new interval
		 */
		size_t forward_search(unsigned char c, searchstate &s);

		searchstate fs(unsigned char c, const searchstate &original_s);
		searchstate bs(unsigned char c, const searchstate &original_s);
		
		/*!
		 * Search for character c in backward direction, starting with a given searchstate
		 * \param[in] c character to search for
		 * \param[in,out] s searchstate (will be overwritten)
		 * \return size of the new interval
		 */
		size_t backward_search(unsigned char c, searchstate &s);

		/*!
		 * Extract the SA-values of the interval $[i..j]$
		 * \param[in] i begin of the interval (including)
		 * \param[in] j end of the interval (excluding)
		 * \param[out] occ result
		 */
		void extract_sa(size_t i, size_t j, int_vector<> &occ);

		/*!
		 * Serializes this data structure to an output stream
		 * \param[in,out] out output stream to write to
		 * \return number of bytes written to out
		 */
		size_t serialize(std::ostream &out) const;

		/*!
		 * Load this data structure from an input stream.
		 * \param[in,out] in input stream to read from
		 */
		void load(std::istream &in);

		/*!
		 * Calculate the used memory of this data structure by serializing into /dev/null
		 * \return size of this data structure in bytes
		 */
		size_t used_bytes() const;

		/*!
		 * Extract the SA-values of the interval $[i..j]$
		 * \param[in] i begin of the interval (including)
		 * \param[in] j end of the interval (excluding)
		 */
		size_t extract(size_t i, size_t j);

		void extract_match(const searchstate &s) ;

		

};	

template<class WaveletTree, uint32_t SampleDens>
void index_bidirectional_waveletindex<WaveletTree, SampleDens>::init_search_state(searchstate& s) {
	s.occ_begin = 0;
	s.occ_end = m_size;
	s.occ_forward_begin = 0;
	s.occ_forward_end = m_size;
	s.length = 0;

}

template<class WaveletTree, uint32_t SampleDens>
size_t index_bidirectional_waveletindex<WaveletTree, SampleDens>::forward_search(unsigned char c, searchstate &s) {
//	c = sa.char2comp[c];

	if (s.occ_begin==0 && s.occ_end==m_size) {
		s.occ_begin = s.occ_forward_begin = m_C[c];
		s.occ_end = s.occ_forward_end = m_C[c+1];
		s.length = 1;
		return s.occ_end - s.occ_begin;
	}

	size_t lesser, greater;
	
	forward_index.getBounds(s.occ_forward_begin, s.occ_forward_end, c, lesser, greater);
	s.occ_begin += lesser;
	s.occ_end -= greater;

	size_t c_begin = m_C[c];  // interval begin of character c (inclusive)
	size_t c_end   = m_C[c+1];// interval end of character c (exclusive)

	if( c_begin >= c_end ){
		s.occ_begin = s.occ_end = s.occ_forward_begin = s.occ_forward_end = c_end = c_begin;
		return 0;
	}    

	size_t c_before = forward_index.occ(c, s.occ_forward_begin);
	size_t c_after = m_C[c+1] - m_C[c] - forward_index.occ(c, s.occ_forward_end);

	s.occ_forward_begin = c_begin + c_before;
	s.occ_forward_end = c_end - c_after;

	s.length += 1;
	return s.occ_forward_end-s.occ_forward_begin;
}
template<class WaveletTree, uint32_t SampleDens>
typename index_bidirectional_waveletindex<WaveletTree, SampleDens>::searchstate index_bidirectional_waveletindex<WaveletTree, SampleDens>::fs(unsigned char c, const searchstate &original_s) {
//	c = sa.char2comp[c];
	searchstate s;
	if (original_s.occ_begin==0 && original_s.occ_end==m_size) {
		s.occ_begin = s.occ_forward_begin = m_C[c];
		s.occ_end = s.occ_forward_end = m_C[c+1];
		s.length = 1;
		return s;
	}

	size_t lesser, greater;
	forward_index.getBounds(original_s.occ_forward_begin, original_s.occ_forward_end, c, lesser, greater);
	//s.occ_begin += lesser;
	//s.occ_end -= greater;

	s.occ_begin = original_s.occ_begin + lesser;
	s.occ_end = original_s.occ_end - greater;

	size_t c_begin = m_C[c];  // interval begin of character c (inclusive)
	size_t c_end   = m_C[c+1];// interval end of character c (exclusive)

	if( c_begin >= c_end ){
		s.occ_begin = s.occ_end = s.occ_forward_begin = s.occ_forward_end = c_end = c_begin;
		return s;
	}    

	size_t c_before = forward_index.occ(c, original_s.occ_forward_begin);
	size_t c_after = m_C[c+1] - m_C[c] - forward_index.occ(c, original_s.occ_forward_end);

	s.occ_forward_begin = c_begin + c_before;
	s.occ_forward_end = c_end - c_after;

	s.length = original_s.length + 1;
	//return s.occ_forward_end-s.occ_forward_begin;
	// s.pMatch = original_s.pMatch;
	return s;
}
template<class WaveletTree, uint32_t SampleDens>
typename index_bidirectional_waveletindex<WaveletTree, SampleDens>::searchstate index_bidirectional_waveletindex<WaveletTree, SampleDens>::bs(unsigned char c, const searchstate &original_s) {  
	
	searchstate s;
	if (original_s.occ_begin==0 && original_s.occ_end==m_size) {
		s.occ_begin = s.occ_forward_begin = m_C[c];
		s.occ_end = s.occ_forward_end = m_C[c+1];
		s.length = 1;
		return s;
	}

	size_t lesser, greater;
	backward_index.getBounds(original_s.occ_begin, original_s.occ_end, c, lesser, greater);

	s.occ_forward_begin = original_s.occ_forward_begin + lesser;
	s.occ_forward_end = original_s.occ_forward_end - greater;

	size_t c_begin = m_C[c];  // interval begin of character c (inclusive)
	size_t c_end   = m_C[c+1];// interval end of character c (exclusive)

	if( c_begin >= c_end ){
		s.occ_begin = s.occ_end = s.occ_forward_begin = s.occ_forward_end = c_end = c_begin;
		return s;
	}

	size_t c_before = backward_index.occ(c, original_s.occ_begin);
	size_t c_after = m_C[c+1] - m_C[c] - backward_index.occ(c, original_s.occ_end);

	s.occ_begin = c_begin + c_before;
	s.occ_end = c_end - c_after;

	s.length = original_s.length + 1;
	//return s.occ_forward_end-s.occ_forward_begin;
	// s.pMatch = original_s.pMatch;

	return s;
}
template<class WaveletTree, uint32_t SampleDens>
size_t index_bidirectional_waveletindex<WaveletTree, SampleDens>::backward_search(unsigned char c, searchstate &s) {
//	c = sa.char2comp[c];

	if (s.occ_begin==0 && s.occ_end==m_size) {
		s.occ_begin = s.occ_forward_begin = m_C[c];
		s.occ_end = s.occ_forward_end = m_C[c+1];
		s.length = 1;
		return s.occ_end - s.occ_begin;
	}

	size_t lesser, greater;
	backward_index.getBounds(s.occ_begin, s.occ_end, c, lesser, greater);
	s.occ_forward_begin += lesser;
	s.occ_forward_end -= greater;

	size_t c_begin = m_C[c];  // interval begin of character c (inclusive)
	size_t c_end   = m_C[c+1];// interval end of character c (exclusive)

	if( c_begin >= c_end ){
		s.occ_begin = s.occ_end = s.occ_forward_begin = s.occ_forward_end = c_end = c_begin;
		return 0;
	}    

	size_t c_before = backward_index.occ(c, s.occ_begin);
	size_t c_after = m_C[c+1] - m_C[c] - backward_index.occ(c, s.occ_end);

	s.occ_begin = c_begin + c_before;
	s.occ_end = c_end - c_after;

	s.length += 1;
	return s.occ_end-s.occ_begin;
}


template<class WaveletTree, uint32_t SampleDens>
void index_bidirectional_waveletindex<WaveletTree, SampleDens>::extract_sa(size_t occ_begin, size_t occ_end, int_vector<> &occ) {
	occ.resize(occ_end-occ_begin);
	for (size_t k=0; k<occ_end-occ_begin; k++) {
		size_t i = occ_begin+k;
		size_t off = 0;
		while( (size()-1-i) % SampleDens != 0 ){// while SA[i] is not sampled
			// go to the position where SA[i]-1 is located
			// i := LF(i)
			uint64_t bwt_result = backward_index.extract(i);
			i = m_C[bwt_result] + backward_index.occ(bwt_result, i); 
			// add 1 to the offset
			++off;
		}   
		occ[k] = (m_sa_sample[i / SampleDens]+off) % size();
	}
}

template<class WaveletTree, uint32_t SampleDens>
size_t index_bidirectional_waveletindex<WaveletTree, SampleDens>::extract(size_t occ_begin, size_t occ_end) {
	size_t i = occ_begin;
	size_t off = 0;
	while( (size()-1-i) % SampleDens != 0 ){// while SA[i] is not sampled
		// go to the position where SA[i]-1 is located
		// i := LF(i)
		uint64_t bwt_result = backward_index.extract(i);
		i = m_C[bwt_result] + backward_index.occ(bwt_result, i); 
		// add 1 to the offset
		++off;
	}
	return (m_sa_sample[i / SampleDens]+off) % size();
}

template<class WaveletTree, uint32_t SampleDens>
void index_bidirectional_waveletindex<WaveletTree, SampleDens>::extract_match(const searchstate &s) {
	
	uint64_t bwt_result = backward_index.extract(s.occ_begin);
	size_t i = m_C[bwt_result] + forward_index.occ(bwt_result, s.occ_forward_begin) - 1;
	std::cout << comp2char[bwt_result];

	for(size_t counter = 0; counter < s.length ; counter++){
		uint64_t bwt_result = forward_index.extract(i);
		i = m_C[bwt_result] + forward_index.occ(bwt_result, i); 
		std::cout << comp2char[bwt_result];
	}

	std::cout << std::endl;
}

template<class WaveletTree, uint32_t SampleDens>
size_t index_bidirectional_waveletindex<WaveletTree, SampleDens>::serialize(std::ostream &out) const {
	size_t written_bytes = 0;
	out.put('6');
	written_bytes += sizeof(char);
	written_bytes += superserialize(out);
	written_bytes += backward_index.serialize(out);
	written_bytes += forward_index.serialize(out);
	written_bytes += m_sa_sample.serialize(out);
	return written_bytes;
}

template<class WaveletTree, uint32_t SampleDens>
void index_bidirectional_waveletindex<WaveletTree, SampleDens>::load(std::istream &in) {
	if (in.get()!='6') {
		std::cerr << "wrong index!!";
		throw("wrong index!!");
	}
	superload(in);
	backward_index.load(in);
	forward_index.load(in);
	m_sa_sample.load(in);
}

template<class WaveletTree, uint32_t SampleDens>
size_t index_bidirectional_waveletindex<WaveletTree, SampleDens>::used_bytes() const {
	boost::onullstream dev_null;
   	return serialize(dev_null);
}

#endif
