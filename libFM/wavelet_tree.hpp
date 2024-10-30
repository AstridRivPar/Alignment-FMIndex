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
/*! \file wavelet_tree.hpp
    \brief A Wavelet Tree implementation supporting occ, extract 
	and getBounds operations
	\author Thomas Schnattinger
*/
#ifndef INCLUDED_SDSL_WAVELET_TREE
#define INCLUDED_SDSL_WAVELET_TREE

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
// #include <sdsl/algorithms.hpp>

using namespace sdsl;

/*! A class for the Wavelet Tree. The Wavelet Tree is parameterized with 
 * an RankSupport and the data type of the supported vector
 */
template <typename DataType = unsigned char, class BitVector = bit_vector, class RankSupport = rank_support_v<> >
class wavelet_tree {
	private:

		BitVector *B;
		RankSupport *rank_B;
		// m_size = 2^(ceil(log2(m_sigma)))
		size_t m_size;
		uint8_t m_sigma;

	public:
		//! Default Constructor
		wavelet_tree() : m_sigma(0), m_size(0), B(NULL), rank_B(NULL) {}
		//! Standard Constructor
		/*!
		 *  \param bwt the sequence this tree is built for
		 *  \param len the length of bwt
		 *  \param sigma the size of the alphabet
		 */
		wavelet_tree(DataType bwt[], size_t len, uint8_t sigma) : m_sigma(sigma) {
			size_t logsigma = 1;
			// m_size = 2^(ceil(log2(m_sigma)))
			m_size = 1;
			while (m_size<m_sigma) {
				m_size <<= 1;
				logsigma++;
			}

			B = new BitVector[m_size];
			rank_B = new RankSupport[m_size];

			DataType *temp = new DataType[2*len];
			DataType* const t[2] = {temp, temp+len};

			// fill the first array initially
			for (size_t i=0; i<len; i++) t[0][i] = bwt[i];

			// create stack
			size_t treepos[logsigma], start[logsigma], end[logsigma];
			uint8_t which_t[logsigma], sigma_begin[logsigma], sigma_end[logsigma];
			// push root onto the stack
			treepos[0] = 1;
			start[0] = 0; // left border of the interval (inclusive)
			end[0] = len; // right border of the interval (exclusive)
			which_t[0] = 0;
			sigma_begin[0] = 0; // first character of alphabet (inclusive)
			sigma_end[0] = m_sigma; // last character of alphabet (exclusive)
			size_t pointer = 1;
			// while stack not empty
			while (pointer>0) {
				--pointer;
				size_t s_treepos=treepos[pointer], s_start=start[pointer], s_end=end[pointer];
				uint8_t s_which_t=which_t[pointer],
						s_sigma_begin=sigma_begin[pointer], 
						s_sigma_end=sigma_end[pointer];

				// continue if node consists of only one character
			   	if (s_sigma_end - s_sigma_begin <= 1) 
					continue;

				uint8_t sigma_mid = (s_sigma_begin+s_sigma_end) / 2;

				DataType *from = t[s_which_t];
				DataType *to = t[1-s_which_t];

				DataType *left = to+s_start;
				DataType *right = to+s_end;

				// resize bit vector
				bit_vector b;
				b.resize(s_end-s_start);
				size_t size_of_left_child = 0;
				if (s_treepos==1 || s_treepos%2==0) { // left-to-right
					for (size_t k=0, i=s_start; i<s_end; i++, k++) {
						//assert(from[i]>=s_sigma_begin && from[i]<s_sigma_end);
						if (from[i] < sigma_mid) {
							// -> left child
							*(left++) = from[i];
							b[k] = 0;
							size_of_left_child++;
						} else {
							// -> right child
							*(--right) = from[i];
							b[k] = 1;
						}
					}
				} else { // right-to-left
					for (size_t k=0, i=s_end-1; i!=(size_t)-1 && i>=s_start; i--, k++) {
						//assert(from[i]>=s_sigma_begin && from[i]<s_sigma_end);
						if (from[i] < sigma_mid) {
							// -> left child
							*(left++) = from[i];
							b[k] = 0;
							size_of_left_child++;
						} else {
							// -> right child
							*(--right) = from[i];
							b[k] = 1;
						}
					}
				}
				assert(left==right);

				B[s_treepos] = BitVector(b);

				// init rank support
				rank_B[s_treepos].init(&B[s_treepos]);

				// push right child
				treepos[pointer] = s_treepos*2+1;
				start[pointer] = s_start + size_of_left_child;
				end[pointer] = s_end;
				which_t[pointer] = 1-s_which_t;
				sigma_begin[pointer] = sigma_mid;
				sigma_end[pointer] = s_sigma_end;
				pointer++;

				// push left child
				treepos[pointer] = s_treepos*2;
				start[pointer] = s_start;
				end[pointer] = s_start + size_of_left_child;
				which_t[pointer] = 1-s_which_t;
				sigma_begin[pointer] = s_sigma_begin;
				sigma_end[pointer] = sigma_mid;
				pointer++;

			}
			delete[] temp;
			temp = NULL;
		}

		//! Default Destructor
		~wavelet_tree() {
			if (B!=NULL) {
				delete[] B;
				B = NULL;
			}
			if (rank_B!=NULL) {
				delete[] rank_B;
				rank_B = NULL;
			}
		}

		//! Copy constructor
		wavelet_tree(wavelet_tree<RankSupport, DataType> &wt) {
			m_sigma = wt.m_sigma;
			m_size = wt.m_size;
			B = new BitVector[m_size];
			rank_B = new RankSupport[m_size];
			for (size_t i=0; i<m_size; i++) {
				B[i] = BitVector(wt.B[i]);
				rank_B = RankSupport(wt.rank_B[i]);
			}
		}
			   
		//! Swap method for wavelet_tree
  		/*! \param wt wavelet_tree to swap.
		 *  Required for the Assignable Conecpt of the STL.
		 */
		void swap(wavelet_tree &wt);

		//! Serialize to a stream.
		/*! \param out Outstream to write the data structure.
		 *  \return The number of written bytes.
		 */
		size_t serialize(std::ostream &out) const;

		//! Load from a stream.
		/*! \param in Inputstream to load the data structure from.
		 */
		void load(std::istream &in);

		//! Assignment Operator.
		/*!
		 *  Required for the Assignable Concept of the STL.
		 */
		wavelet_tree& operator=(const wavelet_tree &wt);

		//! Equality Operator
		/*! Two Instances of wavelet_tree are equal if all bit vectors are equal.
		 *  \par Required for the Equality Comparable Concept of the STL.
		 *  \sa operator!=
		 */
		bool operator==(const wavelet_tree &wt)const;

		//! Unequality Operator
		/*! Two Instances of wavelet_tree are equal if all bit vectors are equal.
		 *  \par Required for the Equality Comparable Concept of the STL.
		 *  \sa operator==
		 */
		bool operator!=(const wavelet_tree &wt)const;

	  	//! Number of elements in the wavelet tree.
		size_t size() const { return m_size;}
	
		//! Character-Rank Query
		/*! Returns the number of occurrences of character c up to position pos 
		 *  in the underlying sequence
		 *  \sa operator()
		 */
		size_t occ(DataType c, size_t pos) const;

		//! Same as occ(DataType c, size_t pos)
		size_t operator() (DataType c, size_t pos) const {
			return occ(c, pos);
		}

		//! Extract character from position pos
		/*! Returns the character at position pos of the underlying sequence
		 *  \sa operator[]
		 */
		DataType extract(size_t pos) const;

		//! Same as extract(size_t pos)
		DataType operator[] (size_t pos) const {
			return extract(pos);
		}

		//! Wie in der Ausarbeitung...
		/*! Returns the number of characters in bwt[i..j] lesser than c and greater than c
	     */
		void getBounds(size_t i, size_t j, DataType c, size_t &lesser, size_t &greater);
};

template <typename DataType, class BitVector, class RankSupport>
void wavelet_tree<DataType, BitVector, RankSupport>::swap(wavelet_tree &wt) {
	std::swap(m_sigma, wt.m_sigma);
	std::swap(m_size, wt.m_size);
	std::swap(B, wt.B);
	std::swap(rank_B, wt.rank_B);
}

template <typename DataType, class BitVector, class RankSupport>
size_t wavelet_tree<DataType, BitVector, RankSupport>::serialize(std::ostream &out) const {
	size_t written_bytes = 0;
	out.write((char*)&m_sigma, sizeof(m_sigma));
	written_bytes += sizeof(m_sigma);
	out.write((char*)&m_size, sizeof(m_size));
	written_bytes += sizeof(m_size);
	for (size_t i=0; i<m_size; i++) {
		written_bytes += B[i].serialize(out);
		written_bytes += rank_B[i].serialize(out);
	}
	return written_bytes;
}

template <typename DataType, class BitVector, class RankSupport>
void wavelet_tree<DataType, BitVector, RankSupport>::load(std::istream &in){
	in.read((char*)&m_sigma, sizeof(m_sigma));
	in.read((char*)&m_size, sizeof(m_size));
   	B = new BitVector[m_size];
	rank_B = new RankSupport[m_size];
	for (size_t i=0; i<m_size; i++) {
		B[i].load(in);
		rank_B[i].load(in, B+i);
	}
}

template <typename DataType, class BitVector, class RankSupport>
wavelet_tree<DataType, BitVector, RankSupport>& wavelet_tree<DataType, BitVector, RankSupport>::operator=(const wavelet_tree<DataType, BitVector, RankSupport> &wt) {
	if(this != &wt){
		m_size = wt.m_size;
		m_sigma = wt.m_sigma;
		if (B!=NULL){
			delete[] B;
			B = NULL;
		}
		if (rank_B!=NULL){
			delete[] rank_B;
			rank_B=NULL;
		}
		B = new BitVector[m_size];
		rank_B = new RankSupport[m_size];
		for (size_t i=0; i<m_size; i++) {
			B[i] = wt.B[i];
			rank_B[i] = wt.rank_B[i];
			rank_B[i].init(B+i);
		}
	}
	return *this;
}

template <typename DataType, class BitVector, class RankSupport>
bool wavelet_tree<DataType, BitVector, RankSupport>::operator==(const wavelet_tree<DataType, BitVector, RankSupport> &wt) const {
	if (this != &wt) {
		if (m_sigma!=wt.m_sigma) return false;
		for (size_t i=0; i<m_size; i++) if (B[i]!=wt.B[i]) return false;
	} else return true;
}

template <typename DataType, class BitVector, class RankSupport>
bool wavelet_tree<DataType, BitVector, RankSupport>::operator!=(const wavelet_tree<DataType, BitVector, RankSupport> &wt) const {
	return (*this != wt);
}

template <typename DataType, class BitVector, class RankSupport>
size_t wavelet_tree<DataType, BitVector, RankSupport>::occ(DataType c, size_t pos) const {
	size_t pointer = 1;
	uint8_t sigma_start = 0;
	uint8_t sigma_end = m_sigma;
	while (sigma_end-sigma_start>1) {
		if (sigma_end-sigma_start==1) break;
		uint8_t sigma_mid = (sigma_start+sigma_end)/2;
		if (c<sigma_mid) {
			pos -= rank_B[pointer](pos);
			pointer = pointer * 2;
			sigma_end = sigma_mid;
		} else {
			pos = rank_B[pointer](pos);
			pointer = pointer * 2 + 1;
			sigma_start = sigma_mid;
		}
	}
	return pos;
}

template <typename DataType, class BitVector, class RankSupport>
DataType wavelet_tree<DataType, BitVector, RankSupport>::extract(size_t pos) const {
	size_t pointer = 1;
	uint8_t sigma_start = 0;
	uint8_t sigma_end = m_sigma;
	while (sigma_end-sigma_start>1) {
		if (sigma_end-sigma_start==1) break;
		uint8_t sigma_mid = (sigma_start+sigma_end)/2;
		if (B[pointer][pos]==0) {
			pos -= rank_B[pointer](pos);
			pointer = pointer * 2;
			sigma_end = sigma_mid;
		} else {
			pos = rank_B[pointer](pos);
			pointer = pointer * 2 + 1;
			sigma_start = sigma_mid;
		}
	}
	return (DataType) (sigma_end+sigma_start)/2;
}

template <typename DataType, class BitVector, class RankSupport>
void wavelet_tree<DataType, BitVector, RankSupport>::getBounds(size_t i, size_t j, DataType c, size_t &lesser, size_t &greater) {
	size_t pointer = 1;
	uint8_t sigma_start = 0;
	uint8_t sigma_end = m_sigma;

	lesser = greater = 0;

	while (sigma_end-sigma_start>1) {
		if (sigma_end-sigma_start==1) break;
		uint8_t sigma_mid = (sigma_start+sigma_end)/2;

		size_t a1 = rank_B[pointer](i);
		size_t b1 = rank_B[pointer](j);
		size_t a0 = i-a1;
		size_t b0 = j-b1;

		if (c<sigma_mid) {
			greater += b1-a1;
			i = a0;
			j = b0;
			pointer = pointer * 2;
			sigma_end = sigma_mid;
		} else {
			lesser += b0-a0;
			i = a1;
			j = b1;
			pointer = pointer * 2 + 1;
			sigma_start = sigma_mid;
		}
	}
}

#endif

