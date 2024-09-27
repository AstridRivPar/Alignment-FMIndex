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
/*! \file indexes.hpp
    \brief The base class of the index classes.
	\author Thomas Schnattinger
*/
#ifndef INCLUDED_BIDIRECTIONAL_INDEX
#define INCLUDED_BIDIRECTIONAL_INDEX

#include <iostream>
#include <cassert>
#include <utility>
//! The base class of the index classes.
class index {
	protected:

		//! Super Constructor. Invoked by all subclass constructors
		index() : sigma(m_sigma), C(m_C), char2comp(m_char2comp), comp2char(m_comp2char) {}

		//! Initialize sigma, char2comp, comp2char from the given text
		void setText(const unsigned char *str, size_t len);
	
		//! Alphabet size
		uint8_t m_sigma;
		
		/*! Array specifying at position i the number of characters in 
		 *  the text which are strictly smaller than character i 
		 */
		size_t m_C[257];

		//! Convert between ASCII and working alphabet
		unsigned char m_char2comp[256];

		//! Convert between working and ASCII alphabet
		unsigned char m_comp2char[256];

		//! The lenth of the original text (including final null character)
		size_t m_size;

		//! Serialize the attributes of the base class to an output stream
		size_t superserialize(std::ostream &out) const;

		//! Load the attributes of the base class from an input stream
		void superload(std::istream &in);

	

	public:

		//! Alphabet size 
		const uint8_t &sigma;

		const size_t *C;

		
		
		//! Convert between ASCII and working alphabet
		const unsigned char *char2comp;
		
		//! Convert between working and ASCII alphabet
		const unsigned char *comp2char;

		//! The lenth of the original text (including final null character)
		size_t size() const;

		std::pair<unsigned char*, uint8_t> get_alphabet();

};

std::pair<unsigned char*,uint8_t> index::get_alphabet(){
	
	unsigned char * p = m_comp2char;
	return std::make_pair(p, m_sigma);
}

size_t index::superserialize(std::ostream &out) const {
	size_t written_bytes = 0;
	size_t wb   = sizeof(m_char2comp[0])*256;
	out.write((char*)m_char2comp, wb);
	written_bytes += wb; 
	wb             = sizeof(m_comp2char[0])*256;
	out.write((char*)m_comp2char, wb);
	written_bytes += wb; 
	wb             = sizeof(m_C[0])*257;
	out.write((char*)m_C, wb);
	written_bytes += wb; 
	out.write((char*)&m_sigma, sizeof(m_sigma));
	wb += sizeof(m_sigma);
	out.write((char*)&m_size, sizeof(m_size));
	wb += sizeof(m_size);
	return written_bytes;
}

void index::superload(std::istream &in) {
	in.read((char*)m_char2comp, sizeof(m_char2comp[0])*256);
	in.read((char*)m_comp2char, sizeof(m_comp2char[0])*256);
	in.read((char*)m_C, sizeof(m_C[0])*257);
	in.read((char*)&m_sigma, sizeof(m_sigma));
	in.read((char*)&m_size, sizeof(m_size));
}

// written by Simon Gog
void index::setText(const unsigned char *str, size_t len) {
	m_size = len;

	for(uint16_t i=0; i<256; ++i)
		m_char2comp[i] = m_comp2char[i] = m_C[i] = 0;
	m_C[256] = 0;
	if( str == NULL or len ==0 )
		return;
	const unsigned char *p = str;
	for(size_t i=0;i<len-1;++i){
		m_char2comp[*p++] = 1;
	}   
	assert(m_char2comp[0] == 0); 
	m_char2comp[*p] = 1;
	assert(m_char2comp[0] == 1); 
	size_t value = 0;
	for(uint16_t i=0; i<256; ++i)
		if(m_char2comp[i])
			m_char2comp[i] = value++;
	m_sigma = value;
	for(uint16_t i=0; i<256; ++i)
		if(m_char2comp[i])
			m_comp2char[m_char2comp[i]] = i;
	for(uint16_t i=0; i<257; ++i) m_C[i]=0;
	for(size_t i=0; i<len; ++i) ++m_C[m_char2comp[str[i]]];
	for(uint16_t i=256;i>0; --i) m_C[i] = m_C[i-1];
	m_C[0] = 0;
	for(uint16_t i=1; i<257; ++i) m_C[i] += m_C[i-1];   
	if(m_C[256]!=len){
		std::cerr<<"m_C[256]="<<m_C[256]<<" "<<len<<std::endl;
	}   
	assert(m_C[256]==len);
	assert(m_C[m_sigma+1]==len);
}

size_t index::size() const {
	return m_size;
}

// includes...
#include "index_bidirectional_waveletindex.hpp"



#endif // end file 
