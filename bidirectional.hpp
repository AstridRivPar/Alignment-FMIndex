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
/*! \file bidirectional.hpp
  	\brief This file contains several methods that search bidirectionally
	in a given Index.
	\author Thomas Schnattinger
*/
#ifndef INCLUDED_BIDIRECTIONAL_HPP
#define INCLUDED_BIDIRECTIONAL_HPP

#include <sdsl/suffixarrays.hpp>
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <cstring>
#include <stack>

/*! \mainpage Index Page - Bidirectional search algorithms

\section intro_sec Introduction

This is a set of files written for the Diploma Thesis
"Bidirektionale indexbasierte Suche in Texten"
written by Thomas Schnattinger (thomas.schnattinger@uni-ulm.de)
at University of Ulm, 12.01.2010

\section install_sec Installation

To compile you need
\li c++ boost library >= 1.3.5
  http://www.boost.org/
\li succinct data structure library (sdsl)
  http://www.uni-ulm.de/in/theo/research/sdsl
  (library used for the tests in the thesis is 0.7.5)

  Then adjust the Makefile to your working directory

  For the dynamic libraries you have to set
\verbatim export LD_LIBRARY_PATH=/path/to/the/sdsl/libraries \endverbatim

\section using_sec Using the algorithms

  Three executables are produced by "make"

  \li test_indexes. Small program where each index method is used once.
  \li createindex. Create an index.
  \verbatim
Creates an index from the given file

usage: ./createindex inputfile #index(1..6) indexfile

 # | Index
---+-----------------------------------
 1 | index_sa_text_psi
 2 | index_sa_isa
 3 | index_sa_text_occ<>
 4 | index_csa_psi_text<>
 5 | index_csa_fmi_text<>
 6 | index_bidirectional_waveletindex<>
\endverbatim
  \li starttestcases. Test case
  \verbatim
Starts the test cases mentioned in the thesis.

usage: ./starttestcases indexfile [durchgaenge (default=1)]
\endverbatim

\section example_sec Example
\li Create the Bidirectional Wavelet Index for the file "yeast.raw"
\verbatim
user@pc ~ $ ./createindex yeast.raw 6 yeast.index
\endverbatim
\li Start test cases for the created index
\verbatim
user@pc ~ $ ./testcase yeast.index
Found index_bidirectional_waveletindex<>
Speicherplatz: 9.867.588 Bytes
* Suche nach "hairpin1 = (stem:=N{20,50}) (loop:=NNN) ^stem"
  Muster 4 mal gefunden...
  Zeit zum Suchen: 924 ms.
* Suche nach "hairpin2 = (stem:=N{10,50}) (loop:=GGAC) ^stem"
  Muster 2 mal gefunden...
  Zeit zum Suchen: 0 ms.
* Suche nach "hairpin4 = (stem:=N{10,15}) (loop:=GGAC[0,0,1]) ^stem"
  Muster 24 mal gefunden...
  Zeit zum Suchen: 12 ms.
* Suche nach "hloop(5) = (stem:=N{15,20}) (loop:=N{5}) ^stem"
  Muster 62 mal gefunden...
  Zeit zum Suchen: 852 ms.
* Suche nach "acloop(5) = (stem:=N{15,20}) (loop:=(A|C){5}) ^stem"
  Muster 1 mal gefunden...
  Zeit zum Suchen: 28 ms.
* Suche nach "acloop(10) = (stem:=N{15,20}) (loop:=(A|C){10}) ^stem"
  Muster 0 mal gefunden...
  Zeit zum Suchen: 12 ms.
* Suche nach "acloop(15) = (stem:=N{15,20}) (loop:=(A|C){15}) ^stem"
  Muster 0 mal gefunden...
  Zeit zum Suchen: 8 ms.
\endverbatim
*/

using namespace std;

/*!
 * \brief class to parameterize a bidirectional search
 * specifiying which nucleotides should be complementary
 */
class complement {
	private:
		bool data[256][256];
	public:
		/*!
		 * Construct the object from the given rules.
		 *
		 * \param[in] rules Comma separated rules. One rule has the form 
		 * "X-Y", where X and Y are compelementary bases.
		 * \param[in] char2comp Must point to the char2comp of the index.
		 */
		complement(const char* rules, const unsigned char* char2comp) {
			memset(data, 0, sizeof(data));
			size_t n = strlen(rules);
			for (int pos = 0; pos<n; pos += 4) {
				uint8_t a = char2comp[rules[pos]], b = char2comp[rules[pos+2]];
				data[a][b] = data[b][a] = true;
			}
		}

		/*!
		 * operator() so that the object can be used like a method.
		 * \param[in] a nucleotide (in the range [1..Sigma])
		 * \param[in] b nucleotide (in the range [1..Sigma])
		 * \return true if the nucleotides are complementary, according to the given rules
		 */
		bool operator() (unsigned char a, unsigned char b) const {
			return data[a][b];
		}
};

/*!
 * Print the result of a search. The purpose of this method is to measure the time 
 * used by the extract_sa method of the index
 * \param[in] index index
 * \param[in] s search state to report
 * \param out output stream to write the results to
 */
template<class Index>
void report_hit(Index &index, typename Index::searchstate s, std::ostream &out) {
	int_vector<> sa;
	index.extract_sa(s.occ_begin, s.occ_end, sa);
	out << "[" << (s.occ_begin+1) << ".." << s.occ_end << "] = {";
	for (size_t i=0; i<sa.size(); i++) {
		if (i>0) out << ",";
		out << (uint64_t)sa[i];
	}
	out << "}" << endl;
}

/*!
 * This method is searches for a given pattern. It is not mentioned in the thesis. 
 * It searches for a given pattern and starts bidirectional search afterwards. 
 * It prints the maximal loops to stdout
 *
 * \param[in] index Search in this Index 
 * \param[in] pattern pattern to search for
 * \param[in] pattern_len length of pattern
 */
template<class Index>
void locate_hairpin(Index &index, const unsigned char *pattern, size_t pattern_len) {

	typedef typename Index::searchstate searchstate;

	/* specify the complementary nukleotides */
	complement complement("A-T,C-G,G-T", index.char2comp);

	/* search the initial pattern by forward search */
	searchstate s;
	index.init_search_state(s);
	for (size_t i=0; i<pattern_len; i++) index.forward_search(index.char2comp[pattern[i]], s);

	const uint8_t sigma(index.sigma);

	/* the stack used by the depth-first search */
	stack<searchstate> stack;
	stack.push(s);

	while (!stack.empty()) {

		searchstate s = stack.top();
		stack.pop();

		/* don't procede an empty interval */
		if (s.occ_begin==s.occ_end) continue;

		bool matched = false;

		/* search backward */
		// try characters from highest to lowest to get a left-to-right depth-first traversal
		for (uint8_t c=sigma-1; c>0; c--) {
			searchstate t = s;
			index.backward_search(c, t);
			if (t.occ_end > t.occ_begin) {

				/* search forward */
				// try characters from highest to lowest to get a left-to-right depth-first traversal
				for (uint8_t c2=sigma-1; c2>0; c2--) {
					if (complement(c, c2)) {
						searchstate u = t;
						index.forward_search(c2, u);
						if (u.occ_end > u.occ_begin) {
							matched = true;
							stack.push(u);
						}
					}
				}
	
			}
		}

		/* if the given pattern can't be extended to both sides then this is the biggest hairpin loop */
		if (!matched) {
			// not matched any character in forward direction 
			// -> "maximal hairpin loop found"
			cout << "maximal hairpin loop: ";
			cout << "[" << s.occ_begin+1 << ".." << s.occ_end << "]";
			cout << " (Anzahl: "<<s.occ_end-s.occ_begin<<", Länge: "<<s.length<<")\n";
		}
		
	}
}



/*!
 * Search for the following pattern in the given index and print the results to out
 *
 * \verbatim hairpin1 = (stem:=N{20,50}) (loop:=NNN) ^stem\endverbatim
 *
 * \param[in] index Search in this Index 
 * \param[in,out] out Print the result to this OutputStream
 * \return number of hits
 */
template<class Index>
size_t locate_testpattern_hairpin1(Index &index, std::ostream &out) {
	printf("* Suche nach \"hairpin1 = (stem:=N{20,50}) (loop:=NNN) ^stem\"\n");
	typedef typename Index::searchstate searchstate;
	complement complement("A-T,C-G,G-T", index.char2comp);

	stack<searchstate> stack;
	searchstate f;
	index.init_search_state(f);

	// search (loop:=NNN)
	unsigned char nukleotides[4] = {'A', 'C', 'T', 'G'};
	for (int x=0; x<4; x++) {
		searchstate sx = f;
		index.forward_search(index.char2comp[nukleotides[x]], sx);
		if (sx.occ_begin>=sx.occ_end) continue;
		for (int y=0; y<4; y++) {
			searchstate sy = sx;
			index.forward_search(index.char2comp[nukleotides[y]], sy);
			if (sy.occ_begin>=sy.occ_end) continue;
			for (int z=0; z<4; z++) {
				searchstate sz = sy;
				index.forward_search(index.char2comp[nukleotides[z]], sz);
				if (sz.occ_begin>=sz.occ_end) continue;
				stack.push(sz);
			}
		}
	}

	size_t min_stem_length = 20;
	size_t max_stem_length = 50;

	size_t hits = 0;

	const uint8_t sigma(index.sigma);

	while (!stack.empty()) {
		

		searchstate s = stack.top();
		stack.pop();

		if (s.occ_begin==s.occ_end) continue;
		if ((s.length-3)/2 >= max_stem_length) {
			++hits;
			report_hit(index, s, out);
			continue;
		}

		bool matched = false;

		// search backward
		// try characters from highest to lowest to get a left-to-right depth-first traversal
		for (uint8_t c=sigma-1; c>0; c--) {
			searchstate t = s;
			index.backward_search(c, t);
			if (t.occ_end > t.occ_begin) {
				// search forward
				// try characters from highest to lowest to get a left-to-right depth-first traversal
				for (uint8_t c2=sigma-1; c2>0; c2--) {
					if (complement(c, c2)) {
						searchstate u = t;
						index.forward_search(c2, u);
						if (u.occ_end > u.occ_begin) {
							matched = true;
							stack.push(u);
						}
					}
				}
	
			}
		};

		if (!matched && (s.length-3)/2 >= min_stem_length) {
			report_hit(index, s, out);
			++hits;
		}
		
	}
	return hits;
}

/*!
 * Search for the following pattern in the given index and print the results to out
 *
 * \verbatim hairpin2 = (stem:=N{10,50}) (loop:=GGAC) ^stem \endverbatim
 *
 * \param[in] index Search in this Index 
 * \param[in,out] out Print the result to this OutputStream
 * \return number of hits
 */
template<class Index>
size_t locate_testpattern_hairpin2(Index &index, std::ostream &out) {
	printf("* Suche nach \"hairpin2 = (stem:=N{10,50}) (loop:=GGAC) ^stem\"\n");
	typedef typename Index::searchstate searchstate;
	complement complement("A-T,C-G,G-T", index.char2comp);

	stack<searchstate> stack;
	searchstate f;
	index.init_search_state(f);

	// search (loop:=GGAC)
	const unsigned char *pattern = (const unsigned char*)"GGAC";
	for (size_t i=0; i<4; i++) index.forward_search(index.char2comp[pattern[i]], f);
	stack.push(f);

	size_t min_stem_length = 10;
	size_t max_stem_length = 50;

	size_t hits = 0;

	const uint8_t sigma(index.sigma);

	while (!stack.empty()) {

		searchstate s = stack.top();
		stack.pop();

		if (s.occ_begin==s.occ_end) continue;
		if ((s.length-4)/2 >= max_stem_length) {
			++hits;
			report_hit(index, s, out);
			continue;
		}

		bool matched = false;

		// search backward
		// try characters from highest to lowest to get a left-to-right depth-first traversal
		for (uint8_t c=sigma-1; c>0; c--) {
			searchstate t = s;
			index.backward_search(c, t);
			if (t.occ_end > t.occ_begin) {
				// search forward
				// try characters from highest to lowest to get a left-to-right depth-first traversal
				for (uint8_t c2=sigma-1; c2>0; c2--) {
					if (complement(c, c2)) {
						searchstate u = t;
						index.forward_search(c2, u);
						if (u.occ_end > u.occ_begin) {
							matched = true;
							stack.push(u);
						}
					}
				}
	
			}
		}

		if (!matched && (s.length-3)/2 >= min_stem_length) {
			++hits;
			report_hit(index, s, out);
		}
		
	}
	return hits;
}

/*!
 * Search for the following pattern in the given index and print the results to out
 *
 * \verbatim hairpin4 = (stem:=N{10,15}) (loop:=GGAC[0,0,1]) {1 insertion} ^stem\endverbatim
 *
 * \param[in] index Search in this Index 
 * \param[in,out] out Print the result to this OutputStream
 * \return number of hits
 */
template<class Index>
size_t locate_testpattern_hairpin4(Index &index, std::ostream &out) {
	printf("* Suche nach \"hairpin4 = (stem:=N{10,15}) (loop:=GGAC[0,0,1]) ^stem\"\n");
	typedef typename Index::searchstate searchstate;
	complement complement("A-T,C-G,G-T", index.char2comp);

	stack<searchstate> stack;
	searchstate f;
	index.init_search_state(f);

	// search (loop:=GGAC)
	const unsigned char *pattern = (const unsigned char*)"GGAC";
	for (size_t i=0; i<4; i++) index.forward_search(index.char2comp[pattern[i]], f);
	stack.push(f);

	unsigned char nukleotides[4] = {'A', 'C', 'T', 'G'};
	for (int x=0; x<4; x++) {
		for (int k=0; k<5; k++) {
			index.init_search_state(f);
			for (int i=0; i<k; i++) index.forward_search(index.char2comp[pattern[i]], f);
			index.forward_search(index.char2comp[nukleotides[x]], f);
			for (int i=k; i<4; i++) index.forward_search(index.char2comp[pattern[i]], f);
			stack.push(f);
		}
	}

	size_t min_stem_length = 10;
	size_t max_stem_length = 15;

	size_t hits = 0;

	const uint8_t sigma(index.sigma);

	while (!stack.empty()) {

		searchstate s = stack.top();
		stack.pop();

		if (s.occ_begin==s.occ_end) continue;
		if ((s.length-4)/2 >= max_stem_length) {
			++hits;
			report_hit(index, s, out);
			continue;
		}

		bool matched = false;

		// search backward
		// try characters from highest to lowest to get a left-to-right depth-first traversal
		for (uint8_t c=sigma-1; c>0; c--) {
			searchstate t = s;
			index.backward_search(c, t);
			if (t.occ_end > t.occ_begin) {
				// search forward
				// try characters from highest to lowest to get a left-to-right depth-first traversal
				for (uint8_t c2=sigma-1; c2>0; c2--) {
					if (complement(c, c2)) {
						searchstate u = t;
						index.forward_search(c2, u);
						if (u.occ_end > u.occ_begin) {
							matched = true;
							stack.push(u);
						}
					}
				}
	
			}
		}

		if (!matched && (s.length-3)/2 >= min_stem_length) {
			++hits;
			report_hit(index, s, out);
		}
		
	}
	return hits;
}

/*!
 * Search for the following pattern in the given index and print the results to out
 *
 * \verbatim sequence = CAGUAGAAA \endverbatim
 *
 * \param[in] index Search in this Index 
 * \param[in,out] out Print the result to this OutputStream
 * \return number of hits
 */
template<class Index>
size_t locate_testpattern_sequence(Index &index, std::ostream &out) {
	printf("* Suche nach \"sequence = CAGUAGAAA\"\n");
	typedef typename Index::searchstate searchstate;
	searchstate f;
	index.init_search_state(f);
	const unsigned char* pattern = (const unsigned char*) "CAGUAGAAA";

	for (size_t i=0; f.occ_begin<f.occ_end && pattern[i]!=0; i++) {
		index.forward_search(index.char2comp[pattern[i]], f);
	}

	return f.occ_end-f.occ_begin;
}

/*! Help function. Fills the stack recursively with the searchstates of every possible
 * loop combination of length looplength specified by characters&characters_len
 * \param[in] index
 * \param[out] stack Stack where the results are pushed on
 * \param[in] f original searchstate
 * \param[in] characters characters in the loop
 * \param[in] characters_len size of characters
 * \param[in] looplength length of the loop
 */
template<class Index>
void fillstack_hloop(
		Index &index, 
		stack<typename Index::searchstate> &stack, 
		typename Index::searchstate f, 
		const unsigned char* characters, 
		size_t characters_len, 
		size_t looplength) 
{

	if (looplength==0) {

		stack.push(f);

	} else {

		for (int x=0; x<characters_len; x++) {
			typename Index::searchstate s = f;
			index.forward_search(index.char2comp[characters[x]], s);
			if (s.occ_begin>=s.occ_end) continue;
			fillstack_hloop(index, stack, s, characters, characters_len, looplength-1);
		}

	}

}

/*!
 * Search for the following pattern in the given index and print the results to out
 *
 * \verbatim hloop(looplength) = (stem:=N{15,20}) (loop:=N{looplength}) ^stem \endverbatim
 *
 * \param[in] index Search in this Index 
 * \param[in] looplength Length of the loop
 * \param[in,out] out Print the result to this OutputStream
 * \return number of hits
 */
template<class Index>
size_t locate_testpattern_hloop(Index &index, size_t looplength, std::ostream &out) {
	printf("* Suche nach \"hloop(%d) = (stem:=N{15,20}) (loop:=N{%d}) ^stem\"\n", looplength, looplength);
	typedef typename Index::searchstate searchstate;
	complement complement("A-T,C-G,G-T", index.char2comp);

	stack<searchstate> stack;
	searchstate f;
	index.init_search_state(f);

	// search (loop:=NNN)
	const unsigned char nukleotides[4] = {'A', 'C', 'T', 'G'};
	fillstack_hloop(index, stack, f, nukleotides, 4, looplength);

//	printf("  stack size = %d\n", stack.size());
	size_t min_stem_length = 15;
	size_t max_stem_length = 20;

	size_t hits = 0;

	const uint8_t sigma(index.sigma);

	while (!stack.empty()) {
		

		searchstate s = stack.top();
		stack.pop();

		if (s.occ_begin==s.occ_end) continue;
		if ((s.length-3)/2 >= max_stem_length) {
			++hits;
			report_hit(index, s, out);
			continue;
		}

		bool matched = false;

		// search backward
		// try characters from highest to lowest to get a left-to-right depth-first traversal
		for (uint8_t c=sigma-1; c>0; c--) {
			searchstate t = s;
			index.backward_search(c, t);
			if (t.occ_end > t.occ_begin) {
				// search forward
				// try characters from highest to lowest to get a left-to-right depth-first traversal
				for (uint8_t c2=sigma-1; c2>0; c2--) {
					if (complement(c, c2)) {
						searchstate u = t;
						index.forward_search(c2, u);
						if (u.occ_end > u.occ_begin) {
							matched = true;
							stack.push(u);
						}
					}
				}
	
			}
		};

		if (!matched && (s.length-3)/2 >= min_stem_length) {
			++hits;
			report_hit(index, s, out);
		}
		
	}
	return hits;
}

/*!
 * Search for the following pattern in the given index and print the results to out
 *
 * \verbatim acloop(looplength) = (stem:=N{15,20}) (loop:=N{looplength}) ^stem \endverbatim
 *
 * \param[in] index Search in this Index 
 * \param[in] looplength Length of the loop
 * \param[in,out] out Print the result to this OutputStream
 * \return number of hits
 */
template<class Index>
size_t locate_testpattern_acloop(Index &index, size_t looplength, std::ostream &out) {
	printf("* Suche nach \"acloop(%d) = (stem:=N{15,20}) (loop:=(A|C){%d}) ^stem\"\n", looplength, looplength);
	typedef typename Index::searchstate searchstate;
	complement complement("A-T,C-G,G-T", index.char2comp);

	stack<searchstate> stack;
	searchstate f;
	index.init_search_state(f);

	// search (loop:=NNN)
	const unsigned char nukleotides[4] = {'A', 'C'};
	fillstack_hloop(index, stack, f, nukleotides, 2, looplength);

//	printf("  stack size = %d\n", stack.size());
	size_t min_stem_length = 15;
	size_t max_stem_length = 20;

	size_t hits = 0;

	const uint8_t sigma(index.sigma);

	while (!stack.empty()) {
		

		searchstate s = stack.top();
		stack.pop();

		if (s.occ_begin==s.occ_end) continue;
		if ((s.length-3)/2 >= max_stem_length) {
			++hits;
			report_hit(index, s, out);
			continue;
		}

		bool matched = false;

		// search backward
		// try characters from highest to lowest to get a left-to-right depth-first traversal
		for (uint8_t c=sigma-1; c>0; c--) {
			searchstate t = s;
			index.backward_search(c, t);
			if (t.occ_end > t.occ_begin) {
				// search forward
				// try characters from highest to lowest to get a left-to-right depth-first traversal
				for (uint8_t c2=sigma-1; c2>0; c2--) {
					if (complement(c, c2)) {
						searchstate u = t;
						index.forward_search(c2, u);
						if (u.occ_end > u.occ_begin) {
							matched = true;
							stack.push(u);
						}
					}
				}
	
			}
		};

		if (!matched && (s.length-3)/2 >= min_stem_length) {
			++hits;
			report_hit(index, s, out);
		}
		
	}
	return hits;
}

#endif
