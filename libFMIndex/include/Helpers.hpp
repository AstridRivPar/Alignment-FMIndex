#pragma once
#include <stdio.h>
#include <sstream>
#include <memory>



class Range{
        
    public:
        size_t begin;
        size_t end;
        /**
         * Default constructor, initializes an empty range
         */
        Range() : begin(0), end(0) {}
        /**
         * Constructor
         * @param b, the beginning of the range
         * @param e, the end of the range (non-inclusive)
         */
        Range(size_t b, size_t e) : begin(b), end(e) {
        }
        /**
         * Constructor
         * @param size, size of the whole interval
         */
        Range(size_t size) : begin(0), end(size) {
        }

        size_t getBegin() const {
        return begin;
        }
        size_t getEnd() const {
            return end;
        }
        /**
         * Check if this range is empty
         * @returns true if the range is empty, false otherwise
         */
        bool empty() const {
            return end <= begin;
        }
        /**
         * Gets the width of the range (end - begin)
         * @returns the width of this range
         */
        size_t width() const {
            return (empty()) ? 0 : end - begin;
        }

        /**
         * Operator overloading, two ranges are equal if their begin and end field
         * are equal
         */
        bool operator==(const Range& o) const {
            return o.getBegin() == begin && o.getEnd() == end;
        }
        // friend std::ostream& operator<<(std::ostream& os, const Range& r);
        friend std::ostream& operator<<(std::ostream& os, const Range& r) {
            os << "[" << r.begin << ", " << r.end << ")";
            return os;
        }

};

class SARangePair {
    
  public:
    Range rangeSA;    // the range over the suffix array
    Range rangeSARev; // the range over the suffix array of the reversed text
    std::string pMatch = ""; //only for debug
	size_t length = 0; // how many characters have been matched
 	 
    /**
     * Default constructor, creates two empty ranges
     */
    SARangePair() : rangeSA(Range()), rangeSARev(Range()) {
    }

    SARangePair(Range rangeSA, Range rangeSARev, std::string pMatch, size_t length)
        : rangeSA(rangeSA), rangeSARev(rangeSARev), pMatch(pMatch), length(length){
    }

    SARangePair(const SARangePair &range)
        :rangeSA(range.rangeSA), rangeSARev(range.rangeSARev), pMatch(range.pMatch), length(range.length){}

    SARangePair(size_t size)
        : rangeSA(size), rangeSARev(size) {
    }

    const Range& getRangeSA() const {
        return rangeSA;
    }

    const Range& getRangeSARev() const {
        return rangeSARev;
    }
    /**
     * @returns true if the ranges are empty, false otherwise
     */
    bool empty() const {
        return rangeSA.empty();
    }

    size_t width() const {
        return rangeSA.width();
    }

    bool isCorrect() const{
        return rangeSA.width() == rangeSARev.width() && !empty();
    }
    size_t get_length() const {
        return length;
    }
    /**
     * Operator overloading
     * @returns true if this is equal to rhs
     */
    bool operator==(const SARangePair& o) const {
        // only the first range matters as the ranges imply each other
        return o.getRangeSA() == rangeSA;
    }
    friend std::ostream& operator<<(std::ostream& os, const SARangePair& r) {
        os << "BW:" << r.getRangeSA()<<" FW:"<< r.getRangeSARev();
        return os;
    }
}; 

class Result{
    public:
        SARangePair range;
        size_t row = 0;
        std::vector<std::string> alignments;
        std::string match;
        int k = 0;
        size_t idx  = 0; // location in index
        // Result (Result res){

        // }
        Result(const SARangePair &r, std::vector<std::string> alignments, int k): range(r), k(k),
        alignments(std::move(alignments)){}

        Result(std::vector<std::string> alignments, std::string &match,int k): alignments(std::move(alignments)), 
        k(k), match(std::move(match)){}

        Result(const SARangePair &r): range(r){
            // alignments = "";
            alignments.reserve(10);
            match = "";
        }

        Result(const SARangePair &r, std::vector<std::string> alignments, std::string match,int k): range(r), alignments(std::move(alignments)), 
        k(k), match(std::move(match)){}

        friend std::ostream& operator<<(std::ostream& os, const Result &res) {
            os << res.match << ",{";
            size_t al_size = res.alignments.size() - 1;
            for (size_t i = 0; i< al_size; i++){
                os << res.alignments[i] <<",";
            }
            os << res.alignments[al_size] << "}";
            return os;
        }
        bool operator == (const Result &res){ 
            if (range == res.range) 
                return true; 
            return false; 
        } 
        
        
};

struct Node{
    Node(const SARangePair &range_, char c_, const int &row_): 
    range(range_),
    row(row_),
    c(c_){}

    public:
        SARangePair range;
        int row = 0;
        char c;
};

template<typename Iterator>
struct State {
    int row;
    int col;
    Iterator model_it;
    Iterator log_it;
    std::string path;

    State(int r, int c, Iterator m_it, Iterator l_it, std::string p)
        : row(r), col(c), model_it(m_it), log_it(l_it), path(p) {}
};