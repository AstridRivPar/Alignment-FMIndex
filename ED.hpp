#include "index_bidirectional_waveletindex.hpp"
#include "bandmatrix.h"
#include "customtypedefs.h"
#include <unordered_map>
#include <queue>
#include <deque>
#include <algorithm> 
#include <limits>
#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

class ED{
    private:
       
        typedef index_bidirectional_waveletindex<> Index; 
        typedef Index::searchstate searchstate;

        Index *FMIndex;
        const unsigned char * text;
        
        //direction*
        size_t (Index::*extraChar)(unsigned char, searchstate&); //// returns |f-g|
        searchstate (Index::*fbsPtr)(unsigned char, const searchstate&); //returns a new searchstate(SA pair) 
        bool fw = true;
        vector<Node> nodesToCheck;
        Matrix M;
        //alignment with replace or without replace
        // void (Matrix::*alignFunc)(int, int, std::unique_ptr<std::string>&);
       
        vector<Result> *sols;

        vector<int> foundPerK;

        

    public:
        std::map<std::string, std::vector<Result>> results;

        ED():FMIndex(NULL), text(NULL){};
        /**
         * Constructor 
         * @param index FM-Index
         * @param txt Text
         */
        ED(Index *index, const unsigned char* txt):FMIndex(index), text(txt){}
        /**
         * Match exactly a trace or a part of the trace
         * @param res stores the alignment 
         * @param offset for a partial alignment
         * @param st Iterator or reverse iterator (start)
         * @param ed ITerator or reverse iterator (end)
         */

        
        template <typename It>
        bool exactMatching(Result &res, int offset, It st, It ed);

        void setIdx(Index *index, const unsigned char*txt);
        /**
        * Match a trace with maxED = 0
        */
        bool matchExactly(const std::string &query);

        int CheckConf(const std::string &trace, int maxED, bool complete);

        std::map<std::string, std::vector<Result>> CompleteAligment(const vector<std::string> &queries, int maxED, bool isFw, bool complete, bool replace);

        void setDirection( bool fw);
        void BFSearch(const std::string &query, int maxED, bool complete, const searchstate &s);
        void BuildSolution(searchstate &s, int row, int col);
        /**
         * Returns the optimal value
         */
        int TakeAllOptimal();
        void BuildSolution(Result &res, int row, int col);
        void pushChildren(const searchstate &s, int row);
        /**
         * Get number of optimal traces for a given trace
         * @param query Query string for which we want to return the number of optimal traces
         */
        int getNoOpt(const std::string &query);
        void getBest(std::ostream& of);
        vector<int> getFoundPerK();
};

void ED::setIdx(Index *index, const unsigned char*txt){
    this->FMIndex = index;
    this->text = txt;
}
int ED::getNoOpt(const std::string &trace){
    if (results.find(trace) != results.end()){
        sols = &results[trace];
        if(!sols->empty())
            return sols->size();
    }
    return -1;
}
vector<int> ED::getFoundPerK(){
    return foundPerK;
}
// }
void ED::setDirection( bool isFw){
    fw = isFw;
    if (fw){
        extraChar = &Index::forward_search;
        fbsPtr = &Index::fs;

    }else{
        extraChar = &Index::backward_search;
        fbsPtr = &Index::bs;

    }

}

void ED::pushChildren(const searchstate &s, int row){
    const auto &[alphabet_begin, alphabet_length] =  FMIndex->get_alphabet(); 
        
    for (int i = alphabet_length; i >=2; --i){
        searchstate sp = (FMIndex->*fbsPtr)(i, s);
        if (sp.occ_end - sp.occ_begin > 0){
            nodesToCheck.emplace_back(sp, *(alphabet_begin + i), row + 1);
        }
    }
}

int ED::TakeAllOptimal(){
    if (sols->empty()) return -1;
    
    auto minK = sols->back().k;
    for (auto it = sols->begin();it < sols->end(); it++){
        
        it->match = std::unique_ptr<std::string> (new std::string());

        auto firstOcc = FMIndex ->extract(it->s.occ_begin, it->s.occ_end);
        for(size_t j = 1; j < it->s.length - 1; j++){ //revisar esto
            it->match->push_back(text[firstOcc + j]);
        }

        if (it->k > minK) return minK;
    }
    return minK;
    
}

int ED::CheckConf(const std::string &query, int maxED, bool complete){
    //init
    nodesToCheck.clear();
    nodesToCheck.reserve(100);
    
   
    //if trace is already aligned       
    
    results[query] = vector<Result>();
    sols = &results[query];
    //match exactly
    if (maxED == 0) {
        bool isFound = matchExactly(query);
        
        return isFound? 0: -1;
    }
    //initialize matrix
    M.init(query.size(), maxED);
    

    searchstate init;
    FMIndex->init_search_state(init);

    if(complete){
        (FMIndex->*extraChar)(FMIndex->char2comp[','], init);
    }
    BFSearch(query, maxED, complete, init);
    
    return TakeAllOptimal();

    
}
std::map<std::string, std::vector<Result>> ED::CompleteAligment(const vector<std::string> &queries, int maxED, bool isFw, bool complete, bool replace){

    setDirection(isFw);
    M.setAligmentFunc(replace);
    foundPerK.resize(maxED + 2);
    
    for (const auto & query: queries){
        int pK = -1;
        if (results.find(query) != results.end()){
            sols = &results[query];
            if(!sols->empty())
                pK = sols->back().k;
        }else{
            pK = CheckConf(query, maxED, complete);
        }
        
        foundPerK[pK != -1 ? pK : maxED + 1]++;
        
    }

    // // results[trace] = std::move(sols);
    return std::move(results);

    
}
void ED::BuildSolution(searchstate &s, int row, int col){
    std::unique_ptr<std::string> al (new std::string());
    
    M.getAlign(row, col, al);
    if (fw) std::reverse(al->begin(), al->end());
    // (M.*alignFunc)(row, col, al);
    sols->emplace_back(s,al,M(row, col));
    
}
void ED::BuildSolution(Result &res, int row, int col){
    // (M.*alignFunc)(row, col, res.alignment);
    M.getAlign(row, col, res.alignment);
    if (fw) std::reverse(res.alignment->begin(), res.alignment->end());
    sols->push_back(std::move(res));
}
void ED::BFSearch(const std::string &query, int maxED, bool complete, const searchstate &s){
   
    int qSize = query.size();
    pushChildren(s, 0);
    while (!nodesToCheck.empty()){
        auto [sp, row, c] = nodesToCheck.back();
        nodesToCheck.pop_back();
        
        int minimalEDOfRow = M(row, 0);
        for(int col = 1; col <= qSize; col++){ 
            bool notMatch = fw ? query[col - 1] != c : query[qSize - col] != c;
            M.updateMatrix(notMatch, row, col);
            minimalEDOfRow = std::min(minimalEDOfRow, M(row, col));

        }
        if (minimalEDOfRow > maxED) continue;
        if (row < qSize  && minimalEDOfRow == maxED){
            //exact match
            for(int j = row - maxED; j < M.cols; j++){
                if (M(row, j) == minimalEDOfRow){
                    Result res(sp);
                    bool isFound = fw? exactMatching(res, j, query.begin(), query.end()):exactMatching(res, j, query.rbegin(), query.rend());
                    if (complete){
                        isFound = isFound && (FMIndex->*extraChar)(FMIndex->char2comp[','], res.s) > 0;
                    } 
                    if (isFound){ //get exact and if found push to Sols{}
                        res.k = minimalEDOfRow;                
                        BuildSolution(res, row, j);                        
                        
                    }
                    
                }

            }
        }
        else {
            if (M(row, qSize) <= maxED){
                //search separator for complete trace
                searchstate spp(sp);
                //true for all partial
                bool found = complete?(FMIndex->*extraChar)(FMIndex->char2comp[','], spp) > 0:true;  
                if (found){
                    if(M(row, qSize) < maxED){ // updated maxED
                        maxED = M(row, qSize);
                        sols->clear();
                    }  
                    BuildSolution(spp, row, qSize);
                }
                            
            }
            if (row < qSize + maxED) 
                pushChildren(sp, row);
        }

        
    }
}


template<typename It>
bool ED:: exactMatching(Result &res, int offset, It st, It ed){
    for(auto it = st + offset; it!=ed; it++){
        size_t result = (FMIndex->*extraChar)(FMIndex->char2comp[*it], res.s);
        if (result < 1) return false;
        res.alignment->push_back('-');
    }
    return true;
}

bool ED::matchExactly(const std::string &trace){
    searchstate s;
    FMIndex->init_search_state(s);
    Result res(s);
    bool isFound = fw ? exactMatching(res,0,trace.begin(), trace.end()) : exactMatching(res,0,trace.rbegin(), trace.rend());
    if (isFound){
        sols->push_back(std::move(res));
    }
    return isFound;
    
}
