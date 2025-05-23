#include "ED.hpp"
Poco::Logger &ED::log = Poco::Logger::get("ED");

void ED::setDirection(bool isFw){
    fw = isFw;
    extraChar = fw ? &Bi_FM_Index::forward_search : &Bi_FM_Index::backward_search;
}

std::vector<int> ED::getFoundPerK(){
    return foundPerK;
}


void ED::pushChildren(const SARangePair &s, int row){
    uint8_t alphabet_length = FMIndex->get_alphabet_size(); 
    // std::cout << "Starting" << std::endl;
    for (uint8_t c = alphabet_length - 1; c >= 2; c--){

        SARangePair rp(s);
        char cc = FMIndex->comp2char(c);
        // std::cout << "real: " << cc << std::endl;
        if ((FMIndex->*extraChar)(rp, c)){
            // rp.pMatch += c;
            nodesToCheck.emplace_back(rp, cc, row + 1);
        }
    }
}

int ED::TakeAllOptimal(){
    if (sols->empty()) return -1;
    
    return sols->back().k;    
}

int ED::CheckConf(const std::string &query, int maxED, bool complete){
    //init
    nodesToCheck.clear();
    nodesToCheck.reserve(100);
    
    log.trace("Initiating alignment", __FILE__, __LINE__);
    //if trace is already aligned, meaning there are repeated traces      
    
    results[query] = std::vector<Result>();
    sols = &results[query];
    //match exactly
    if (maxED == 0) {
        std::string al = "";
        bool isFound = matchExactly(query, al);
        return isFound? 0: -1;
    }
    //initialize matrix
    log.trace("Init matrix", __FILE__, __LINE__);
    M.init(query.size(), maxED);
    
    
    SARangePair init(FMIndex->size());
    

    if(complete){
        (FMIndex->*extraChar)(init, FMIndex->char2comp(','));
    }
    BFSearch(query, maxED, complete, init);
    return TakeAllOptimal();    
}
void ED::CompleteAligment(const std::vector<std::string> &queries, int maxED, bool isFw, bool complete){

    setDirection(isFw);
    foundPerK.resize(maxED + 2);
    int counter = 0;
    for (const auto & query: queries){
        if (counter % 20 == 0){
            log.information("Query " + std::to_string(counter) + ": "+ query, __FILE__, __LINE__);
        }
        
        int pK = -1;
        if (results.find(query) != results.end()){
            sols = &results[query];
            if(!sols->empty())
                pK = sols->back().k;
        }else{
            
            pK = CheckConf(query, maxED, complete);
        }
        
        foundPerK[pK != -1 ? pK : maxED + 1]++;
        counter++;
    }
    

    
}

void ED::BuildSolution(SARangePair & range, int row, int col, const std::string &query, bool complete){
    std::string al;
    //extract only first occurrence
    auto firstOcc = FMIndex ->get_occ(range);
    int offset = complete ? 1 : 0;
    
    const std::string m(text->substr(firstOcc + offset, row));
    // std::cout << "App " << query << " Match " << m << std::endl;
    std::vector<std::string> alignments = fw ? M.getAlignment(row,col, query.rbegin(), m.rbegin(), al, fw)
        : M.getAlignment(row,col, query.begin(), m.begin(), al, fw);

    // if (fw) std::reverse(al.begin(), al.end());
    sols->emplace_back(range, std::move(alignments), std::move(m), M(row, col));
    
}

void ED::BuildSolution(Result &res, int row, int col, const std::string &query, bool complete, int count, std::string &al){
   
    auto firstOcc = FMIndex ->get_occ(res.range);
    int init_offset = complete? 1 : 0; //including first separator
    int final_offset = complete? 2 : 0;

    res.match = text->substr(firstOcc + init_offset, res.range.get_length() - final_offset);
    const std::string &m = res.match;
    // std::cout << "Exact " << query << " Match " << m << std::endl;
    // std::vector<std::string> getAlignment(int row, int col, Iterator q_beg, Iterator m_beg, std::string &al)
    std::vector<std::string> alignments = fw ? M.getAlignment(row, col, query.rbegin() + count, m.rbegin() + count, al, fw)
        : M.getAlignment(row,col,query.begin() + count, m.begin() + count, al, fw);
    
    res.alignments = alignments;
    //Add exact part
    // if (fw) std::reverse(res.alignment.begin(), res.alignment.end());
    sols->push_back(std::move(res));
}
void ED::BFSearch(const std::string &query, int maxED, bool complete, const SARangePair &s){
    
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
                    std::string al = "";
                    int count = fw? exactMatching(res, j, query.begin(), query.end(), al):exactMatching(res, j, query.rbegin(), query.rend(), al);
                    
                    bool isFound = count > -1;

                    if (complete){
                        isFound = isFound && (FMIndex->*extraChar)(res.range, FMIndex->char2comp(','));

                    } 
                    if (isFound){ //get exact and if found push to Sols{}
                        res.k = minimalEDOfRow;  
                        BuildSolution(res, row, j, query, complete, count, al);                        
                        
                    }
                    
                }

            }
        }
        else {
            if (M(row, qSize) <= maxED){
                SARangePair spp(sp);
                //true for all partial
                //search separator for complete trace
                 
                bool found = complete?(FMIndex->*extraChar)(spp, FMIndex->char2comp(',')):true;  
                if (found){
                    if(M(row, qSize) < maxED){ // updated maxED
                        maxED = M(row, qSize);
                        sols->clear();
                    }  
                    BuildSolution(spp, row, qSize, query, complete);
                }
                            
            }
            if (row < qSize + maxED) 
                pushChildren(sp, row);
        }

        
    }
}


template<typename It>
int ED:: exactMatching(Result &res, int offset, It st, It ed, std::string &al){
    int count = 0;
    for(auto it = st + offset; it!=ed; it++){
        bool result = (FMIndex->*extraChar)(res.range, *it);
        if (!result) return -1;
        // res.range.pMatch += *it;
        al.push_back('-');
        count ++;
    }
    return count;
}

bool ED::matchExactly(const std::string &trace, std::string &al){
    SARangePair init_range(FMIndex->size());
    Result res(init_range);
    int isFound = fw ? exactMatching(res,0,trace.begin(), trace.end(), al) : exactMatching(res,0,trace.rbegin(), trace.rend(), al);
    if (isFound >= 0){
        res.alignments.push_back(al);
        sols->push_back(std::move(res));
    }
    return isFound;
    
}

void ED::printAll(std::ofstream& of, const Reader&ccf){
    of << "Trace,NoTraces,K,NoOptAl,Alignments" << std::endl;
    const auto & queries = ccf.get_log();
    const auto & counter_vec = ccf.get_no_traces();
    for (size_t q_idx = 0; q_idx < queries.size(); q_idx++){
        const auto & q = queries[q_idx];
        auto &res = results.at(q);


        if(!res.empty()){
            of << q << ","<<counter_vec[q_idx]<<","<< res[0].k << "," << res.size() << ",";

                        
            for (size_t i = 0; i < res.size(); i++){
                const Result &element = res[i];
                of <<"{"<<element<< "};";

            }
            
           
        }else{
            of << q << ",No matches";

        }
        of << std::endl;
    }
}

