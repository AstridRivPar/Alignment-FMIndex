
#include <iostream>
#include "Bi_FM_Index.hpp"
#include <Poco/Logger.h>
#include "Helpers.hpp"
#include "Matrix.hpp"
#include "Reader.hpp"
class ED{
    static Poco::Logger &log;


    private:
       
    

        Bi_FM_Index *FMIndex;
        const std::string *text;
        bool (Bi_FM_Index::*extraChar)(SARangePair&, uint8_t);
        bool fw = true;
        std::vector<Node> nodesToCheck;
        Matrix M;
        std::vector<Result> *sols;
        std::vector<int> foundPerK;

        

    public:
        std::map<std::string, std::vector<Result>> results;
    
        ED():FMIndex(NULL){};
        /**
         * Constructor 
         * @param index FM-Index
         * @param txt Text
         */
        ED(Bi_FM_Index *index, const std::string *text):FMIndex(index), text(text){};
        /**
         * Match exactly a trace or a part of the trace
         * @param res stores the alignment 
         * @param offset for a partial alignment
         * @param st Iterator or reverse iterator (start)
         * @param ed ITerator or reverse iterator (end)
         */

        
        template <typename It>
        int exactMatching(Result &res, int offset, It st, It ed);

        /**
        * Match a trace with maxED = 0
        */
        bool matchExactly(const std::string &query);

        int CheckConf(const std::string &trace, int maxED, bool complete);

        void CompleteAligment(const std::vector<std::string> &queries, int maxED, bool isFw, bool complete);

        void setDirection( bool fw);
        void BFSearch(const std::string &query, int maxED, bool complete, const SARangePair &s);
        /**
         * Returns the optimal value
         */
        int TakeAllOptimal();
        void BuildSolution(Result &res, int row, int col, const std::string &query, bool complete, int count);
        void BuildSolution(SARangePair &range, int row, int col, const std::string &query, bool complete);

        void pushChildren(const SARangePair &s, int row);
        /**
         * Get number of optimal traces for a given trace
         * @param query Query string for which we want to return the number of optimal traces
         */

        std::vector<int> getFoundPerK();

        const std::map<std::string, std::vector<Result>> & get_results()  const  { return results;}
        void printAll(std::ofstream& of, const Reader&ccf);
// }
};
