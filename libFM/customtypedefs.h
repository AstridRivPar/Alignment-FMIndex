#include <vector>
#include <string>
#include "indexes.hpp"

typedef index_bidirectional_waveletindex<> Index; 
typedef index_bidirectional_waveletindex<>::searchstate searchstate;
    
struct Node{
    // typedef index_bidirectional_waveletindex<> Index; 
    // typedef index_bidirectional_waveletindex<>::searchstate searchstate;
    Node(const searchstate &ns, char cc, const int &r): 
        s(ns),
        row(r),
        c(cc){}

    searchstate s;
    int row = 0;
    char c;
};
class Result{
    public:
        searchstate s;
        std::unique_ptr<std::string> alignment;
        std::unique_ptr<std::string> match;
        int k = 0;
        // Result (Result res){

        // }
        Result(const searchstate &s, std::unique_ptr<std::string> &alignment, int k): s(s), k(k),
        alignment(std::move(alignment)){}

        Result(std::unique_ptr<std::string> &alignment, std::unique_ptr<std::string> &match,int k): alignment(std::move(alignment)), 
        k(k), match(std::move(match)){}

        Result(const searchstate &ns):s(ns){
            alignment = std::unique_ptr<std::string>(new std::string());
        }

        std::tuple<int, int, int> countOperations() const{
            int d = 0, i = 0, r = 0;
            for (const auto &c: *this->alignment){
                switch (c){
                    case 'd':
                        d++;
                        break;
                    case 'i':
                        i++;
                        break;
                    case 'r':
                        r++;
                        break;
                    default:
                        break;
                    }

                }
            return std::make_tuple(d,i,r);
    }
};