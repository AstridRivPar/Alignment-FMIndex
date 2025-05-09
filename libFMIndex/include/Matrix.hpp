
#ifndef BANDMATRIX_H
#define BANDMATRIX_H


#include <iostream>

#include <Poco/Logger.h>

class Matrix{
    static Poco::Logger &log;

    private:
        std::vector<int> matrix;
        
    public:
        int cols;
        int rows;
        int cost = 2; //without replace
        
        Matrix();
          /**
         * @param m Number of rows
         * @param W Number of mismatches allowed
         * @param startValue minED for the row
         */
        Matrix(int l, int W, int startValue);

        /**
         * @param m Number of rows
         * @param W number of mismatches allowed
         */
        Matrix(int l, int W);

        /**
         * Initializes the matrix and sets the number of cols and rows
         * @param m Number of rows
         * @param W number of mismatches allowed
         */
        void init(int l, int W);
        /**
         * updateMatrix without replace
         * @param notMatch Elements not equal = true
         * @param row 
         * @param column
         */

        void initializeMatrix(int startValue);

        int& operator()(int i, int j) {
            return matrix[i * cols + j];
        }

        void printMatrix(int row);

        void printMatrix();

        void updateMatrix(bool notMatch, int row, int column);


        int getRows();

        // void getAlignment(int i, int j, std::unique_ptr<std::string> &alignment, std::unique_ptr<std::string> & match,  const std::string & query);
        
        template <typename Iterator>
        void getAlignment(int i, int j, Iterator q_beg, Iterator m_beg, std::string & alignment){
            
            // std::cout << "Q: "<< query << " M: "<< match_ << " " << i << " " << j << std::endl;
            // printMatrix(match->size());
            
            
            // if (query == "bgfhijkstACEFHPOIJKLDMQRNuvuvuwxyz" && match_ == "acdeijkjkjkltAECFHPOIJKLDMQRNuwxyz"){
                // printMatrix();

            // }
            int counter = 0, i_ = 0, j_= 0;
            while (!(i == 0 && j == 0)){
                
                
                int diag = operator()(i - 1, j - 1); // b + cost
                int left = operator()(i, j - 1); // c + 1
                int up = operator()(i - 1, j); // a + 1
                int r = operator()(i, j);
        
                // std::cout << i << " " << j << " " <<i_ << " "<< j_<< " " << *m_beg << " "<< *q_beg << std::endl;
        
                // std::cout << diag << " " << left << " " << up << " " <<r  << std::endl;
                i_ = i; j_ = j; 
                if (*m_beg == *q_beg && diag == r){
                    alignment.push_back('-');
                                   
                    j--; i--;
                    q_beg++;m_beg++;
                }
                else if (i> 0 && up + 1 == r){
                    alignment.push_back('i'); //vertical gap move in the model
                    i--;
                    m_beg++;
                }
                else if (j > 0 && left + 1 == r){
                    alignment.push_back('d'); //vertical gap move in the model
                    j--;
                    q_beg++;
                }
        
                if (i_ == i && j_ == j){
                    log.error("ERROR IN ALIGNMENT!", __FILE__, __LINE__);
                    std::abort();
                }
                
            }
            
        }
        

};

#endif