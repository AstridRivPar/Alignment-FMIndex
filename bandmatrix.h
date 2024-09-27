

#ifndef BANDMATRIX_H
#define BANDMATRIX_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <deque>
#include <optional>
#include <memory>

typedef uint32_t length_t;

// ============================================================================
// CLASS EDIT DISTANCE
// ============================================================================

// The edit Distance class has two fields, an unsigned integer that stores the
// actual edit distance and a boolean flagged, signifiying if this ED originated
// from a flagged path

// the data is a 32 unsigned bit and the most significant bit is reserved for
// the flag while the other 31 bits are for the actual edit distance
class EditDistance {
  private:
    length_t data;

  public:
    EditDistance() {
        data = 0u;
    }

    EditDistance(int edit) {
        data = edit;
    }

    EditDistance(int edit, bool flag) {
        data = (flag) ? (1u << 31) | edit : edit;
    }

    void setED(int edit) {
        data = (data & 31) | edit;
    }

    unsigned int getED() const {
        unsigned int ret = (data << 1) >> 1;
        return ret;
    }

    bool isFlagged() const {
        return ((data >> 31) > 0);
    }

    /**
     * Operator overloading creates a new edit distance with the added distance
     * @param toAdd the value to add to this edit distance
     * @returns a new EditDistance instance with the increased edit distance and
     * same flag as the current one
     */
    EditDistance operator+(unsigned int toAdd) const {

        return EditDistance(getED() + toAdd, isFlagged());
    }

    /**
     * Operator overloading. Comparing two edit distances goes as follows: if
     * their ED's are equal then a flagged instance is smaller than a non
     * flagged. Else the one with the smallest ED is the smallest
     * @param argument the EditDistance to compare to this
     * @returns a bool indicating wheter this is smaller
     */
    bool operator<=(const EditDistance& argument) const {
        // first, look at equal edit distance

        if (getED() == argument.getED()) {
            // both flagged or not flagged => equal => true
            // if this flagged and other not => true
            // if this is not flagged and other is => false
            return (isFlagged() || !argument.isFlagged());
        }
        return getED() < argument.getED();
    }

    /**
     * Operator overloading.
     * This returns a bool indicating wheter this is smaller then some integer
     * @param argument the integer to compare to
     * @returns true if this is not flagged and smaller or equal, else it will
     * return false
     */
    bool operator<=(unsigned int argument) const {
        // a flagged reported instance is never smaller then an argument
        return !isFlagged() && getED() <= argument;
    }

    /**
     * Makes a string representation of this.
     * @returns the string representation of this EditDistance
     */
    std::string to_string() {
        std::string flagString = (isFlagged() ? "*" : "");
        return std::to_string(getED()) + flagString;
    }
};
class Matrix{
    private:
        std::vector<int> matrix;
        
    public:
        int cols;
        int rows;
        int cost = 1; //without replace
        
        Matrix(){}
          /**
         * Constructor
         * @param m Number of rows
         * @param W Number of mismatches allowed
         * @param startValue minED for the row
         */
        Matrix(int l, int W, int startValue) : cols(l + 1), rows(l + 1 + W) {
            matrix.resize(cols * rows);
            initializeMatrix(startValue);
        }

        /**
         * Constructor without replace
         * @param m Number of rows
         * @param W number of mismatches allowed
         */
        Matrix(int l, int W) : cols(l + 1), rows(l + 1 + W) { 
            matrix.resize(cols * rows);
        }

        /**
         * Initializes the matrix and sets the number of cols and rows
         * @param m Number of rows
         * @param W number of mismatches allowed
         */
        void init(int l, int W){
            cols = l + 1;
            rows = l + 1 + W; 
            matrix.resize(cols*rows);
            initializeMatrix(0);
        }
        /**
         * updateMatrix without replace
         * @param notMatch Elements not equal = true
         * @param row 
         * @param column
         */

        void initializeMatrix(int startValue) {
            for(int i = 0; i < cols; i++){
                matrix[i] = i + startValue;
            }
            for (int i = 1; i < rows; i++){
                matrix[ i * cols] = i + startValue;
        
            }
        }

        int& operator()(int i, int j) {
            return matrix[i * cols + j];
        }
        void printMatrix() {
            for(int i = 0; i  < rows; i ++){
                std::string row = "";
                for (int j = 0; j < cols; j ++){
                    row += std::to_string(matrix[i * cols + j]) + "\t";
                }
                std::cout << row << std::endl;
            }
        }
        void updateMatrix(bool notMatch, int row, int column) {
            int diag = operator()(row - 1, column - 1) + cost*(notMatch);
            int left = operator()(row, column - 1) + 1;
            int up = operator()(row - 1, column) + 1;
            operator()(row, column) = std::min<int>(diag, std::min<int>(left, up));
            //std::min<int>({diag, left, up});
        }


        int getRows(){
            return rows;
        }
        /**
         * without replace
        */
        void getAlignment(int i, int j, std::unique_ptr<std::string> &alignment){
            while (!(i == 0 && j == 0)){
                if (i > 0 && operator()(i, j) == (operator()(i - 1, j) + 1))
                {
                    alignment->push_back('i'); //vertical gap move in the model
                    i--;
                } else if (j > 0 && operator()(i, j) == (operator()(i, j - 1) + 1) ){ //horizontal gap move in the log
                    alignment->push_back('d');
                    j--;
                }else{
                    alignment->push_back('-');                
                    j--; i--;
                }
                
                
            }
            
        }

        /**
         * With replace
        */
        void getAlignmentWR(int i, int j, std::unique_ptr<std::string> &alignment){

            while (!(i == 0 && j == 0)){
                if (i > 0 && operator()(i, j) == (operator()(i - 1, j) + 1))
                {
                    alignment->push_back('d');
                    i--;
                } else if (j > 0 && operator()(i, j) == (operator()(i, j - 1) + 1) ){
                    alignment->push_back('i');
                    j--;
                }else{
                    if (operator()(i, j) == (operator()(i - 1, j - 1) + 1)){
                        alignment->push_back('r');
                    }else{
                        alignment->push_back('-');
                    }
                        
                    j--; i--;
                }
                
                
            }
        }

        void (Matrix::*getAlign(bool replace))(int, int, std::unique_ptr<std::string>&){
            if (replace) {
                return &Matrix::getAlignmentWR;
            }
            this->cost = 2;
            return &Matrix::getAlignment;
            
        }

    

};

// ============================================================================
// CLASS BANDED MATRIX
// ============================================================================

// The band matrix class can best be understood as m horizontal bands each
// with a width of 2W+1 elements. It is allocated as a single m*(2W+1) array.
// For example, for W = 2 and m = 6.
// XX|XXX...|
//  X|XXXX..|
//   |XXXXX.|
//   |.XXXXX|
//   |..XXXX|X
//   |...XXX|XX
// The actual matrix is between |.| The Xs left and right of the |.| are
// allocated but should in principle not be addressed.
// The storage order is row-major.

class BandMatrix {
  private:
    std::vector<int> matrix;
    int W;
    int m;
    int cols;
    int cost = 1;

  public:
    /**
     * Constructor
     * @param m Number of rows
     * @param W Number of off-diagonal elements (one sided)
     */
    BandMatrix(int m, int W, int startValue) : W(W), m(m) {
        matrix.resize(m * (2 * W + 1));
        initializeMatrix(startValue);
    }
    /**
     * Constructor
     * @param m Number of rows
     * @param W Number of off-diagonal elements (one sided)
     */
    BandMatrix(int m, int W) : W(W), m(m) {
        
        cols = 2 * W + 1;

        matrix.resize(m * (cols), W + 1);
        initializeMatrix(0);
    }


    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Reference to element at position (i, j)
     */
    int& operator()(int i, int j) {
        return matrix[i * (cols) + j - i + W];
    }

    /**
     * Get the band width
     * @return The band width
     */
    int getWidth() const {
        return W;
    }
    int getCols() const {
        return cols;
    }
    int getRows() const {
        return m;
    }

    void printMatrix() {
        for (int i = 0; i < m; i++) {
            int firstCol = std::max<int>(1, i - W);
            int lastCol = std::min<int>(m, i + W);
            std::string row = "";
            for (int j = 0; j < firstCol; j++) {
                row += ".\t";
            }
            for (int j = firstCol; j <= lastCol; j++) {
                int number = operator()(i, j);
                row += std::to_string(number) + "\t";
            }
            for (int j = lastCol + 1; j < 2 * W + 1; j++) {
                row += ".\t";
            }
            std::cout << row << std::endl;
        }
    }

    void initializeMatrix(int startValue) {
        for (int i = 0; i <= W; i++) {
            matrix[i * cols - i + W] = i + startValue;
            matrix[i + W] = i + startValue;
        }
    }

    /**
     * Update the matrix by calculating the element at position row, column.
     * @param match whether the character was a match
     * @param row the row of the elment to update
     * @param collumn the collumn of the element to update
     */
    void updateMatrix(bool notMatch, unsigned int row, unsigned int column) {
        int diag = operator()(row - 1, column - 1) + (cost * notMatch);

        // watch out for integer overflow
        unsigned int diffRowAndW = (row >= W) ? row - W : 0;
        int gapX =
            (column > diffRowAndW) ? operator()(row, column - 1) + 1 : diag + 1;
        int gapY =
            (column < row + W) ? operator()(row - 1, column) + 1 : diag + 1;

        operator()(row, column) =
            std::min<int>(diag, std::min<int>(gapX, gapY));
    }

    void getAlignment(int i, int j, std::unique_ptr<std::string> &alignment){
        while (!(i == 0 && j == 0)){
            if (i > 0 && operator()(i, j) == (operator()(i - 1, j) + 1))
            {
                alignment->push_back('d');
                i--;
            } else if (j > 0 && operator()(i, j) == (operator()(i, j - 1) + 1) ){
                alignment->push_back('i');
                j--;
            }else{
                alignment->push_back('-');                
                j--; i--;
            }
            
            
        }
            
    }

    /**
     * With replace
    */
    void getAlignmentWR(int i, int j, std::unique_ptr<std::string> &alignment){

        while (!(i == 0 && j == 0)){
            if (i > 0 && operator()(i, j) == (operator()(i - 1, j) + 1))
            {
                alignment->push_back('d');
                i--;
            } else if (j > 0 && operator()(i, j) == (operator()(i, j - 1) + 1) ){
                alignment->push_back('i');
                j--;
            }else{
                if ( operator()(i, j) == (operator()(i - 1, j - 1) + 1)){
                    alignment->push_back('r');
                }else{
                    alignment->push_back('-');
                }
                    
                j--; i--;
            }
            
            
        }
    }

    void (BandMatrix::*getAlign(bool replace))(int, int, std::unique_ptr<std::string>&){
        if (replace) {
            return &BandMatrix::getAlignmentWR;
        }
        this->cost = 2;
        return &BandMatrix::getAlignment;
        
    }
};

// ============================================================================
// CLASS EDIT MATRIX
// ============================================================================

// This is a band matrix but it elements are of class
// EditDistance. This checks whether an EditDistance originates from a path that
// has a leading gap and will flag as such

class EditMatrix {
  private:
    std::vector<EditDistance> matrix; // the matrix
    length_t W;                       // the width of the matrix

    /** Helper function for update matrix, it checks which path to take, the
     * diagonal, horizontal or verical path.
     * @param diag the value for the element if the diagonal path is chosen
     * @param gapX the value for the element if the horizontal path is chosen
     * @param gapY the value for the elemnt if the vertical path is chosen
     *
     * @returns the best element out of the three provided.
     */
    EditDistance chooseBestElement(EditDistance& diag, EditDistance& gapX,
                                   EditDistance& gapY) {
        // need to check if for the optimal value one of the paths is flagged
        // so we need to find the smallest

        if ((diag <= gapX) && (diag <= gapY)) {
            return diag;
        }
        if (gapX <= gapY) {
            return gapX;
        }
        return gapY;
    }

    /**
     * Helper function for initializing the matrix. The top row and left most
     * collumn are filled in with as value their row/column number plus a
     * startvalue. If the noLeadingGaps argument is true than the leftmost
     * column will be flagged (except the origin).
     * @param noLeadingGaps a bool to indicate if leading gaps are allowed
     * @param startValue the value of the origin, default is 0
     */
    void initializeMatrix(bool noLeadingGapsText, bool noLeadingGapsPattern,
                          int startValue = 0) {
        for (length_t i = 0; i <= W; i++) {
            operator()(i, 0) = EditDistance(i + startValue, noLeadingGapsText);
            operator()(0, i) =
                EditDistance(i + startValue, noLeadingGapsPattern);
        }
        operator()(0, 0) = EditDistance(startValue, false);
    }

  public:
    /**
     * Constructor
     * @param m Number of rows
     * @param W Number of off-diagonal elements (one sided)
     */
    EditMatrix(length_t m, int W, int startValue, bool leadingGapsAllowedText,
               bool leadingGapsAllowedPattern)
        : W(W) {
        matrix.resize(m * (2 * W + 1));
        initializeMatrix(!leadingGapsAllowedText, !leadingGapsAllowedPattern,
                         startValue);
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Element at position (i, j)
     */
    EditDistance operator()(length_t i, int j) const {
        return matrix[i * (2 * W + 1) + j - i + W];
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Reference to element at position (i, j)
     */
    EditDistance& operator()(length_t i, int j) {
        return matrix[i * (2 * W + 1) + j - i + W];
    }

    /**
     * Get the band width
     * @return The band width
     */
    int getWidth() const {
        return W;
    }

    /**
     * Prints the matrix, for debugging purposes
     */
    void printMatrix() {
        length_t m = matrix.size() / (2 * W + 1);

        for (length_t i = 0; i < m; i++) {
            length_t firstCol = std::max<int>(1, i - W);
            length_t lastCol = std::min<int>(m, i + W);
            std::string row = "";
            for (length_t j = 0; j < firstCol; j++) {
                row += ".\t";
            }
            for (length_t j = firstCol; j <= lastCol; j++) {
                row += (operator()(i, j)).to_string() + "\t";
            }
            for (length_t j = lastCol + 1; j < 2 * W + 1; j++) {
                row += ".\t";
            }
            std::cout << row << std::endl;
        }
    }

    void printRow(int i) {
        length_t m = matrix.size() / (2 * W + 1);
        length_t firstCol = std::max<int>(1, i - W);
        length_t lastCol = std::min<int>(m, i + W);
        std::string row = "";
        for (length_t j = 0; j < firstCol; j++) {
            row += ".\t";
        }
        for (length_t j = firstCol; j <= lastCol; j++) {
            row += (operator()(i, j)).to_string() + "\t";
        }
        for (length_t j = lastCol + 1; j < 2 * W + 1; j++) {
            row += ".\t";
        }
        std::cout << row << std::endl;
    }

    /**
     * Update the matrix by calculating the element at position row, column.
     * @param match whether the character was a match
     * @param row the row of the elment to update
     * @param collumn the collumn of the element to update
     */
    void updateMatrix(bool notMatch, unsigned int row, unsigned int collumn) {
        EditDistance diag = operator()(row - 1, collumn - 1) + notMatch;

        // watch out for integer overflow
        unsigned int diffRowAndW = (row >= W) ? row - W : 0;
        EditDistance gapX =
            (collumn > diffRowAndW) ? operator()(row, collumn - 1) + 1
                                    : diag + 1;
        EditDistance gapY =
            (collumn < row + W) ? operator()(row - 1, collumn) + 1 : diag + 1;

        operator()(row, collumn) = chooseBestElement(diag, gapX, gapY);
    }
};

#endif
