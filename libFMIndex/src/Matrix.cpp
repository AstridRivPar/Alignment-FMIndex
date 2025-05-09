#include <Matrix.hpp>

Poco::Logger &Matrix::log = Poco::Logger::get("Matrix");

Matrix::Matrix(){}

Matrix::Matrix(int l, int W, int startValue) : cols(l + 1), rows(l + 1 + W) {
    matrix.resize(cols * rows);
    initializeMatrix(startValue);
}

Matrix::Matrix(int l, int W) : cols(l + 1), rows(l + 1 + W) { 
    matrix.resize(cols * rows);
}

void Matrix::init(int l, int W){
    cols = l + 1;
    rows = l + 1 + W; 
    matrix.resize(cols*rows);
    initializeMatrix(0);
}

void Matrix::initializeMatrix(int startValue) {
    for(int i = 0; i < cols; i++){
        matrix[i] = i + startValue;
    }
    for (int i = 1; i < rows; i++){
        matrix[ i * cols] = i + startValue;

    }
}

void Matrix::printMatrix(int row){
    for(int i = 0; i  < row; i ++){
        std::string row = "";
        for (int j = 0; j < cols; j ++){
            row += std::to_string(matrix[i * cols + j]) + "\t";
        }
        std::cout << row << std::endl;
    }
}

void Matrix::printMatrix(){
    for(int i = 0; i  < rows; i ++){
        std::string row = "";
        for (int j = 0; j < cols; j ++){
            row += std::to_string(matrix[i * cols + j]) + "\t";
        }
        std::cout << row << std::endl;
    }
}
void Matrix::updateMatrix(bool notMatch, int row, int column) {
    int diag = operator()(row - 1, column - 1) + cost*(notMatch);
    int left = operator()(row, column - 1) + 1;
    int up = operator()(row - 1, column) + 1;
    operator()(row, column) = std::min<int>(diag, std::min<int>(left, up));
}


int Matrix::getRows(){
    return rows;
}



// void Matrix::getAlignment(int i, int j, std::unique_ptr<std::string> &alignment, std::unique_ptr<std::string> & match,  const std::string & query){
//     const std::string & match_ = *match;
//     std::cout << "Q: "<< query << " M: "<< match_ << " " << i << " " << j << std::endl;
//     printMatrix(match->size());
    
    
//     // if (query == "bgfhijkstACEFHPOIJKLDMQRNuvuvuwxyz" && match_ == "acdeijkjkjkltAECFHPOIJKLDMQRNuwxyz"){
//     //     printMatrix();
//     // }
//     int counter = 0, i_ = 0, j_= 0;
//     while (!(i == 0 && j == 0)){
        
        
//         int diag = operator()(i - 1, j - 1); // b + cost
//         int left = operator()(i, j - 1); // c + 1
//         int up = operator()(i - 1, j); // a + 1
//         int r = operator()(i, j);

//         std::cout << i << " " << j << " " <<i_ << " "<< j_<< " " << match_[i - 1] << " "<< query[j - 1] << std::endl;

//         std::cout << diag << " " << left << " " << up << " " <<r  << std::endl;
//         i_ = i; j_ = j; 
//         if (match_[i - 1]== query[j - 1] && diag == r){
//             alignment->push_back('-');
                           
//             j--; i--;
//         }
//         else if (i> 0 && up + 1 == r){
//             alignment->push_back('i'); //vertical gap move in the model
//             i--;
//         }
//         else if (j > 0 && left + 1 == r){
//             alignment->push_back('d'); //vertical gap move in the model
//             j--;
//         }

//         if (i_ == i && j_ == j){
//             log.error("ERROR IN ALIGNMENT!", __FILE__, __LINE__);
//             std::abort();
//         }
        
//     }
    
// }

