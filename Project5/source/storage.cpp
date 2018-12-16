#include "storage.h"
#include <armadillo>
#include <string>
#include <sstream>
#include <vector>


void storage::bin(vec Array)
{
    int ArrayDim = Array.size();
    double * tempArray = new double[ArrayDim];
    for(int i = 0; i < ArrayDim; i++)
    {
        tempArray[i] = Array[i];
    }
    std::ofstream file
        (
            m_fileName, std::ofstream::binary
        );
    file.write
        (
            reinterpret_cast<const char*> (tempArray),
            ArrayDim*sizeof(double)
        );
    file.close();
    delete [] tempArray;
} // end of bin


void storage::in( std::vector<double> *importVec )
{
    std::cout << "opening file" << std::endl;
    m_iFile.open(m_fileName);
    std::string str;
    while( std::getline( m_iFile,  str ) )
    {
        if( sizeof(str) > 0 )
        {
            std::istringstream strstream(str);
            double strDouble;
            strstream >> strDouble;
            importVec -> push_back(strDouble);
        }
    }
    m_iFile.close();
} // end if inVector
