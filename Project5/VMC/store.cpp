#include "store.h"


void store::bin(vec Array)
{
    int ArrayDim = Array.size();
    double * tempArray = new double[ArrayDim];
    for(int i = 0; i < ArrayDim; i++)
    {
        tempArray[i] = Array[i];
    }
    std::ofstream file
        (
            m_filename, std::ofstream::binary
        );
    file.write
        (
            reinterpret_cast<const char*> (tempArray),
            ArrayDim*sizeof(double)
        );
    file.close();
    delete [] tempArray;
} // end of bin
