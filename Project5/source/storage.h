#ifndef STORAGE_G
#define STORAGE_G

#include <fstream>
#include <string>
#include <armadillo>
#include <iostream>

using std::string;
using arma::mat;
using arma::vec;
using std::cout;
using std::endl;

class storage
{
    public:
        storage() {}  // construtor

        storage(string fileName)
        { // constructor w/ fileName
            m_fileName = fileName;
        }


       //-----------------------------------------------
       // writes
       //-----------------------------------------------
       // open and closefile
       void out()
       {
           m_oFile.open(
                   m_fileName, std::ofstream::out
               );
           if( !m_oFile )
           {
                std::cout << "couldn't open file" << m_fileName << std::endl;
           }
       }
       void dat()
       {
           m_oFile << m_line << "\n";
       }
       void dat(string line)
       {
           m_oFile << line << "\n";
       }
       void bin(vec Array);
       //-----------------------------------------------
       // setters and getters?
       //-----------------------------------------------
       void name(string fileName)
       {
           m_fileName = fileName;
       }
       string getLine() { return m_line; }
       string getName() { return m_fileName; }



       //-----------------------------------------------
       // Reads
       //-----------------------------------------------
       void in( std::vector<double> *importVec );

       //-----------------------------------------------
       // Lines
       //-----------------------------------------------
       void lineAdd(string element)
       {
           m_line += element+" ";
       }
       void lineClean() { m_line = ""; }

       //-----------------------------------------------
       // Close
       //-----------------------------------------------
       void close()
       {
           m_oFile.close();
       }
//
    private:
        std::ofstream m_oFile;
        std::ifstream m_iFile;
        string m_fileName;
        string m_line;
};

#endif // end STORAGE_G
