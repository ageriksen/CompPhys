#ifndef STORE_H
#define STORE_H

#include <fstream>
#include <string>
#include <armadillo>

using std::string;
using arma::mat;
using arma::vec;

class store
{
    public:
        store();  // construtor

        store(string fileName)
        { // constructor w/ fileName
            m_fileName = fileName;
        }

       // setters and getters?
       void setFile(string fileName)
       {
           m_fileName = fileName;
       }
       string getLine() { return m_line; }

        // open file
        void open()
        {
            m_file.open(
                    m_fileName, std::ofstream::out
                );
        }

        // writes
        void dat()
        {
            m_file << m_line << "\n";
        }
        void dat(string line)
        {
            m_file << line << "\n";
        }

        void bin(vec Array);

        // Lines
        void lineAdd(string element)
        {
            m_line += element+" ";
        }
        void lineClean() { m_line = ""; }

        // close file
        void close()
        {
            m_file.close();
        }
//
    private:
        std::ofstream m_file;
        string m_fileName;
        string m_line;
};

#endif // end of store header
