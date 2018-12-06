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
        ~store(); // destructor

        store(string filename)
        { // constructor w/ filename
            m_filename = filename;
        }

       // setters and getters?
       void setFILE(string filename)
       {
           m_filename = filename;
       }

        // open file
        void open()
        {
            m_file.open(
                    m_filename, std::ofstream::out
                );
        }

        // writes
        void dat()
        {
            m_file << m_line << "\n";
        }

        void bin(vec Array);

        // Lines
        void lineadd(string element)
        {
            m_line += element;
        }

        // close file
        void close()
        {
            m_file.close();
        }
//
    private:
        std::ofstream m_file;
        string m_filename, m_line;
};

#endif // end of store header
