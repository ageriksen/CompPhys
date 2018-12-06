#ifndef RANDOM_H
#define RANDOM_H

#include <random>

class random
{
    public:

        random()
        {
            std::random_device rd;
            m_engine = std::mt19937_64(rd());
        }

        random( int seed ) : m_seed(seed)
        {
            m_engine = std::mt19937_64(m_seed);
        }


        void UniformDistribution()
        {
            m_uniform_real = std::uniform_real_distribution<double>(0.0, 1.0);
        }

        double getUniformReal()
        {
            return m_uniform_real(m_engine);
        }


    private:
        std::mt19937_64 m_engine;
        std::uniform_real_distribution<double> m_uniform_real;

        int m_seed = 0;


};
#endif // end of header
