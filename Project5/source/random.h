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
        void UniformDistribution(int min, int max)
        {
            m_uniform_real = std::uniform_real_distribution<double>(min, max);
        }
        void UniformDistribution(double min, double max)
        {
            m_uniform_real = std::uniform_real_distribution<double>(min, max);
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

/*
 * Feedback on header for future work:
 *
 * hmm a few personal notes, I would use initialization list, delegating ctors, put the public stuff first
 * (private stuff is impl detail you typically "dont need to see"), not sure why you initialize the
 * distribution in a method as opposed to ctor? (especially since it takes 0 params):
 *
 * http://coliru.stacked-crooked.com/a/8eda9cd3013fdc7e
 *
 * #include <iostream>
 * #include <random>
 *
 * template <typename Engine>
 * struct Random
 * {
 *     using Seed = typename Engine::result_type;
 *     Random(Seed seed, double min, double max) : m_seed(seed), m_engine(seed), m_distribution(min, max) {}
 *     Random() : Random(std::random_device()(), 0.0, 1.0) {}
 *
 *     double operator()() { return m_distribution(m_engine); }
 *
 * private:
 *     Seed m_seed;
 *     Engine m_engine;
 *     std::uniform_real_distribution<double> m_distribution;
 * };
 *
 * int main()
 * {
 *     std::cout << Random<std::mt19937>()() << std::endl;
 *     std::cout << Random<std::mt19937>(1337, 2.0, 3.0)() << std::endl;
 *     return 0;
 * }
 */
