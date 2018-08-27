
#include <stdio.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <stdint.h>
#include <iostream>
#include <cmath>
#include <cstdarg>
#include "printing.h"

using namespace std;

void readPrimes(vector<uint64_t> &primes, uint64_t L)
{
    ifstream input("./primes/total.txt");
    if (!input)
    {
        cout << "Read primes failed";
        exit(EXIT_FAILURE);
    }

    uint64_t p;
    while (input >> p)
    {

        if (p > L)
            break;
        primes.push_back(p);
    }
    cout << "Finished reading primes" << endl;
    input.close();
}

// modular exponentiation : b ^ e % mod
uint64_t mod_exp(uint64_t b, uint64_t e, uint64_t mod)
{
    uint64_t r = 1;
    while (e > 0)
    {
        if ((e & 1) == 1)
            r = (r * b) % mod;
        b = (b * b) % mod;
        e >>= 1;
    }
    return r;
}

struct Cache
{
    uint64_t invfac;
    uint64_t invfac2;
};

std::ostream &operator<<(std::ostream &os, const Cache &v)
{
    os << "{" << v.invfac << "," << v.invfac2 << "}";
    return os;
}

void getFactorials(vector<uint64_t> &facts, uint64_t mod)
{
    cout << "Started factorials" << endl;
    facts.push_back(1);
    for (auto i = 2; i < facts.capacity() - 1; i++)
    {
        facts.push_back((i * facts[i - 1]) % mod);
    }
    cout << "Finished factorials." << endl;
}

void getInverses(unordered_map<uint64_t, Cache> &map, vector<uint64_t> &facts, vector<uint64_t> &primes, uint64_t mod)
{
    cout << "Inverses: " << endl;
    int counter = 0;
    for (auto p : primes)
    {
        auto c = Cache();
        c.invfac = mod_exp(facts[p], mod - 2, mod);
        c.invfac2 = mod_exp(facts[2 * p], mod - 2, mod);
        if(c.invfac == 0 || c.invfac2 == 0)
            cout << c << " on "<< p << " facts=" << facts[p] << " f2=" << facts[2*p] << endl;
        map[p] = c;
        counter += 1;
        if (counter % 10000 == 0)
            cout << "\rArgument number\t" << counter << "/" << primes.size() << endl;
    }
    cout << "Finished memoization" << endl;
}

uint64_t modprod(uint64_t mod, uint64_t c1 = 1, uint64_t c2 = 1, uint64_t c3 = 1)
{
    return (c1 * c2 % mod) * c3 % mod;
}

uint64_t counter = 0;
uint64_t A_pre_computed(vector<uint64_t> &facts, vector<uint64_t> &primes, unordered_map<uint64_t, Cache> &map, uint64_t q, uint64_t n, uint64_t mod)
{
    counter += 1;
    if (counter % 10000 == 0)
        cout << "\rArgument number\t" << counter << "/" << primes.size() << endl;
    // 1 / n((qn choose n) +(n - 1) q(mod p)
    //(nq) !+(n - 1) qn !((q - 1) n) !
    //                                / n * n !((q - 1) n) !
    //((nq) !+(n - 1) qn !((q - 1) n) !) *(n - 1) !)
    //                                      *(inv(n !) ^ 2 * inv(((q - 1) n) !))

    uint64_t num, den;
    if (q == 2)
    {
        num =
            (facts[n * q] + modprod(mod, (n - 1) * q, facts[n], facts[(q - 1) * n])) * facts[n - 1] % mod;
        den = modprod(mod, map.at(n).invfac, map.at(n).invfac, map.at(n).invfac);
    }
    else
    {
        num =
            (facts[n * q] + modprod(mod, (n - 1) * q, facts[n], facts[(q - 1) * n])) * facts[n - 1] % mod;
        den = modprod(mod, map.at(n).invfac, map.at(n).invfac, map.at(n).invfac2);
    }
    // cout<<"num="<<num<<" den"<<den<<endl;
    return num * den % mod;
}

uint64_t S_pre_computed(vector<uint64_t> &facts, vector<uint64_t> &primes, unordered_map<uint64_t, Cache> &map, uint64_t L, uint64_t q, uint64_t mod)
{
    uint64_t sum = 0;
    for (auto n : primes)
    {
        if (n > L)
            break;
        sum = (sum + A_pre_computed(facts, primes, map, q, n, mod)) % mod;
    }
    return sum - q;
}

int main(int argc, char *args[])
{
    uint64_t L = 8e8;
    uint64_t mod = 1000000009; // prime

    vector<uint64_t> primes;
    primes.reserve(round(L / log(L)));
    readPrimes(primes, L);

    vector<uint64_t> facts(1);
    facts.reserve(3 * L);
    cout<<"Allocated";
    getFactorials(facts, mod);

    unordered_map<uint64_t, Cache> map;
    map.reserve(primes.size());
    getInverses(map, facts, primes, mod);
    
    cout<<"S_2(1e8)="<< S_pre_computed(facts, primes, map, L, 2, mod)<< endl;
    uint64_t n = 11;
    // cout << "Total sum:" << S_pre_computed(facts, primes, map, L, 3, mod) << endl;
}
