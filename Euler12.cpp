#include <cmath>
#include <ctime>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;


//the sieve of Eratosthenes, the limit by both memory and time is about 10^10
// older version with step 2
// (10^10 takes about 137 seconds)
vector<unsigned long> Erat(unsigned long N){
    vector<unsigned long> primes;
    if (N == 0 || N == 1) return primes;
    primes.push_back(2);
    if (N == 2) return primes;
    primes.push_back(3);
    if (N == 3) return primes;
    vector<bool> Nbool(N+1, true);
    vector<bool>::iterator itr = Nbool.begin();
    *itr = false;
    *(itr+1) = false;
    for (itr = Nbool.begin() + 4; itr < Nbool.end(); itr += 2) *itr = false;
    unsigned long p = 3;
    while (p*p <= N) {
        for (itr = Nbool.begin() + p*p; itr < Nbool.end(); itr += p) *itr = false;
        p += 2;
        while (! Nbool[p]) p += 2;
    }
    for (itr = Nbool.begin()+5; itr < Nbool.end(); itr += 2) {
        if (*itr) primes.push_back(itr - Nbool.begin());
    }
    return primes;
}

vector<unsigned long> IncompletePrimeFactor(unsigned long N, vector<unsigned long> primes0) 
// this one with externally supplied primes, need them to be at least up to sqrt(N)
{
    vector<unsigned long>::iterator itr, idx;
    vector<unsigned long> factors;
    idx = lower_bound(primes0.begin(), primes0.end(), (unsigned long)(floor(sqrt((double) N))));
    // determine the itr below the value for sqrt(N)
    for (itr = primes0.begin(); itr != idx + 1; itr ++) // add +1 in case sqrt(N) is a prime
    {
        if (N % *itr == 0) factors.push_back(*itr);
    }
    return factors;
} 

unsigned long FactorMultiplicityExtPndiv(unsigned long N, vector<unsigned long> primes0)
// this one with externally supplied primes and only counting divisors
{
    if (N == 1) return 1; // corner case
    vector<unsigned long> factors = IncompletePrimeFactor(N, primes0);
    // vector<pair <unsigned long,int>> mfactors;
    // only counting divisors here, turn off all the factor/multiplicity collection 
    unsigned long ndiv = 1;
    if (factors.empty()) // prime
    {
        // mfactors.push_back(make_pair(N,1));
        // return mfactors;
        return 2; 
    }
    vector<unsigned long>::iterator itr;
    int m;
    for (itr = factors.begin(); itr != factors.end(); itr++) {
        m = 0;
        while (N % *itr == 0) {
            m++;
            N /= *itr;
        }
        // mfactors.push_back(make_pair(*itr,m));
        ndiv *= m + 1;    
    }
    if (N != 1) {
        ndiv *= 2; // case of one additional prime factor > sqrt(N)
        // mfactors.push_back(make_pair(N, 1));
    }
    return ndiv;
}

int main() {
    int limdiv = 500;
    int ndiv = 0;
    unsigned long n = 1;
    clock_t c_start = clock();
    vector<unsigned long> primes = Erat(1000);
    while (ndiv < limdiv) {
        if (n % 2 == 0) ndiv = FactorMultiplicityExtPndiv(n/2, primes)*FactorMultiplicityExtPndiv(n+1, primes);
        else ndiv = FactorMultiplicityExtPndiv(n, primes)*FactorMultiplicityExtPndiv((n+1)/2, primes);
        n++;
    }
    n--;
    clock_t c_end = clock();
    cout << endl;
    cout << "n = " << n << " T = " << (n+1)*n/2 << " number of divisors = " << ndiv << endl;
    cout << 1000.0*(c_end - c_start) / CLOCKS_PER_SEC << " ms\n";
    return 0;
}