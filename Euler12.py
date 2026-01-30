import numpy as np
from functools import reduce
from time import perf_counter

def EratMa(N):
    """
    The Sieve of Eratosthenes implemented with numpy arrays
    Returns the list of primes lower or equal than N
    """
    if N < 11:
        first_primes = ([],[],[2],[2,3],[2,3],[2,3,5],[2,3,5],[2,3,5,7],[2,3,5,7],[2,3,5,7],[2,3,5,7])
        return(first_primes[N]) 
    N_bool = np.array([True]*(N+1))
    N_bool[0] = False
    N_bool[1] = False
    N_bool[2**2::2] = False # eiminating even numbers
    N_bool[3**2::3] = False # eiminating multiples of 3
    p = 5
    while p*p <= N:
        # start from p**2
        # because all the smaller composites have factors <p
        # and are already eliminated in previous steps
        if N_bool[p]:
            N_bool[p**2::p] = False # python doesn't care if ::p goes beyond existing array
        if N_bool[p+2]:
            N_bool[(p+2)**2::(p+2)] = False
        p = p + 6 # we can move in steps of 6
    # returns indices of nonzero elements, which in this case
    # ARE the correcponding natural numbers, which were not eliminated
    # Have to use index [0], because for technical reasons it produces a 2D array
    return np.nonzero(N_bool)[0] # in this version returns np.array


def IncompletePrimeFactorExtP(N, primes0): 
    """
    Returns the list of all prime factors of N with a possible exception
    of a single factor > sqrt(N). If N is prime, it returns an empty list.
    
    In this version primes must be proveded as an external np.array
    that goes at least up to sqrt(N)
    """
    # create the list of primes lower or equal to sqrt(N)
    
    idx = np.searchsorted(primes0, int(np.floor(np.sqrt(N))), side="left") 
    # np.searchsorted(a, v, side = "left/right")
    # returns index i corresponding to:
    # left: a[i-1] < v <= a[i]
    # right: a[i-1] <= v < a[i]
    
    # it does not matter if we include a signle prime larger than sqrt(N)
    # so there is no need to fine tune to idx or idx-1
    
    primes = primes0[:idx+1]
    
    # if N is divisible by a prime, N%p = 0
    # so when we create an array of N%p for all primes
    # it has 0 entries for factors and non-zero one for non-factors
    # recasting it as bool turns 0 into False and non-zero to True
    # after we invert them, we get True for factors and False for non-factors
    prime_mask = np.invert(np.array(N%primes, dtype = bool))
    # now we only need to apply this mask to return only prime factors
    # but since we limited ourselves to p < sqrt(N)
    # we might miss a factor > sqrt(N), like e.g. 33 = 3*11
    # so be careful when using this output
    return list(primes[prime_mask])

def FactorMultiplicityExtPndiv(N, primes0):
    """
    Returns two lists, the first with all the factors of N
    the second with corresponding multimlicities of each of them.
    And also it calculates and returns the number of divisors
    
    In this version primes must be proveded as an external np.array
    that goes at least up to sqrt(N)
    """
    factors = IncompletePrimeFactorExtP(N, primes0)
    if not factors: return [N],[1],2 # N is prime
    multiplicity = []
    # just cycle through all the factors 
    # and check how many time N is divisible by each
    for p in factors:
        m = 0
        while N%p == 0:
            m += 1
            N = N//p
        multiplicity.append(m)
    if N != 1: # the case of a single additional factor > sqrt(N)
        factors.append(N)
        multiplicity.append(1)
        
    # must recast them as python int, because
    # since we used numpy array before, their
    # type was changed to long_scalars, the default
    # int type of numpy. If we don't do this
    # we quickly get overflow when multiplying these
    # numbers
    
    factors = [int(x) for x in factors]
    # number of divisors for multiplicities m1, m2, m3... is (m1+1)*(m2+1)*(m3+1)...
    ndiv = reduce(lambda x, y: x*y, [multiplicity[x]+1 for x in range(len(multiplicity))], 1)

    return factors, multiplicity, ndiv


start = perf_counter()  
limdiv = 500
primes = EratMa(1000) # only need up to sqrt(n)
n = 1
ndiv = 0
while ndiv < limdiv: # n up to 1 000 000 means T up to 500 000*1 000 001 = 500 000 500 000 ~ 5*10**11

    # T = n(n+1)/2
    # either n or n+1 is even
    # and they don't have common divisors
    
    if n % 2 == 0:
        ndiv = FactorMultiplicityExtPndiv(n+1, primes)[2] * FactorMultiplicityExtPndiv(n//2, primes)[2]
        # there are no common prime factors in numbers n and n+1 
        # (can be proved by the sieve of Eratosthenes logic)
        # so the multiplicities are independent
        # and the numbers of divisors are simply multiplied by each other
        # to get the number of divisors for the composite
    else:
        ndiv = FactorMultiplicityExtPndiv(n, primes)[2] * FactorMultiplicityExtPndiv((n + 1)//2, primes)[2]
    
    n += 1;
    
end = perf_counter()

print("n =", n, "T =", n*(n+1)//2, "number of divisors =", ndiv)
print(end - start,'sec')

