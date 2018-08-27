#!/usr/bin/env python3

from primes import *
# modular exponentiation: b^e % mod


def mod_exp(b, e, mod):
    r = 1
    while e > 0:
        if (e & 1) == 1:
            r = (r*b) % mod
        b = (b*b) % mod
        e >>= 1
    return r


def precompute(maximal, mod):
    print("Started factorials")
    facts = [0]*3*maximal
    facts[1] = 1
    for i in range(2, len(facts)-1):
        facts[i] = (i * facts[i-1]) % mod
    print("Finished factorials. Primes: ")
    counter = 0

    primes = primes_dict_from_file(maximal)
    for prim, val in primes.items():
        counter += 1
        if(counter % 10000 == 0):
            print("\rArgument number\t", counter, "/", len(primes))
        val.invfac = mod_exp(facts[prim], mod-2, mod)
        val.invfac2 = mod_exp(facts[2 * prim], mod-2, mod)

    print("Finished memoization")
    return facts, primes


def modprod(mod, *wargs):
    ans = 1
    for a in wargs:
        ans = ans * a % mod
    return ans


counter = 0


def A_pre_computed(facts, primes, n, q, mod):
    counter += 1
    if(counter % 10000 == 0):
        print("\rArgument number\t", counter, "/", len(primes))
    # 1/n ( (qn choose n) + (n-1)q (mod p)
    # (nq)! + (n-1)qn! ((q-1)n)!
    #                           / n * n! ((q-1)n)!
    # ((nq)! + (n-1)qn! ((q-1)n)!) * (n-1)! )
    #                                     * ( inv(n!)^2 * inv(((q-1)n)!) )
    if(q == 2):
        num = (facts[n*q] + modprod(mod, (n-1)*q, facts[n],
                                    facts[(q-1)*n])) * facts[n-1] % mod
        den = (primes[n].invfac**2 % mod) * primes[n].invfac % mod
    else:
        num = (facts[n*q] + modprod(mod, (n-1)*q, facts[n],
                                    facts[(q-1)*n])) * facts[n-1] % mod
        den = (primes[n].invfac**2 % mod) * primes[n].invfac2 % mod
    return num * den % mod


def S_pre_computed(facts, primes, L, q, mod):
    sum = 0
    for n in sorted(primes):
        if n > L:
            break
        sum = (sum + A_pre_computed(facts, primes, n, q, mod)) % mod
    return sum - q


if __name__ == '__main__':
    L = 10**8
    mod = 1000000009  # prime

    facts, primes = precompute(L, mod)
    print("Memoization finished")
    # print(A_pre_computed(facts, primes, 5, 3, mod))
    n = 11
    print("Total sum:", S_pre_computed(facts, primes, L, 3, mod))

    # pre-compute factorials and inverses
    # facts, invfacts = fermat_compute(n,mod)

    # print (950 choose 100) mod 1000000007 (with pre-computing)
    # print(binom_pre_computed(facts,invfacts,950,100,mod)) # should be 640644226

    # (950 choose 100) mod 1000000007 (without pre-computing)
    # print(fermat_binom(950,100,mod)) # should be 640644226
