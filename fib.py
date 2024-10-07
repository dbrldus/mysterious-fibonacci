import random
import math

def getFib(n):
    fib = [1,1]
    for i in range(2,n):
        fib.append(fib[i-2] + fib[i-1])
    return fib

def primesbelow(N):
    # http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    #""" Input N>=6, Returns a list of primes, 2 <= p < N """
    correction = N % 6 > 1
    N = {0:N, 1:N-1, 2:N+4, 3:N+3, 4:N+2, 5:N+1}[N%6]
    sieve = [True] * (N // 3)
    sieve[0] = False
    for i in range(int(N ** .5) // 3 + 1):
        if sieve[i]:
            k = (3 * i + 1) | 1
            sieve[k*k // 3::2*k] = [False] * ((N//6 - (k*k)//6 - 1)//k + 1)
            sieve[(k*k + 4*k - 2*k*(i%2)) // 3::2*k] = [False] * ((N // 6 - (k*k + 4*k - 2*k*(i%2))//6 - 1) // k + 1)
    return [2, 3] + [(3 * i + 1) | 1 for i in range(1, N//3 - correction) if sieve[i]]

smallprimeset = set(primesbelow(100000))
_smallprimeset = 100000
def isprime(n, precision=7):
    # http://en.wikipedia.org/wiki/Miller-Rabin_primality_test#Algorithm_and_running_time
    if n < 1:
        raise ValueError("Out of bounds, first argument must be > 0")
    elif n <= 3:
        return n >= 2
    elif n % 2 == 0:
        return False
    elif n < _smallprimeset:
        return n in smallprimeset


    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    for repeat in range(precision):
        a = random.randrange(2, n - 2)
        x = pow(a, d, n)
    
        if x == 1 or x == n - 1: continue
    
        for r in range(s - 1):
            x = pow(x, 2, n)
            if x == 1: return False
            if x == n - 1: break
        else: return False

    return True

# https://comeoncodeon.wordpress.com/2010/09/18/pollard-rho-brent-integer-factorization/
def pollard_brent(n):
    if n % 2 == 0: return 2
    if n % 3 == 0: return 3

    y, c, m = random.randint(1, n-1), random.randint(1, n-1), random.randint(1, n-1)
    g, r, q = 1, 1, 1
    while g == 1:
        x = y
        for i in range(r):
            y = (pow(y, 2, n) + c) % n

        k = 0
        while k < r and g==1:
            ys = y
            for i in range(min(m, r-k)):
                y = (pow(y, 2, n) + c) % n
                q = q * abs(x-y) % n
            g = gcd(q, n)
            k += m
        r *= 2
    if g == n:
        while True:
            ys = (pow(ys, 2, n) + c) % n
            g = gcd(abs(x - ys), n)
            if g > 1:
                break

    return g

smallprimes = primesbelow(1000) # might seem low, but 1000*1000 = 1000000, so this will fully factor every composite < 1000000
def primefactors(n, sort=False):
    factors = []

    for checker in smallprimes:
        while n % checker == 0:
            factors.append(checker)
            n //= checker
        if checker > n: break

    if n < 2: return factors

    while n > 1:
        if isprime(n):
            factors.append(n)
            break
        factor = pollard_brent(n) # trial division did not fully factor, switch to pollard-brent
        factors.extend(primefactors(factor)) # recurse to factor the not necessarily prime factor returned by pollard-brent
        n //= factor

    if sort: factors.sort()

    return factors

def factorization(n):
    factors = {}
    for p1 in primefactors(n):
        try:
            factors[p1] += 1
        except KeyError:
            factors[p1] = 1
    return factors

totients = {}
def totient(n):
    if n == 0: return 1

    try: return totients[n]
    except KeyError: pass

    tot = 1
    for p, exp in factorization(n).items():
        tot *= (p - 1)  *  p ** (exp - 1)

    totients[n] = tot
    return tot

def gcd(a, b):
    if a == b: return a
    while b > 0: a, b = b, a % b
    return a

def lcm(a, b):
    return abs((a // gcd(a, b)) * b)



#################################
#########위쪽 코드는 퍼옴#########
#################################

from copy import deepcopy
from time import time

t_fac = 0; t_eigen = 0
# M, N = map(int, input().split())
# N = int(input())
# M = 1, N = 160
N = 300

fib = getFib(N)


import numpy as np
import math

AAAAA = []  # 점화식 계수
s = []  # s[k][n] = f_(k+1)(n+1)/f_(n+1)
for k in range(1,int((N/2) ** 0.5)):
    s.append([])
    for i in range(1,N // k):
        s[k-1].append((fib[k*i - 1] // fib[i - 1]))

# 확정된 값
AAAAA.append([1])
AAAAA.append([1,1])
AAAAA.append([-1,2,2])
AAAAA.append([-1,-3,6,3])
AAAAA.append([1,-5,-15,15,5])
AAAAA.append([1,8,-40,-60,40,8])
AAAAA.append([-1,13,104,-260,-260,104,13])
AAAAA.append([-1,-21,273,1092,-1820,-1092,273,21])
AAAAA.append([1,-34,-714,4641,12376,-12376,-4641,714,34])
AAAAA.append([1,55,-1870,-19635,85085,136136,-85085,-19635,1870,55])
AAAAA.append([-1,89,4895,-83215,-582505,1514513,1514513,-582505,-83215,4895,89])

# #추측한 규칙으로 계산
# for k in range(1,int((N/2)**0.5) + 1):
#     coef = [1] * k
#     for i in range(1,k//2 + 1):
#         coef[i] = int(np.prod(fib[k-i:k])) // int(math.factorial(i-1))
#     for i in range(k//2 + 1, k):
#         coef[i] = coef[k-i]
#     for i in range(0, k):
#         if (i // 2) % 2 == 1:
#             coef[k-i-1] *= -1
#     AAAAA.append(coef)

# AAAAA 출력
AAAAAstr = []
for i in AAAAA:
    sss = ""
    for a in i:
        if a % 1 == 0:
            # s += "int "
            sss += format(a,'.0f') + ","
        else:
            sss += format(a,'.1f') + ","
    sss += "\b"
    AAAAAstr.append(sss)
for i in range(0, len(AAAAAstr)):
    AAAAAstr[i] = f"{AAAAAstr[i]:^{len(AAAAAstr[-1])}}"

for sss in AAAAAstr:
    print(sss)

# 오차 계산, 
while(True):
    uinput = int(input())
    print(s[uinput-1][0:2*uinput+2])
    print(AAAAA[uinput-1])
    for k in range(0,uinput):
        r = 0
        for i in range(0,uinput):
            r += int(AAAAA[uinput-1][i]) * int(s[uinput-1][k+i])
        print(int(r))
        print(s[uinput-1][k+uinput])
        print(abs(int(r)-s[uinput-1][k+uinput]), "\n")


####################################################################
##################### Fibonacci Factorization ######################
####################################################################

# fibFact = [{},{},{2:1}]
# fibFactOld = [{},{},{}]
# fibFactEigen = deepcopy(fibFact)
# for i in range(3,N):
#     print(f'{i+1:>3}', end = '')

#     #optimization(old, eigen) timer
#     t_start = time()
#     fibFactOld.append({})
#     for j in factorization(i+1).keys():
#         for p,r in fibFact[(i+1)//j - 1].items():
#             try:
#                 fibFactOld[i][p] = max(fibFactOld[i][p], r)
#             except KeyError:
#                 fibFactOld[i][p] = r
#     k = fib[i]
#     for p,r in fibFactOld[i].items():
#         k //= p ** r
#     t_end = time()
#     t_eigen += t_end - t_start

#     #factorization timer
#     t_start = time()
#     fibFactEigen.append(factorization(int(k)))
#     t_end = time()
#     t_fac += t_end - t_start

#     fibFact.append(deepcopy(fibFactOld[i]))
#     for p,r in fibFactEigen[i].items():
#         try:
#             fibFact[i][p] += r
#         except KeyError:
#             fibFact[i][p] = r
#     print('\b\b\b', end = '')


# #################출력부##################
# print("Factorization of Fibonacci series")
# for i in range(M-1,N):
#     print(f"-{i+1:>3}: {fib[i]:>28} \n {fibFact[i]}")
# print('-'*60)

# print("lcm of fibonacci for all divisors of index of fibonacci")
# for i in range(M-1,N):
#     print(f"{i+1:>03}: {fib[i]:>28} \n {fibFactOld[i]}")
# print('-'*60)

# print("Factorization of Eigen factors")
# for i in range(M-1,N):
#     print(f"{i+1:>3}| {'*'*sum(fibFactEigen[i].values()):>4}| {fibFactEigen[i]}")
#     # print(f"{i+1:>3}| {'*'*sum(fibFactEigen[i].values()):>4}| {int(math.prod([p**r for p,r in fibFactEigen[i].items()]))}")


# print(f"factorization time taken: {t_fac:.5f}s")
# print(f"optimization time taken : {t_eigen:.5f}s")


# for p in range(2, N):
#     s = 0
#     for i in range(0,N):
#         if(fib[i] % p == 0):
#             s = i+1
#             break
#     if s == 0: 
#         print(f"{p:>3}| too small data")
#         continue
#     for i in range(s-1, N-s, s):
#         if(fib[i] % p == 0 and fib[i+1] % p == 1):
#             print(f"{p:>3}| {(i+1) // s}")
#             break
