from modules.fib_tools import *
from copy import deepcopy
from time import time

####################################################################
##################### Fibonacci Factorization ######################
####################################################################

M, N = 1, 100
t_fac = 0; t_eigen = 0
fib = getFib(N)

fibFact = [{},{},{2:1}]
fibFactOld = [{},{},{}]
fibFactEigen = deepcopy(fibFact)
for i in range(3,N):
    print(f'{i+1:>3}', end = '')

    #optimization(old, eigen) timer
    t_start = time()
    fibFactOld.append({})
    for j in factorization(i+1).keys():
        for p,r in fibFact[(i+1)//j - 1].items():
            try:
                fibFactOld[i][p] = max(fibFactOld[i][p], r)
            except KeyError:
                fibFactOld[i][p] = r
    k = fib[i]
    for p,r in fibFactOld[i].items():
        k //= p ** r
    t_end = time()
    t_eigen += t_end - t_start

    #factorization timer
    t_start = time()
    fibFactEigen.append(factorization(int(k)))
    t_end = time()
    t_fac += t_end - t_start

    fibFact.append(deepcopy(fibFactOld[i]))
    for p,r in fibFactEigen[i].items():
        try:
            fibFact[i][p] += r
        except KeyError:
            fibFact[i][p] = r
    print('\b\b\b', end = '')


#################출력부##################
print("Factorization of Fibonacci series")
for i in range(M-1,N):
    print(f"-{i+1:>3}: {fib[i]:>28} \n {fibFact[i]}")
print('-'*60)

print("lcm of fibonacci for all divisors of index of fibonacci")
for i in range(M-1,N):
    print(f"{i+1:>03}: {fib[i]:>28} \n {fibFactOld[i]}")
print('-'*60)

print("Factorization of Eigen factors")
for i in range(M-1,N):
    print(f"{i+1:>3}| {'*'*sum(fibFactEigen[i].values()):>4}| {fibFactEigen[i]}")
    # print(f"{i+1:>3}| {'*'*sum(fibFactEigen[i].values()):>4}| {int(math.prod([p**r for p,r in fibFactEigen[i].items()]))}")


print(f"factorization time taken: {t_fac:.5f}s")
print(f"optimization time taken : {t_eigen:.5f}s")


for p in range(2, N):
    s = 0
    for i in range(0,N):
        if(fib[i] % p == 0):
            s = i+1
            break
    if s == 0: 
        print(f"{p:>3}| too small data")
        continue
    for i in range(s-1, N-s, s):
        if(fib[i] % p == 0 and fib[i+1] % p == 1):
            print(f"{p:>3}| {(i+1) // s}")
            break