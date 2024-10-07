# searching the coefficient of recurrence relation.

from modules.fib_tools import *
import numpy as np
import math

N = 1800

fib = getFib(N)

AAAAA = []  # 점화식 계수
s = []  # s[k][n] = f_(k+1)(n+1)/f_(n+1)
for k in range(1,int((N/2) ** 0.5)):
    s.append([])
    for i in range(1,N // k):
        s[k-1].append((fib[k*i - 1] // fib[i - 1]))

# 확정된 값
# AAAAA.append([1])
# AAAAA.append([1,1])
# AAAAA.append([-1,2,2])
# AAAAA.append([-1,-3,6,3])
# AAAAA.append([1,-5,-15,15,5])
# AAAAA.append([1,8,-40,-60,40,8])
# AAAAA.append([-1,13,104,-260,-260,104,13])
# AAAAA.append([-1,-21,273,1092,-1820,-1092,273,21])
# AAAAA.append([1,-34,-714,4641,12376,-12376,-4641,714,34])
# AAAAA.append([1,55,-1870,-19635,85085,136136,-85085,-19635,1870,55])
# AAAAA.append([-1,89,4895,-83215,-582505,1514513,1514513,-582505,-83215,4895,89])


#추측한 규칙으로 계산
for k in range(1,int((N/2)**0.5) + 1):
    coef = [1] * k
    for i in range(1,k//2 + 1):
        coef[i] = int(np.prod(fib[k-i:k])) // int(np.prod(fib[:i]))
    for i in range(k//2 + 1, k):
        coef[i] = coef[k-i]
    for i in range(0, k):
        if (i // 2) % 2 == 1:
            coef[k-i-1] *= -1
    AAAAA.append(coef)

# AAAAA 출력
# AAAAAstr = []
# for i in AAAAA:
#     sss = ""
#     for a in i:
#         if isinstance(a, int):
#             # s += "int "
#             sss += format(a,'.0f') + ","
#         else:
#             sss += format(a,'.1f') + ","
#     sss += "\b"
#     AAAAAstr.append(sss)
# for i in range(0, len(AAAAAstr)):
#     AAAAAstr[i] = f"{AAAAAstr[i]:^{len(AAAAAstr[-1])}}"

# for sss in AAAAAstr:
#     print(sss)

# 오차 계산
while(True):
    u = int(input())
    # print(s[u-1][0:2*u+2])
    print(AAAAA[u-1])
    for k in range(0,u):
        r = 0
        for i in range(0,u):
            r += int(AAAAA[u-1][i]) * int(s[u-1][k+i])
        print(int(r))
        print(s[u-1][k+u])
        print(abs(int(r)-s[u-1][k+u]), "\n")
