# RobinSearch

Superabundant numbers take the form:

N = \prod_{i=0}^{m-1} p_i^{k_i}

Where p_0=2, p_1=3, etc., and k_i \ge k_{i+1}.

Additionally, it is known that p_{m-1} \approx log(N) and for large numbers p_{m-1} < p_0^{k_0} = 2^{k_0}.

From theorem 5 of Alaoglu, roughly, p log p (1/(2^(k+1)-2)) > log 2, so p log p > log 2 (2^(k+1)-2), so say p log p \gtrapprox (log 2) 2^(k+1).

We represent exponents k_i using uint8_t and prime factors using uint64_t.  Based on the above two inequalities it is clear that p<2^64 is the more constraining of the two.

When p_{m-1} \approx 2^64 then N \approx exp(2^64) \approx 2^(1.44 * 2^64).  We rely on mpfr to compute log(N) and the maximum exponent is 2^62-1 so the largest possible mpfr value is (1-epsilon)*2^(2^62-1).  So the limiting factor is actually the size that can be held by mpfr. If we compute sigma(N)/exp(gamma) on one hand and N log log N on the other, then we are limited by roughly N log log N < 2^(2^62-1), or roughly log N \approx p_{m-1} < 2^62 log 2.  Still this is good enough.


References:
Alaoglu, L., and P. Erdos. “On Highly Composite and Similar Numbers.” Transactions of the American Mathematical Society 56, no. 3 (1944): 448–69. https://doi.org/10.2307/1990319.
