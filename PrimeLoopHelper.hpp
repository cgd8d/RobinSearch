// This file can be included multiple times and enables
// a critical loop to be optimized separately for hot and
// cold instantiations.

FastBigFloat<3> lhs_update_rndd_test = lhs_update_rndd;
FastBigFloat<3> lhs_update_rndu_test = lhs_update_rndu;
FastBigFloat<3> rhs_update_rndd_test = rhs_update_rndd;
FastBigFloat<3> rhs_update_rndu_test = rhs_update_rndu;

while(NextPrimeIdx + BunchSize - 1 <= MaxBunchIdx)
{
    cnt_FastBunchMul++;

    _Pragma(ROBINSEARCH_PRAGMA_HINT)
    for(size_t i = NextPrimeIdx;
        i < NextPrimeIdx + BunchSize;
        i++)
    {
        lhs_update_rndd_test.mul_ui_rndd(PrimeQueue[i]+1);
        lhs_update_rndu_test.mul_ui_rndu(PrimeQueue[i]+1);
        rhs_update_rndd_test.mul_ui_rndd(PrimeQueue[i]);
        rhs_update_rndu_test.mul_ui_rndu(PrimeQueue[i]);
    }

    // Check if test values indicate possible violation of bound.
    lhs_update_rndu_test.get_rndu(mpfr_helper.a);
    mpfr_mul(mpfr_helper.a, mpfr_helper.a, LHS_rndu, MPFR_RNDU);
    rhs_update_rndd_test.get_rndd(mpfr_helper.c);
    mpfr_mul(mpfr_helper.b, mpfr_helper.c, NloglogN_rndd, MPFR_RNDD);

    if(mpfr_less_p(mpfr_helper.a, mpfr_helper.b))
    {
        // LHS < RHS is guaranteed.
        // Save current progress and keep going.
        lhs_update_rndd = lhs_update_rndd_test;
        lhs_update_rndu = lhs_update_rndu_test;
        rhs_update_rndd = rhs_update_rndd_test;
        rhs_update_rndu = rhs_update_rndu_test;
        NextPrimeIdx += BunchSize;
        cnt_NumUniquePrimeFactors += BunchSize;
        Number_factors.back().PrimeHi = PrimeQueue[NextPrimeIdx-1];
        cnt_NumPrimeFactors += BunchSize;
        cnt_FastBunchMul_keep++;
    }
    else
    {
        // Possibly LHS >= RHS.
        // We need to drop that last bunch and go more carefully.
        // Reduce MaxBunchIdx so we won't try the same thing again.
        MaxBunchIdx = NextPrimeIdx + BunchSize - 2;
        break;
    }
}
