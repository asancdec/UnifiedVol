// SPDX-License-Identifier: Apache-2.0

namespace uv::models::svi::detail
{

template <opt::nlopt::Algorithm Algo> void fillCalendarMContext(
    CalendarMContext& mctx,
    Vector<double>& logKFBuf,
    Vector<double>& prevWkBuf,
    std::size_t numStrikes,
    const Params<double>& prevParams,
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    std::span<const double> logKF,
    const SliceData& sliceData
) noexcept
{
    const std::size_t m = numStrikes + 2;

    logKFBuf.resize(m);
    prevWkBuf.resize(m);

    const double eps{optimizer.eps()};
    const double atm{sliceData.atmTotalVariance};
    constexpr double DELTA{0.15};

    for (std::size_t i = 0; i < numStrikes; ++i)
    {
        const double k = logKF[i];
        logKFBuf[i] = k;

        prevWkBuf[i] = totalVariance(
            prevParams.a,
            prevParams.b,
            prevParams.rho,
            prevParams.m,
            prevParams.sigma,
            k
        );
    }

    {
        const double kL = sliceData.logKFMin - DELTA;
        logKFBuf[numStrikes] = kL;
        prevWkBuf[numStrikes] = totalVariance(
            prevParams.a,
            prevParams.b,
            prevParams.rho,
            prevParams.m,
            prevParams.sigma,
            kL
        );
    }

    {
        const double kR = sliceData.logKFMax + DELTA;
        logKFBuf[numStrikes + 1] = kR;
        prevWkBuf[numStrikes + 1] = totalVariance(
            prevParams.a,
            prevParams.b,
            prevParams.rho,
            prevParams.m,
            prevParams.sigma,
            kR
        );
    }

    mctx.logKF = logKFBuf;
    mctx.prevWk = prevWkBuf;
    mctx.eps = eps;
    mctx.atmTotalVariance = atm;
}

} // namespace uv::models::svi::detail
