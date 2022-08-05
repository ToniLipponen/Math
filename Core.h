#pragma once
#include <limits>

namespace tml
{
    namespace Impl
    {
        constexpr double sqrt(double x, double curr, double prev) noexcept
        {
            return curr == prev ? curr : Impl::sqrt(x, 0.5 * (curr + x / curr), curr);
        }

        constexpr double pow(double x, double curr, unsigned long y) noexcept
        {
            return y > 0 ? Impl::pow(x, curr * x, y - 1) : curr;
        }
    }

    constexpr double sqrt(double x) noexcept
    {
        if(x < 0)
        {
            x *= -1;
        }

        return Impl::sqrt(x, x, 0);
    }

    constexpr double pow(double x, unsigned long y) noexcept
    {
        return Impl::pow(x, 1, y);
    }
}