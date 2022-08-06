#pragma once
#include <type_traits>

namespace tml
{
    namespace Impl
    {
        constexpr double sqrt(double x, double current, double previous) noexcept
        {
            return current == previous ? current : Impl::sqrt(x, 0.5 * (current + x / current), current);
        }

        constexpr double pow(double x, double current, unsigned long y) noexcept
        {
            return y > 0 ? Impl::pow(x, current * x, y - 1) : current;
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

    template<unsigned int N, typename T>
    class Vector
    {
    public:
        constexpr Vector() noexcept
        {
            static_assert(std::is_arithmetic<T>::value, "T has to be an arithmetic type");
            static_assert(N > 1, "N has to be greater than 1");
        }

        template<typename ... A>
        constexpr explicit Vector(const A& ... args) noexcept
        : m_data{args ...}
        {
            static_assert(std::is_arithmetic<T>::value, "T has to be an arithmetic type");
            static_assert(N > 1, "N has to be greater than 1");
        }

        constexpr T& operator[](unsigned int index) noexcept
        {
            return m_data[index];
        }

        constexpr Vector<N, T> operator+(const Vector<N, T>& other) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[i] + other.m_data[i];
            }

            return result;
        }

        constexpr Vector<N, T> operator+(const T& scalar) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[i] + scalar;
            }

            return result;
        }

        constexpr Vector<N, T> operator-(const Vector<N, T>& other) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[i] - other.m_data[i];
            }

            return result;
        }

        constexpr Vector<N, T> operator-(const T& scalar) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[i] - scalar;
            }

            return result;
        }

        constexpr Vector<N, T> operator*(const Vector<N, T>& other) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[i] * other.m_data[i];
            }

            return result;
        }

        constexpr Vector<N, T> operator*(const T& scalar) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[i] * scalar;
            }

            return result;
        }

        constexpr Vector<N, T> operator/(const Vector<N, T>& other) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[i] / other.m_data[i];
            }

            return result;
        }

        constexpr Vector<N, T> operator/(const T& scalar) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[i] / scalar;
            }

            return result;
        }

        constexpr double Length() const noexcept
        {
            double length = 0;

            for(unsigned int i = 0; i < N; ++i)
            {
                length += pow(m_data[i], 2);
            }

            return sqrt(length);
        }

        constexpr Vector<N, T> Normalized() const noexcept
        {
            return *this / Length();
        }

        constexpr double Dot(const Vector<N, T>& other) const noexcept
        {
            double result = 0;

            for(unsigned int i = 0; i < N; ++i)
            {
                result += m_data[i] * other.m_data[i];
            }

            return result;
        }

        constexpr Vector<N, T> Cross(const Vector<N, T>& other) const noexcept
        {
            Vector<N, T> result;

            for(unsigned int i = 0; i < N; ++i)
            {
                result.m_data[i] = m_data[(1+i) % N] * other.m_data[(2+i) % N] - m_data[(2+i) % N] * other.m_data[(1+i) % N];
            }

            return result;
        }

    protected:
        T m_data[N];
    };
}