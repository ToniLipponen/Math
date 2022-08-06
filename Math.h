#pragma once
#include <type_traits>
#include <cmath>

namespace tml
{
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

        constexpr T operator[](unsigned int index) const noexcept
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

    using Vector2f = tml::Vector<2, float>;
    using Vector3f = tml::Vector<3, float>;
    using Vector4f = tml::Vector<4, float>;

    using Vector2i = tml::Vector<2, int>;
    using Vector3i = tml::Vector<3, int>;
    using Vector4i = tml::Vector<4, int>;

    using Vector2d = tml::Vector<2, double>;
    using Vector3d = tml::Vector<3, double>;
    using Vector4d = tml::Vector<4, double>;

    template<unsigned int R, unsigned int C, typename T>
    class Matrix
    {
    public:
        constexpr Vector<C, T>& operator[](unsigned int index) noexcept
        {
            return m_rows[index];
        }

        constexpr static Matrix<R,C,T> Identity() noexcept
        {
            Matrix<R,C,T> result{};

            for(unsigned int i = 0; i < C; ++i)
            {
                result.m_rows[i][i] = 1;
            }

            return result;
        }

        constexpr static Matrix<R,C,T> Scale(const Vector<C, T>& scale) noexcept
        {
            Matrix<R,C,T> result{};

            for(unsigned int i = 0; i < C; ++i)
            {
                result[i][i] = scale[i];
            }

            return result;
        }

        static Matrix<R,C,T> Rotate(const Vector<C,T>& axis, double r) noexcept
        {
            const double sinr = sin(r);
            const double cosr = cos(r);

            return {};
        }

        constexpr static Matrix<R,C,T> Translate(const Vector<C,T>& offset) noexcept
        {
            Matrix<R,C,T> result = Matrix<R,C,T>::Identity();

            for(unsigned int i = 0; i < C; ++i)
            {
                result[i][C-1] = offset[i];
            }

            return result;
        }

    private:
        Vector<C, T> m_rows[R];
    };

    using Matrix4f = Matrix<4,4,float>;
    using Matrix2f = Matrix<2,2,float>;
}