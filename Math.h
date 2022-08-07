/**
 * MIT License
 *
 * Copyright (c) 2022 Toni Lipponen
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

#pragma once
#include <type_traits>
#include <cmath>

namespace tml
{
#if defined(TML_USE_DOUBLE_PRECISION)
    using float_type = double;
#else
    using float_type = float;
#endif

    template<unsigned int N, typename T>
    class Vector
    {
    public:
        constexpr Vector() noexcept
        {
            static_assert(std::is_arithmetic<T>::value, "T has to be an arithmetic type");
            static_assert(N > 0, "N has to be greater than 0");
        }

        template<typename ... A>
        constexpr explicit Vector(const A& ... args) noexcept
        : m_data{args ...}
        {
            static_assert(std::is_arithmetic<T>::value, "T has to be an arithmetic type");
            static_assert(N > 0, "N has to be greater than 0");
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

    template<unsigned int R, unsigned int C, typename T>
    class Matrix
    {
    public:
        constexpr Matrix() noexcept = default;

        template<typename ... A>
        constexpr explicit Matrix(const A& ... args) noexcept
        : m_rows{args ...}
        {

        }

        constexpr Vector<C, T>& operator[](unsigned int index) noexcept
        {
            return m_rows[index];
        }

        constexpr Vector<C, T> operator[](unsigned int index) const noexcept
        {
            return m_rows[index];
        }

        template<unsigned int R2, unsigned C2>
        constexpr Matrix<R,C2,T> operator*(const Matrix<R2,C2,T>& other) const noexcept
        {
            static_assert(C == R2, "Invalid operands");

            Matrix<R,C2,T> result{};

            for(unsigned int i = 0; i < R; ++i)
            {
                Vector<C, T> row = m_rows[i];

                for(unsigned int j = 0; j < C2; ++j)
                {
                    Vector<R2, T> col{};

                    for(unsigned int y = 0; y < R2; ++y)
                    {
                        col[y] = other[y][j];
                    }

                    result[i][j] = row.Dot(col);
                }
            }

            return result;
        }

        constexpr static Matrix<R,C,T> Identity() noexcept
        {
            static_assert(R == C, "This function is only defined for square matrices");
            Matrix<R,C,T> result{};

            for(unsigned int i = 0; i < C; ++i)
            {
                result.m_rows[i][i] = 1;
            }

            return result;
        }

        constexpr static Matrix<R,C,T> Scale(const Vector<C, T>& scale) noexcept
        {
            static_assert(R == C, "This function is only defined for square matrices");
            Matrix<R,C,T> result{};

            for(unsigned int i = 0; i < C; ++i)
            {
                result[i][i] = scale[i];
            }

            return result;
        }

        static Matrix<R,C,T> Rotate(const Vector<C,T>& axis, double r) noexcept
        {
            static_assert(R == C, "This function is only defined for square matrices");
            static_assert(R > 1 && R < 5, "Invalid matrix dimensions");

            const auto sinr = sin(r);
            const auto cosr = cos(r);

            switch(R)
            {
                case 2:
                {
                    return Matrix<2,2,float_type>(
                            Vector<2,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(-sinr)),
                            Vector<2,float_type>(static_cast<float_type>(sinr), static_cast<float_type>(cosr)));
                }

                case 3:
                {
                    if(axis[0] > 0)
                    {
                        return Matrix<3,3,float_type>(
                                Vector<3,float_type>(static_cast<float_type>(1), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<3,float_type>(static_cast<float_type>(0), static_cast<float_type>(cosr), static_cast<float_type>(-sinr)),
                                Vector<3,float_type>(static_cast<float_type>(0), static_cast<float_type>(sinr), static_cast<float_type>(cosr)));
                    }
                    else if(axis[1] > 0)
                    {
                        return Matrix<3,3,float_type>(
                                Vector<3,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(0), static_cast<float_type>(sinr)),
                                Vector<3,float_type>(static_cast<float_type>(0), static_cast<float_type>(1), static_cast<float_type>(0)),
                                Vector<3,float_type>(static_cast<float_type>(-sinr), static_cast<float_type>(0), static_cast<float_type>(cosr)));
                    }
                    else
                    {
                        return Matrix<3,3,float_type>(
                                Vector<3,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(-sinr), static_cast<float_type>(0)),
                                Vector<3,float_type>(static_cast<float_type>(sinr), static_cast<float_type>(cosr), static_cast<float_type>(0)),
                                Vector<3,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1)));
                    }
                }

                case 4:
                {
                    if(axis[0] > 0)
                    {
                        return Matrix<4,4,float_type>(
                                Vector<4,float_type>(static_cast<float_type>(1), static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(0), static_cast<float_type>(cosr), static_cast<float_type>(-sinr), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(0), static_cast<float_type>(sinr), static_cast<float_type>(cosr), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1))
                                        );
                    }
                    else if(axis[1] > 0)
                    {
                        return Matrix<4,4,float_type>(
                                Vector<4,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(0), static_cast<float_type>(sinr), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(0), static_cast<float_type>(1), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(-sinr), static_cast<float_type>(0), static_cast<float_type>(cosr), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1)));
                    }
                    else
                    {
                        return Matrix<4,4,float_type>(
                                Vector<4,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(-sinr), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(sinr), static_cast<float_type>(cosr), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1), static_cast<float_type>(0)),
                                Vector<4,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1)));
                    }
                }

                default:
                    break;
            }

            return {};
        }

        constexpr static Matrix<R,C,T> Translate(const Vector<C,T>& offset) noexcept
        {
            static_assert(R == C, "This function is only defined for square matrices");
            Matrix<R,C,T> result = Matrix<R,C,T>::Identity();

            for(unsigned int i = 0; i < C; ++i)
            {
                result[i][C-1] = offset[i];
            }

            return result;
        }

        constexpr Matrix<R,C,T> Inverse() const noexcept
        {
            Matrix<R,C,T> result{};

            return result;
        }

        const unsigned int rows = R;
        const unsigned int columns = C;

    protected:
        Vector<C, T> m_rows[R];
    };

    using Vector2 = tml::Vector<2, float_type>;
    using Vector3 = tml::Vector<3, float_type>;
    using Vector4 = tml::Vector<4, float_type>;

    using Matrix4 = Matrix<4, 4, float_type>;
    using Matrix3 = Matrix<3, 3, float_type>;
    using Matrix2 = Matrix<2, 2, float_type>;
}