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

#define TML_DEFINE_VECTOR_OPERATOR_SCALAR(op)                       \
constexpr Vector<N, T> operator op (const T& scalar) const noexcept \
{                                                                   \
    Vector<N, T> result;                                            \
    for(unsigned int i = 0; i < N; ++i)                             \
    {                                                               \
        result.m_data[i] = m_data[i] op scalar;                     \
    }                                                               \
    return result;                                                  \
}

#define TML_DEFINE_VECTOR_OPERATOR_VECTOR(op)                                   \
constexpr Vector<N, T> operator op (const Vector<N, T>& vector) const noexcept  \
{                                                                               \
    Vector<N, T> result;                                                        \
    for(unsigned int i = 0; i < N; ++i)                                         \
    {                                                                           \
        result.m_data[i] = m_data[i] op vector.m_data[i];                       \
    }                                                                           \
    return result;                                                              \
}

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

        TML_DEFINE_VECTOR_OPERATOR_SCALAR(+);
        TML_DEFINE_VECTOR_OPERATOR_SCALAR(-);
        TML_DEFINE_VECTOR_OPERATOR_SCALAR(*);
        TML_DEFINE_VECTOR_OPERATOR_SCALAR(/);

        TML_DEFINE_VECTOR_OPERATOR_VECTOR(+);
        TML_DEFINE_VECTOR_OPERATOR_VECTOR(-);
        TML_DEFINE_VECTOR_OPERATOR_VECTOR(*);
        TML_DEFINE_VECTOR_OPERATOR_VECTOR(/);

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

        constexpr static unsigned int elements = N;
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

        constexpr Matrix<R,C,T> operator*(T scalar) const noexcept
        {
            Matrix<R,C,T> result{};

            for(unsigned int i = 0; i < R; ++i)
            {
                for(unsigned int j = 0; j < C; ++j)
                {
                    result[i][j] = m_rows[i][j] * scalar;
                }
            }

            return result;
        }

        constexpr Vector<R,T> operator*(const Vector<R,T>& vector) const noexcept
        {
            Matrix<R, 1, T> matrix;

            for(unsigned int i = 0; i < R; ++i)
            {
                matrix[i][0] = vector[i];
            }

            matrix = *this * matrix;

            Vector<R,T> result{};

            for(unsigned int i = 0; i < R; ++i)
            {
                result[i] = matrix[i][0];
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

#if defined(TML_USE_DEGREES)
            r *= 0.01745329252;
#endif
            const auto sinr = sin(r);
            const auto cosr = cos(r);

            switch(R)
            {
                case 2:
                {
                    return Matrix<R,C,float_type>(
                            Vector<C,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(-sinr)),
                            Vector<C,float_type>(static_cast<float_type>(sinr), static_cast<float_type>(cosr)));
                }

                case 3:
                {
                    if(axis[0] > 0)
                    {
                        return Matrix<R,C,float_type>(
                                Vector<C,float_type>(static_cast<float_type>(1), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(cosr), static_cast<float_type>(-sinr)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(sinr), static_cast<float_type>(cosr)));
                    }
                    else if(axis[1] > 0)
                    {
                        return Matrix<R,C,float_type>(
                                Vector<C,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(0), static_cast<float_type>(sinr)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(1), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(-sinr), static_cast<float_type>(0), static_cast<float_type>(cosr)));
                    }
                    else
                    {
                        return Matrix<R,C,float_type>(
                                Vector<C,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(-sinr), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(sinr), static_cast<float_type>(cosr), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1)));
                    }
                }

                case 4:
                {
                    if(axis[0] > 0)
                    {
                        return Matrix<R,C,float_type>(
                                Vector<C,float_type>(static_cast<float_type>(1), static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(cosr), static_cast<float_type>(-sinr), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(sinr), static_cast<float_type>(cosr), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1))
                        );
                    }
                    else if(axis[1] > 0)
                    {
                        return Matrix<R,C,float_type>(
                                Vector<C,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(0), static_cast<float_type>(sinr), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(1), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(-sinr), static_cast<float_type>(0), static_cast<float_type>(cosr), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1)));
                    }
                    else
                    {
                        return Matrix<R,C,float_type>(
                                Vector<C,float_type>(static_cast<float_type>(cosr), static_cast<float_type>(-sinr), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(sinr), static_cast<float_type>(cosr), static_cast<float_type>(0), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1), static_cast<float_type>(0)),
                                Vector<C,float_type>(static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(0), static_cast<float_type>(1)));
                    }
                }
                default:
                    break;
            }
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
            static_assert(R == C, "This function is only defined for square matrices");
            static_assert(R <= 4 && R >= 2, "Invalid matrix dimensions");

            return Inverse::Get(*this);
        }

        struct Inverse
        {
            constexpr static Matrix<2,2,T> Get(const Matrix<2,2,T>& m) noexcept
            {
                T oneOverDeterminant = static_cast<T>(1) / (
                         m[0][0] * m[1][1]
                        -m[1][0] * m[0][1]);

                return Matrix<2, 2, T>(
                         m[1][1] * oneOverDeterminant,
                        -m[0][1] * oneOverDeterminant,
                        -m[1][0] * oneOverDeterminant,
                         m[0][0] * oneOverDeterminant);
            }

            constexpr static Matrix<3,3,T> Get(const Matrix<3,3,T>& m) noexcept
            {
                T oneOverDeterminant = static_cast<T>(1) / (
                          m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
                        - m[1][0] * (m[0][1] * m[2][2] - m[2][1] * m[0][2])
                        + m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]));

                Matrix<3, 3, T> inverse;
                inverse[0][0] =  (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * oneOverDeterminant;
                inverse[1][0] = -(m[1][0] * m[2][2] - m[2][0] * m[1][2]) * oneOverDeterminant;
                inverse[2][0] =  (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * oneOverDeterminant;
                inverse[0][1] = -(m[0][1] * m[2][2] - m[2][1] * m[0][2]) * oneOverDeterminant;
                inverse[1][1] =  (m[0][0] * m[2][2] - m[2][0] * m[0][2]) * oneOverDeterminant;
                inverse[2][1] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]) * oneOverDeterminant;
                inverse[0][2] =  (m[0][1] * m[1][2] - m[1][1] * m[0][2]) * oneOverDeterminant;
                inverse[1][2] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]) * oneOverDeterminant;
                inverse[2][2] =  (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * oneOverDeterminant;

                return inverse;
            }

            inline constexpr static Matrix<4,4,T> Get(const Matrix<4,4,T>& m) noexcept
            {
                T coef00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
                T coef02 = m[1][2] * m[3][3] - m[3][2] * m[1][3];
                T coef03 = m[1][2] * m[2][3] - m[2][2] * m[1][3];

                T coef04 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
                T coef06 = m[1][1] * m[3][3] - m[3][1] * m[1][3];
                T coef07 = m[1][1] * m[2][3] - m[2][1] * m[1][3];

                T coef08 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
                T coef10 = m[1][1] * m[3][2] - m[3][1] * m[1][2];
                T coef11 = m[1][1] * m[2][2] - m[2][1] * m[1][2];

                T coef12 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
                T coef14 = m[1][0] * m[3][3] - m[3][0] * m[1][3];
                T coef15 = m[1][0] * m[2][3] - m[2][0] * m[1][3];

                T coef16 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
                T coef18 = m[1][0] * m[3][2] - m[3][0] * m[1][2];
                T coef19 = m[1][0] * m[2][2] - m[2][0] * m[1][2];

                T coef20 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
                T coef22 = m[1][0] * m[3][1] - m[3][0] * m[1][1];
                T coef23 = m[1][0] * m[2][1] - m[2][0] * m[1][1];

                Vector<4, T> fac0(coef00, coef00, coef02, coef03);
                Vector<4, T> fac1(coef04, coef04, coef06, coef07);
                Vector<4, T> fac2(coef08, coef08, coef10, coef11);
                Vector<4, T> fac3(coef12, coef12, coef14, coef15);
                Vector<4, T> fac4(coef16, coef16, coef18, coef19);
                Vector<4, T> fac5(coef20, coef20, coef22, coef23);

                Vector<4, T> vec0(m[1][0], m[0][0], m[0][0], m[0][0]);
                Vector<4, T> vec1(m[1][1], m[0][1], m[0][1], m[0][1]);
                Vector<4, T> vec2(m[1][2], m[0][2], m[0][2], m[0][2]);
                Vector<4, T> vec3(m[1][3], m[0][3], m[0][3], m[0][3]);

                Vector<4, T> inv0(vec1 * fac0 - vec2 * fac1 + vec3 * fac2);
                Vector<4, T> inv1(vec0 * fac0 - vec2 * fac3 + vec3 * fac4);
                Vector<4, T> inv2(vec0 * fac1 - vec1 * fac3 + vec3 * fac5);
                Vector<4, T> inv3(vec0 * fac2 - vec1 * fac4 + vec2 * fac5);

                Vector<4, T> SignA(
                        static_cast<T>( 1),
                        static_cast<T>(-1),
                        static_cast<T>( 1),
                        static_cast<T>(-1)
                );

                Vector<4, T> SignB(
                        static_cast<T>(-1),
                        static_cast<T>( 1),
                        static_cast<T>(-1),
                        static_cast<T>( 1)
                );

                Matrix<4, 4, T> Inverse(inv0 * SignA, inv1 * SignB, inv2 * SignA, inv3 * SignB);

                Vector<4, T> Row0(Inverse[0][0], Inverse[1][0], Inverse[2][0], Inverse[3][0]);

                Vector<4, T> Dot0(m[0] * Row0);
                T Dot1 = (Dot0[0] + Dot0[1]) + (Dot0[2] + Dot0[3]);

                T OneOverDeterminant = static_cast<T>(1) / Dot1;

                return Inverse * OneOverDeterminant;
            }
        };

        constexpr static unsigned int rows = R;
        constexpr static unsigned int columns = C;

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