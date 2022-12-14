#define TML_USE_DEGREES
#include "../Math.h"
#include <iostream>

template<unsigned int R, unsigned int C, typename T>
void PrintMatrix(const tml::Matrix<R,C,T>& m)
{
    for(int i = 0; i < m.rows; i++)
    {
        for(int j = 0; j < m.columns; j++)
        {
            std::cout << m[i][j] << "  ";
        }

        std::cout << "\n";
    }
    std::cout << "\n";
}

template<unsigned int N, typename T>
void PrintVector(const tml::Vector<N,T>& v)
{
    for(int i = 0; i < v.elements; i++)
        std::cout << v[i] << " ";

    std::cout << "\n\n";
}

int main()
{
    // Vector - matrix multiplication
    auto matrix = tml::Matrix4x4f::Rotate(tml::Vector3f(0.f,0.f,1.f), 90);
    tml::Vector4f result = matrix * tml::Vector4f(10.f, 0.f, 0.f, 0.f);

    // Matrix - matrix multiplication
    tml::Matrix<3,2,float> m1;
    m1[0] = tml::Vector2f{1.f, 2.f};
    m1[1] = tml::Vector2f{3.f, 4.f};
    m1[2] = tml::Vector2f{5.f, 1.f};

    tml::Matrix<2,1,float> m2;
    m2[0] = tml::Vector<1,float>{2.f};
    m2[1] = tml::Vector<1,float>{4.f};

    auto m3 = m1 * m2;

    PrintMatrix(matrix);
    PrintVector(result);
    PrintMatrix(m3);

    return 0;
}