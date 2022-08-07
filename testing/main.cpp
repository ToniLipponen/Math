#define TML_USE_DEGREES
#include "../Math.h"
#include <iostream>

int main()
{
    // Vector - matrix multiplication
    auto matrix = tml::Matrix4::Rotate(tml::Vector4(0.f,0.f,1.f,0.f), 90);
    auto result = matrix * tml::Vector4(10.f, 0.f, 0.f, 0.f);

    for(int i = 0; i < result.elements; i++)
        std::cout << static_cast<int>(result[i]) << " ";

    std::cout << "\n\n";

    // Matrix - matrix multiplication
    auto m1 = tml::Matrix<3,2,float>(
            tml::Vector2(1.f, 2.f),
            tml::Vector2(3.f, 4.f),
            tml::Vector2(5.f, 1.f));

    auto m2 = tml::Matrix<2,1,float>(
            tml::Vector<1,float>(2.f),
            tml::Vector<1,float>(4.f));

    auto m3 = m1 * m2;

    for(int i = 0; i < m3.rows; i++)
    {
        for(int j = 0; j < m3.columns; j++)
        {
            std::cout << m3[i][j] << "  ";
        }

        std::cout << "\n";
    }

    return 0;
}