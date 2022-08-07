#include "../Math.h"
#include <iostream>

int main()
{
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