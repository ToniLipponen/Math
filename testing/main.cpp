#include "../Math.h"
#include <iostream>

int main()
{
    tml::Matrix4f scaleMat = tml::Matrix4f::Scale(tml::Vector4f(2.f,3.f,2.f,1.f));

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            std::cout << scaleMat[i][j] << " ";
        }

        std::cout << "\n";
    }

    std::cout << "\n";

    tml::Matrix4f translateMat = tml::Matrix4f::Translate(tml::Vector4f(2.f,3.f,2.f,1.f));

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            std::cout << translateMat[i][j] << " ";
        }

        std::cout << "\n";
    }

    std::cout << "\n";

    return 0;
}