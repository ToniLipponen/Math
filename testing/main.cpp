#include "../Math.h"
#include <iostream>
#include <cmath>

namespace tml
{
    using Vector2f = tml::Vector<2, float>;
    using Vector3f = tml::Vector<3, float>;
    using Vector4f = tml::Vector<4, float>;
}

int main()
{
    tml::Vector4f vec(1.0f, 2.0f, 3.0f, 4.0f);
    tml::Vector4f vec2(5.0f, 6.0f, 7.0f, 8.0f);
    auto result = vec.Cross(vec2);

    for(int i = 0; i < 4; i++)
    {
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;
}