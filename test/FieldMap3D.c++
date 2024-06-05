#include "EFM/FieldMap3D.h++"

#include <complex>
#include <iostream>

auto main() -> int {
    EFM::FieldMap3D<double> phi{"data.root", "phi"};
    std::cout << phi << '\n';
    std::cout << phi(1.5, -3.2, 1.3) << '\n';
    std::cout << '\n';
    EFM::FieldMap3D<std::complex<double>> phi1{"data.root", "phi1"};
    std::cout << phi1 << '\n';
    std::cout << phi1(1, 2, 4) << '\n';
}
