#pragma once

#include <complex>

namespace EFM::detail {

template<typename, typename = void>
struct is_std_complex : std::false_type {};

template<typename T>
struct is_std_complex<std::complex<T>,
                      std::enable_if_t<std::is_arithmetic_v<T>>> :
    std::true_type {};

template<typename T>
inline constexpr bool is_std_complex_v{is_std_complex<T>::value};

} // namespace EFM::detail
