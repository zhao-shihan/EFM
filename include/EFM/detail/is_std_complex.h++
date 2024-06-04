#pragma once

#include <complex>

namespace EFM::detail {

template<typename>
struct is_std_complex : std::false_type {};

template<typename T>
struct is_std_complex<std::complex<T>> : std::true_type {};

template<typename T>
inline constexpr bool is_std_complex_v{is_std_complex<T>::value};

} // namespace EFM::detail
