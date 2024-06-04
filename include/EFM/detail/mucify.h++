#pragma once

#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

namespace EFM {

namespace detail::impl {

template<typename, typename, typename = void>
struct has_binary_arithmetic_with : std::false_type {};

template<typename T, typename U>
struct has_binary_arithmetic_with<
    T, U,
    std::void_t<
        // T +/- T
        decltype(std::declval<const T&>() + std::declval<const T&>()),
        decltype(std::declval<const T&>() - std::declval<const T&>()),
        // k * T
        decltype(+std::declval<const T&>()),
        decltype(-std::declval<const T&>()),
        decltype(1 * std::declval<const T&>()),
        decltype(std::declval<const T&>() * 1),
        decltype(std::declval<const T&>() / 1),
        // U +/- U
        decltype(std::declval<const U&>() + std::declval<const U&>()),
        decltype(std::declval<const U&>() - std::declval<const U&>()),
        // k * U
        decltype(+std::declval<const U&>()),
        decltype(-std::declval<const U&>()),
        decltype(1 * std::declval<const U&>()),
        decltype(std::declval<const U&>() * 1),
        decltype(std::declval<const U&>() / 1),
        // T +/- U
        decltype(std::declval<const T&>() + std::declval<const U&>()),
        decltype(std::declval<const T&>() - std::declval<const U&>()),
        // U +/- T
        decltype(std::declval<const U&>() + std::declval<const T&>()),
        decltype(std::declval<const U&>() - std::declval<const T&>()),
        // T x= X
        decltype(std::declval<T&>() += std::declval<const U&>()),
        decltype(std::declval<T&>() -= std::declval<const U&>()),
        decltype(std::declval<T&>() *= 1), decltype(std::declval<T&>() /= 1)>> :
    std::true_type {};

template<typename T, typename U>
inline constexpr bool has_binary_arithmetic_with_v{
    has_binary_arithmetic_with<T, U>::value};

} // namespace detail::impl

namespace detail {

template<typename T>
struct is_general_arithmetic :
    std::bool_constant<std::is_default_constructible_v<T> and
                       impl::has_binary_arithmetic_with_v<T, T> and
                       impl::has_binary_arithmetic_with_v<
                           T, decltype(std::declval<const T&>() +
                                       std::declval<const T&>())> and
                       impl::has_binary_arithmetic_with_v<
                           T, decltype(std::declval<const T&>() -
                                       std::declval<const T&>())> and
                       impl::has_binary_arithmetic_with_v<
                           T, decltype(1 * std::declval<const T&>())> and
                       impl::has_binary_arithmetic_with_v<
                           T, decltype(std::declval<const T&>() * 1)> and
                       impl::has_binary_arithmetic_with_v<
                           T, decltype(std::declval<const T&>() / 1)>> {};

template<typename T>
inline constexpr bool is_general_arithmetic_v{is_general_arithmetic<T>::value};

} // namespace detail

namespace detail {

/// @brief Performs bilinear interpolation.
/// @tparam T value type, can be a scalar or vector or something.
/// @tparam U scalar type
/// @param c values on square grid. See note for details.
/// @param u interpolation parameter. 0<u<1 implies interpolation, otherwise
/// extrapolation.
/// @param v interpolation parameter. 0<v<1 implies interpolation, otherwise
/// extrapolation.
/// @return interpolated value.
/// @note For parameter c:
///
///      c[2]      c[3]
///          +----+
///          |    |
/// v        +----+
/// ^    c[0]      c[1]
/// |
/// +----> u
template<typename T, typename U,
         std::enable_if_t<is_general_arithmetic_v<U>, bool> = true,
         std::enable_if_t<std::is_floating_point_v<U>, bool> = true>
constexpr auto bilerp(const std::array<T, 4>& c, U u, U v) -> T {
    const T a00{c[0]};
    const T a10{c[1] - a00};
    const T a01{c[2] - a00};
    const T a11{c[3] - c[2] - a10};
    return a00 + a10 * u + (a01 + a11 * u) * v;
}

} // namespace detail

namespace detail {

/// @brief Performs bilinear interpolation.
/// @tparam T value type, can be a scalar or vector or something.
/// @tparam U scalar type
/// @param c values on square grid. See note for details.
/// @param u interpolation parameter. 0<u<1 implies interpolation, otherwise
/// extrapolation.
/// @param v interpolation parameter. 0<v<1 implies interpolation, otherwise
/// extrapolation.
/// @return interpolated value.
/// @note For parameter c:
///
///         c[6]    c[7]
///            +----+
///      c[2] /|   /|
///          +----+ + c[5]
///          |c[4]|/
/// v w      +----+
/// ^ ^  c[0]     c[1]
/// |/
/// +----> u
template<typename T, typename U,
         std::enable_if_t<is_general_arithmetic_v<U>, bool> = true,
         std::enable_if_t<std::is_floating_point_v<U>, bool> = true>
constexpr auto trilerp(const std::array<T, 8>& c, U u, U v, U w) -> T {
    const T a000{c[0]};
    const T a100{c[1] - a000};
    const T a010{c[2] - a000};
    const T a001{c[4] - a000};
    const T b101{c[5] - c[4]};
    const T a110{c[3] - c[1] - a010};
    const T a101{b101 - a100};
    const T a011{c[5] - c[1] - a001};
    const T a111{c[7] - c[6] - b101 - a110};
    return a000 + (a100 + a101 * w + (a110 + a111 * w) * v) * u + a001 * w +
           (a010 + a011 * w) * v;
}

} // namespace detail

} // namespace EFM
