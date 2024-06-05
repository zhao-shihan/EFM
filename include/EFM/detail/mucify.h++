#pragma once

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
        /* decltype(1 * std::declval<const T&>()),
        decltype(std::declval<const T&>() * 1),
        decltype(std::declval<const T&>() / 1), */
        // U +/- U
        decltype(std::declval<const U&>() + std::declval<const U&>()),
        decltype(std::declval<const U&>() - std::declval<const U&>()),
        // k * U
        decltype(+std::declval<const U&>()),
        decltype(-std::declval<const U&>()),
        /* decltype(1 * std::declval<const U&>()),
        decltype(std::declval<const U&>() * 1),
        decltype(std::declval<const U&>() / 1), */
        // T +/- U
        decltype(std::declval<const T&>() + std::declval<const U&>()),
        decltype(std::declval<const T&>() - std::declval<const U&>()),
        // U +/- T
        decltype(std::declval<const U&>() + std::declval<const T&>()),
        decltype(std::declval<const U&>() - std::declval<const T&>()),
        // T x= X
        decltype(std::declval<T&>() += std::declval<const U&>()),
        decltype(std::declval<T&>() -= std::declval<const U&>())//,
        /* decltype(std::declval<T&>() *= 1), decltype(std::declval<T&>() /= 1) */>> :
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
                                       std::declval<const T&>())> /* and
                       impl::has_binary_arithmetic_with_v<
                           T, decltype(1 * std::declval<const T&>())> and
                       impl::has_binary_arithmetic_with_v<
                           T, decltype(std::declval<const T&>() * 1)> and
                       impl::has_binary_arithmetic_with_v<
                           T, decltype(std::declval<const T&>() / 1)> */> {};

template<typename T>
inline constexpr bool is_general_arithmetic_v{is_general_arithmetic<T>::value};

} // namespace detail

namespace detail {

/// @brief Performs linear interpolation.
/// @tparam T value type, can be a scalar or vector or something.
/// @tparam U scalar type
/// @param c values on endpoints
/// @param t interpolation parameter. 0<t<1 implies interpolation, otherwise
/// extrapolation.
/// @return interpolated value.
template<typename T, typename U,
         std::enable_if_t<is_general_arithmetic_v<T>, bool> = true,
         std::enable_if_t<std::is_floating_point_v<U>, bool> = true>
constexpr auto lerp(const T& a, const T& b, U t) -> T {
    if constexpr (std::is_integral_v<T>) {
        return (1 - t) * a + t * b;
    } else {
        return a + t * (b - a);
    }
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
///      c01      c11
///         +----+
///         |    |
/// v       +----+
/// ^    c00      c10
/// |
/// +----> u
template<typename T, typename U,
         std::enable_if_t<is_general_arithmetic_v<T>, bool> = true,
         std::enable_if_t<std::is_floating_point_v<U>, bool> = true>
constexpr auto bilerp(const T& c00, const T& c01, const T& c10, const T& c11,
                      U u, U v) -> T {
    return detail::lerp(detail::lerp(c00, c01, v), detail::lerp(c10, c11, v),
                        u);
}

} // namespace detail

namespace detail {

/// @brief Performs trilinear interpolation.
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
///         c011    c111
///            +----+
///      c010 /|   /|
///          +----+ + c101
///          |c001|/
/// v w      +----+
/// ^ ^  c000     c100
/// |/
/// +----> u
template<typename T, typename U,
         std::enable_if_t<is_general_arithmetic_v<T>, bool> = true,
         std::enable_if_t<std::is_floating_point_v<U>, bool> = true>
constexpr auto trilerp(const T& c000, const T& c001, const T& c010,
                       const T& c011, const T& c100, const T& c101,
                       const T& c110, const T& c111, U u, U v, U w) -> T {
    return detail::bilerp(
        detail::lerp(c000, c001, w), detail::lerp(c010, c011, w),
        detail::lerp(c100, c101, w), detail::lerp(c110, c111, w), u, v);
}

} // namespace detail

} // namespace EFM
