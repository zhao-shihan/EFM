#pragma once

#include <type_traits>
#include <utility>

namespace EFM::detail {

template<typename, typename = void>
struct has_subscript : std::false_type {};

template<typename T>
struct has_subscript<T, std::void_t<decltype(std::declval<T&>()[0]),
                                    decltype(std::declval<const T&>()[0])>> :
    std::true_type {};

template<typename T>
inline constexpr bool has_subscript_v{has_subscript<T>::value};

} // namespace EFM::detail
