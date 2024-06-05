#pragma once

#include "EFM/detail/has_subscript.h++"
#include "EFM/detail/is_std_complex.h++"

#include <type_traits>

namespace EFM {

enum FieldEnum {
    RealScalarField,
    ComplexScalarField,
    RealVectorField,
    ComplexVectorField
};

template<typename, typename = void>
struct FieldType {};

template<typename T>
struct FieldType<T, std::enable_if_t<std::is_arithmetic_v<T>>> :
    std::integral_constant<FieldEnum, RealScalarField> {};

template<typename T>
struct FieldType<T, std::enable_if_t<detail::is_std_complex_v<T>>> :
    std::integral_constant<FieldEnum, ComplexScalarField> {};

namespace detail {
struct dummy {};
} // namespace detail

template<typename T>
struct FieldType<
    T,
    std::enable_if_t<std::is_arithmetic_v<std::decay_t<
        decltype(std::declval<std::conditional_t<detail::has_subscript_v<T>, T&,
                                                 detail::dummy*>>()[0])>>>> :
    std::integral_constant<FieldEnum, RealVectorField> {};

template<typename T>
struct FieldType<
    T,
    std::enable_if_t<detail::is_std_complex_v<std::decay_t<
        decltype(std::declval<std::conditional_t<detail::has_subscript_v<T>, T&,
                                                 detail::dummy*>>()[0])>>>> :
    std::integral_constant<FieldEnum, ComplexVectorField> {};

template<typename, typename = void>
struct FieldTypeQualified : std::false_type {};

template<typename T>
struct FieldTypeQualified<T, std::void_t<typename FieldType<T>::value_type>> :
    std::true_type {};

} // namespace EFM
