#pragma once

#include "EFM/detail/check_root.h++"

#ifndef EFM_FieldMap3D_d257e2ce3a0979995ab8973a34ffab9a11b34d4356af2a193de7791a6a6121b7
#define EFM_FieldMap3D_d257e2ce3a0979995ab8973a34ffab9a11b34d4356af2a193de7791a6a6121b7

#include "EFM/FieldType.h++"
#include "EFM/detail/mucify.h++"

#include "TFile.h"
#include "TNtuple.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#if __cplusplus >= 202002L
#include <span>
#endif

namespace EFM {

/// @brief Interpolate field on a 3-dimensional regular grid. Trilinear
/// interpolation is performed. The coordinates is projected by `Proj`, this can
/// be used to move the field or do symmetry operation.
/// @tparam Proj projection to apply to the coordinates (default: identity)
/// @tparam T field value type
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<typename T, typename Proj = detail::multidentity,
         typename Coord = double, typename Allocator = std::allocator<T>,
         std::enable_if_t<detail::is_general_arithmetic_v<T>, bool> = true,
         std::enable_if_t<FieldTypeQualified<T>::value, bool> = true,
         std::enable_if_t<std::is_floating_point_v<Coord>, bool> = true>
class FieldMap3D : private Proj {
public:
    using ValueType = T;
    using CoordinateType = Coord;

public:
    FieldMap3D(std::string_view fileName, std::string_view nTupleName,
               Proj proj = {}, Coord tolerance = 0.001,
               const Allocator& allocator = {}) :
        Proj{std::move(proj)},
        fFieldDimension{},
        fGrid{},
        fField{allocator} {
        TFile file{std::string{fileName}.c_str()};
        if (not file.IsOpen()) {
            std::stringstream err;
            err << "EFM::FieldMap3D: Cannot open file "
                << std::quoted(fileName);
            throw std::runtime_error{err.str()};
        }
        const auto nTuple{file.Get<TNtuple>(std::string{nTupleName}.c_str())};
        if (nTuple == nullptr) {
            std::stringstream err;
            err << "EFM::FieldMap3D: TNtuple " << std::quoted(nTupleName)
                << " does not exist in file " << std::quoted(fileName);
            throw std::runtime_error{err.str()};
        }
        Initialize(nTuple, tolerance);
    }

    explicit FieldMap3D(TNtuple* nTuple, Proj proj = {},
                        Coord tolerance = 0.001,
                        const Allocator& allocator = {}) :
        Proj{std::move(proj)},
        fFieldDimension{},
        fGrid{},
        fField{allocator} {
        Initialize(nTuple, tolerance);
    }

    auto Projection(Proj proj) -> void {
        static_cast<Proj&>(*this) = std::move(proj);
    }

    auto operator()(Coord x, Coord y, Coord z) const -> T {
        std::tie(x, y, z) = static_cast<const Proj&>(*this)(x, y, z);

        x = std::clamp(x, get<0>(fGrid).min, get<0>(fGrid).max);
        y = std::clamp(y, get<1>(fGrid).min, get<1>(fGrid).max);
        z = std::clamp(z, get<2>(fGrid).min, get<2>(fGrid).max);

        const auto u{(x - get<0>(fGrid).min) / get<0>(fGrid).delta};
        const auto v{(y - get<1>(fGrid).min) / get<1>(fGrid).delta};
        const auto w{(z - get<2>(fGrid).min) / get<2>(fGrid).delta};

        const auto i{static_cast<int>(u)};
        const auto j{static_cast<int>(v)};
        const auto k{static_cast<int>(w)};
        // clang-format off
        return detail::trilerp(Field(i,     j,     k    ),
                               Field(i,     j,     k + 1),
                               Field(i,     j + 1, k    ),
                               Field(i,     j + 1, k + 1),
                               Field(i + 1, j,     k    ),
                               Field(i + 1, j,     k + 1),
                               Field(i + 1, j + 1, k    ),
                               Field(i + 1, j + 1, k + 1),
                               u - i, v - j, w - k); // clang-format on
    }

    auto operator()(std::array<Coord, 3> x) const -> T {
        return (*this)(x[0], x[1], x[2]);
    }

#if __cplusplus >= 202002L
    auto operator()(std::span<Coord, 3> x) const -> T {
        return (*this)(x[0], x[1], x[2]);
    }
#endif

    template<typename Char>
    friend auto operator<<(std::basic_ostream<Char>& os,
                           const FieldMap3D& self) -> decltype(os) {
        const auto& x{get<0>(self.fGrid)};
        const auto& y{get<1>(self.fGrid)};
        const auto& z{get<2>(self.fGrid)};
        for (int i{}; i < ssize(self.fField); ++i) {
            const auto x0{x.min + i / (y.n * z.n) % x.n * x.delta};
            const auto y0{y.min + i / z.n % y.n * y.delta};
            const auto z0{z.min + i % z.n * z.delta};
            os << x0 << '\t' << y0 << '\t' << z0;
            const auto& field{self.fField[i]};
            if constexpr (FieldType<T>::value == RealScalarField) {
                os << '\t' << field;
            } else if constexpr (FieldType<T>::value == ComplexScalarField) {
                os << '\t' << real(field) << '\t' << imag(field);
            } else if constexpr (FieldType<T>::value == RealVectorField) {
                for (int i{}; i < self.fFieldDimension; ++i) {
                    os << '\t' << field[i];
                }
            } else if constexpr (FieldType<T>::value == ComplexVectorField) {
                assert(self.fFieldDimension % 2 == 0);
                for (int i{}; i < self.fFieldDimension / 2; ++i) {
                    os << '\t' << real(field[i]) << '\t' << imag(field[i]);
                }
            }
            os << '\n';
        }
        return os;
    }

private:
    auto Initialize(TNtuple* nTuple, Coord tolerance) {
        if (nTuple == nullptr) {
            throw std::runtime_error{"EFM::FieldMap3D: nTuple == nullptr"};
        }

        const auto nColumn{nTuple->GetNvar()};
        if (nColumn <= 3) {
            std::stringstream err;
            err << "EFM::FieldMap3D: "
                   "There must be >= 4 columns in field data, but "
                << std::quoted(nTuple->GetName()) << " has only " << nColumn
                << " columns";
            throw std::runtime_error{err.str()};
        }
        fFieldDimension = nColumn - 3;
        if (FieldType<T>::value == RealScalarField and fFieldDimension > 1) {
            throw std::runtime_error{
                "EFM::FieldMap3D: Data looks like vector or complex "
                "field but this is a real scalar field"};
        }
        if (FieldType<T>::value == ComplexScalarField and
            fFieldDimension != 2) {
            throw std::runtime_error{
                "EFM::FieldMap3D: this is a complex scalar field but data does "
                "not look like one"};
        }
        if (FieldType<T>::value == ComplexVectorField and
            fFieldDimension % 2 != 0) {
            throw std::runtime_error{
                "EFM::FieldMap3D: this is a complex vector field but data does "
                "not look like one"};
        }

        std::map<std::array<Coord, 3>, T> gridData;
        for (int i{}; nTuple->GetEntry(i); ++i) {
            const auto data{nTuple->GetArgs()}; // clang-format off
            const auto [it, inserted]{
                gridData.insert({{data[0], data[1], data[2]}, T{}})}; // clang-format on
            if (not inserted) {
                std::stringstream err;
                err << "EFM::FieldMap3D: "
                       "point ["
                    << data[0] << ", " << data[1] << ", " << data[2]
                    << "] "
                       "appeared twice";
                throw std::runtime_error{err.str()};
            }
            auto& field{it->second};
            if constexpr (FieldType<T>::value == RealScalarField) {
                field = data[3];
            } else if constexpr (FieldType<T>::value == ComplexScalarField) {
                field = {data[3], data[4]};
            } else if constexpr (FieldType<T>::value == RealVectorField) {
                for (int i{}; i < fFieldDimension; ++i) {
                    field[i] = data[i + 3];
                }
            } else if constexpr (FieldType<T>::value == ComplexVectorField) {
                assert(fFieldDimension % 2 == 0);
                for (int i{}; i < fFieldDimension / 2; ++i) {
                    field[i] = {data[2 * i], data[2 * i + 1]};
                }
            }
            if (i == std::numeric_limits<int>::max()) {
                std::runtime_error{"EFM::FieldMap3D: Too much grid points "
                                   "(should < INT_MAX in total)"};
            }
        }

        std::array<Coord, 3> x0;
        x0.fill(std::numeric_limits<Coord>::lowest());
        std::array<Coord, 3> lastDelta{};
        std::array<int, 3> counter{};
        std::array<int, 3> countChecker{};
        fField.reserve(gridData.size());
        for (auto&& [x, field] : std::as_const(gridData)) {
            for (int i{}; i < 3; ++i) {
                if (x[i] > x0[i]) {
                    // check normal delta
                    if (fGrid[i].n >= 2) {
                        if (std::abs((x[i] - x0[i]) / lastDelta[i] - 1) >
                            tolerance) {
                            throw std::runtime_error{
                                "EFM::FieldMap3D: Irregular grid "
                                "(inconsistent delta)"};
                        }
                    }
                    lastDelta[i] = x[i] - x0[i];
                    // check normal switch count
                    if (fGrid[i].n <= 1) {
                        countChecker[i] = counter[i];
                    } else if (counter[i] != countChecker[i]) {
                        throw std::runtime_error{
                            "EFM::FieldMap3D: Irregular grid"};
                    }
                    counter[i] = 0;
                    // update
                    ++fGrid[i].n;
                    x0[i] = x[i];
                }
                ++counter[i];
            }
            fField.emplace_back(field);
        }
        fField.shrink_to_fit();

        if (ssize(fField) !=
            get<0>(fGrid).n * get<1>(fGrid).n * get<2>(fGrid).n) {
            throw std::runtime_error{
                "EFM::FieldMap3D: Irregular grid (N != Nx * Ny * Nz)"};
        }

        const auto& first{gridData.cbegin()->first};
        const auto& last{gridData.crbegin()->first};
        for (int i{}; i < 3; ++i) {
            fGrid[i].min = first[i];
            fGrid[i].max = last[i];
        }
        for (auto&& grid : fGrid) {
            if (grid.n < 2) {
                throw std::runtime_error{"EFM::FieldMap3D: Too few grid points "
                                         "(should >= 2 in each direction)"};
            }
            if (grid.n >= 1 / tolerance) {
                std::stringstream err;
                err << "EFM::FieldMap3D: Too much grid points (in each "
                       "direction should < 1 / tolerance, tolerance == "
                    << tolerance << ")";
                throw std::runtime_error{err.str()};
            }
            grid.delta = (grid.max - grid.min) / (grid.n - 1);
        }
    }

    auto Field(int i, int j, int k) const -> const T& {
        assert(i < get<0>(fGrid).n);
        assert(j < get<1>(fGrid).n);
        assert(k < get<2>(fGrid).n);
        return fField[(i * get<1>(fGrid).n + j) * get<2>(fGrid).n + k];
    }

private:
    struct GridInfomation {
        int n;
        Coord min;
        Coord max;
        Coord delta;
    };

private:
    int fFieldDimension;
    std::array<GridInfomation, 3> fGrid;
    std::vector<T, Allocator> fField;
};

/// @brief Mirror symmetry operation along x-axis. Flip the x-coordinate from
/// negative to positive.
struct SymmetryX {
    [[nodiscard]] constexpr auto
    operator()(double x, double y,
               double z) const noexcept -> std::array<double, 3> {
        return {std::abs(x), y, z};
    }
};

/// @brief A YZ plane mirror symmetric field.
/// @tparam T field value type
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<typename T, typename Coord = double,
         typename Allocator = std::allocator<T>>
using FieldMap3DSymX = FieldMap3D<T, SymmetryX, Coord, Allocator>;

/// @brief Mirror symmetry operation along y-axis. Flip the y-coordinate from
/// negative to positive.
struct SymmetryY {
    [[nodiscard]] constexpr auto
    operator()(double x, double y,
               double z) const noexcept -> std::array<double, 3> {
        return {x, std::abs(y), z};
    }
};

/// @brief A ZX plane mirror symmetric field.
/// @tparam T field value type
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<typename T, typename Coord = double,
         typename Allocator = std::allocator<T>>
using FieldMap3DSymY = FieldMap3D<T, SymmetryY, Coord, Allocator>;

/// @brief Mirror symmetry operation along z-axis. Flip the z-coordinate from
/// negative to positive.
struct SymmetryZ {
    [[nodiscard]] constexpr auto
    operator()(double x, double y,
               double z) const noexcept -> std::array<double, 3> {
        return {x, y, std::abs(z)};
    }
};

/// @brief A XY plane mirror symmetric field.
/// @tparam T field value type
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<typename T, typename Coord = double,
         typename Allocator = std::allocator<T>>
using FieldMap3DSymZ = FieldMap3D<T, SymmetryZ, Coord, Allocator>;

/// @brief Mirror symmetry operation along x-axis and y-axis. Flip the
/// x-coordinate and y-coordinate from negative to positive.
struct SymmetryXY {
    [[nodiscard]] constexpr auto
    operator()(double x, double y,
               double z) const noexcept -> std::array<double, 3> {
        return {std::abs(x), std::abs(y), z};
    }
};

/// @brief A YZ and ZX plane mirror symmetric field.
/// @tparam T field value type
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<typename T, typename Coord = double,
         typename Allocator = std::allocator<T>>
using FieldMap3DSymXY = FieldMap3D<T, SymmetryXY, Coord, Allocator>;

/// @brief Mirror symmetry operation along x-axis and z-axis. Flip the
/// x-coordinate and z-coordinate from negative to positive.
struct SymmetryXZ {
    [[nodiscard]] constexpr auto
    operator()(double x, double y,
               double z) const noexcept -> std::array<double, 3> {
        return {std::abs(x), y, std::abs(z)};
    }
};

/// @brief A XY and YZ plane mirror symmetric field.
/// @tparam T field value type
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<typename T, typename Coord = double,
         typename Allocator = std::allocator<T>>
using FieldMap3DSymXZ = FieldMap3D<T, SymmetryXZ, Coord, Allocator>;

/// @brief Mirror symmetry operation along y-axis and z-axis. Flip the
/// y-coordinate and z-coordinate from negative to positive.
struct SymmetryYZ {
    [[nodiscard]] constexpr auto
    operator()(double x, double y,
               double z) const noexcept -> std::array<double, 3> {
        return {x, std::abs(y), std::abs(z)};
    }
};

/// @brief A ZX and XY plane mirror symmetric field.
/// @tparam T field value type
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<typename T, typename Coord = double,
         typename Allocator = std::allocator<T>>
using FieldMap3DSymYZ = FieldMap3D<T, SymmetryYZ, Coord, Allocator>;

/// @brief Mirror symmetry operation along x-axis, y-axis and z-axis. Flip the
/// x-coordinate, y-coordinate and z-coordinate from negative to positive.
struct SymmetryXYZ {
    [[nodiscard]] constexpr auto
    operator()(double x, double y,
               double z) const noexcept -> std::array<double, 3> {
        return {std::abs(x), std::abs(y), std::abs(z)};
    }
};

/// @brief A XY, YZ and ZX plane mirror symmetric field.
/// @tparam T field value type
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<typename T, typename Coord = double,
         typename Allocator = std::allocator<T>>
using FieldMap3DSymXYZ = FieldMap3D<T, SymmetryXYZ, Coord, Allocator>;

} // namespace EFM

#endif
