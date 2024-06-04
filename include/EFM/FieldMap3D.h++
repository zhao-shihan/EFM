#pragma once

#include "EFM/detail/has_subscript.h++"
#include "EFM/detail/is_std_complex.h++"
#include "EFM/detail/mucify.h++"
#include "TFile.h"
#include "TNtuple.h"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <limits>
#include <map>
#include <memory_resource>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#if __cplusplus >= 202002L
#include <span>
#endif

namespace EFM {

/// @brief Interpolate field on a 3-dimensional regular grid. Trilinear
/// interpolation is performed.
/// @tparam T field value type.
/// @tparam Coord grid value type (floating-point, default: double)
/// @tparam Allocator the allocator type used by internal field point data
/// member (default: std::allocator<T>).
template<
    typename T, typename Coord = double, typename Allocator = std::allocator<T>,
    std::enable_if_t<detail::is_general_arithmetic_v<T>, bool> = true,
    std::enable_if_t<
        std::is_arithmetic_v<T> or             // real scalar field
            detail::is_std_complex_v<T> or     // complex scalar field
            std::is_arithmetic_v<std::decay_t< // real vector field
                decltype(std::declval<std::conditional_t<
                             detail::has_subscript_v<T>, T&, char*>>()[0])>> or
            detail::is_std_complex_v<std::decay_t< // complex vector field
                decltype(std::declval<std::conditional_t<
                             detail::has_subscript_v<T>, T&, char*>>()[0])>>,
        bool> = true,
    std::enable_if_t<std::is_floating_point_v<Coord>, bool> = true>
class FieldMap3D {
public:
    FieldMap3D(std::string_view fileName, std::string_view nTupleName,
               Coord tolerance = 10 * std::numeric_limits<float>::epsilon(),
               const Allocator& allocator = {}) :
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

    explicit FieldMap3D(TNtuple* nTuple, // TNtuple stores float
                        Coord tolerance = 10 *
                                          std::numeric_limits<float>::epsilon(),
                        const Allocator& allocator = {}) :
        fGrid{},
        fField{allocator} {
        Initialize(nTuple, tolerance);
    }

    auto operator()(Coord x, Coord y, Coord z) const -> T {
        x = std::clamp(x, fGrid[0].min, fGrid[0].max);
        y = std::clamp(y, fGrid[1].min, fGrid[1].max);
        z = std::clamp(z, fGrid[2].min, fGrid[2].max);

        const auto u{(x - fGrid[0].min) / fGrid[0].delta};
        const auto v{(y - fGrid[1].min) / fGrid[1].delta};
        const auto w{(z - fGrid[2].min) / fGrid[2].delta};

        const auto i{static_cast<std::size_t>(u)};
        const auto j{static_cast<std::size_t>(v)};
        const auto k{static_cast<std::size_t>(w)};
        // clang-format off
        return detail::trilerp({fField[Idx(i,     j,     k    )],
                                fField[Idx(i + 1, j,     k    )],
                                fField[Idx(i,     j + 1, k    )],
                                fField[Idx(i + 1, j + 1, k    )],
                                fField[Idx(i,     j,     k + 1)],
                                fField[Idx(i + 1, j,     k + 1)],
                                fField[Idx(i,     j + 1, k + 1)],
                                fField[Idx(i + 1, j + 1, k + 1)]},
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
        const auto fieldDimension{nColumn - 3};
        if (std::is_arithmetic_v<T> and fieldDimension > 1) {
            throw std::runtime_error{
                "EFM::FieldMap3D: Data looks like vector or complex "
                "field but this is a real scalar field"};
        }
        if (detail::is_std_complex_v<T> and fieldDimension != 2) {
            throw std::runtime_error{
                "EFM::FieldMap3D: this is a complex scalar field but data does "
                "not look like one"};
        }
        if (detail::is_std_complex_v<std::decay_t<
                decltype(std::declval<std::conditional_t<
                             detail::has_subscript_v<T>, T&, char*>>()[0])>> and
            fieldDimension % 2 != 0) {
            throw std::runtime_error{
                "EFM::FieldMap3D: this is a complex vector field but data does "
                "not look like one"};
        }

        std::map<std::array<Coord, 3>, T> gridData;
        for (std::size_t i{}; nTuple->GetEntry(i); ++i) {
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
            if constexpr (std::is_arithmetic_v<T>) {
                // real scalar field
                field = data[3];
            } else if constexpr (detail::is_std_complex_v<T>) {
                // complex scalar field
                field = {data[3], data[4]};
            } else if constexpr (std::is_arithmetic_v<std::decay_t<
                                     decltype(std::declval<std::conditional_t<
                                                  detail::has_subscript_v<T>,
                                                  T&, char*>>()[0])>>) {
                // real vector field
                for (int i{}; i < fieldDimension; ++i) {
                    field[i] = data[i + 3];
                }
            } else if constexpr (detail::is_std_complex_v<std::decay_t<
                                     decltype(std::declval<std::conditional_t<
                                                  detail::has_subscript_v<T>,
                                                  T&, char*>>()[0])>>) {
                // complex vector field
                for (int i{}; i < fieldDimension / 2; ++i) {
                    field[i] = {data[2 * i], data[2 * i + 1]};
                }
            }
        }

        fField.reserve(gridData.size());

        std::array<Coord, 3> x0;
        x0.fill(std::numeric_limits<Coord>::lowest());
        std::array<Coord, 3> lastDelta{};
        std::array<std::size_t, 3> counter{};
        std::array<std::size_t, 3> countChecker{};
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

        if (fField.size() != fGrid[0].n * fGrid[1].n * fGrid[2].n) {
            throw std::runtime_error{
                "EFM::FieldMap3D: Irregular grid (N != Nx * Ny * Nz)"};
        }

        const auto& first{gridData.cbegin()->first};
        const auto& last{gridData.crbegin()->first};
        for (int i{}; i < 3; ++i) {
            fGrid[i].min = first[i];
            fGrid[i].max = last[i];
            const auto length{fGrid[i].max - fGrid[i].min};
            if (length < tolerance) {
                std::stringstream err;
                err << "EFM::FieldMap3D: Grid too small (max - min == "
                    << length << " < tolerance == " << tolerance << ')';
                throw std::runtime_error{err.str()};
            }
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

    /// @brief Get nomalized index for fField.
    /// @param i original x index.
    /// @param j original y index.
    /// @param k original z index.
    /// @return nomalized index.
    auto Idx(std::size_t i, std::size_t j, std::size_t k) const -> std::size_t {
        assert(i < fGrid[0].n);
        assert(j < fGrid[1].n);
        assert(k < fGrid[2].n);
        return i + fGrid[1].n * (j + fGrid[2].n * k);
    }

private:
    struct GridInfomation {
        std::size_t n;
        Coord min;
        Coord max;
        Coord delta;
    };

private:
    std::array<GridInfomation, 3> fGrid;
    std::vector<T, Allocator> fField;
};

} // namespace EFM
