#pragma once

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "Manager/MaCh3Logger.h"

namespace M3::Utils {
/// @brief Column alignment
enum class Alignment {
  Left,    ///< Align text to the left.
  Right,   ///< Align text to the right.
  Center   ///< Center text within the column.
};

/// @brief Column specification used by TablePrinter
struct ColumnSpec {
  /// Column header displayed in the table.
  std::string name;
  /// Minimum column width.
  std::size_t minWidth = 0;
  /// Maximum column width. Use 0 for no maximum.
  std::size_t maxWidth = 0;
  /// Alignment of values within the column.
  Alignment alignment = Alignment::Left;
};

/// @brief Join all values in a container into a single string.
template <typename Container>
std::string join(const Container& values, std::string_view separator = " ") {
  return fmt::format("{}", fmt::join(values, separator));
}

/// @brief Format a value using a runtime fmt format string.
/// @param value  Value to format.
/// @param format fmt-style format string, for example "{:.2f}".
template <typename T>
static std::string Cell(T&& value, std::string_view format) {
  return fmt::format(fmt::runtime(format), std::forward<T>(value));
}

/// @brief Creates fancy formatted text tables for terminal output.
/// @author Kamil Skwarczynski
class TablePrinter {
 public:
 /// @brief Construct a table printer with the specified columns.
 explicit TablePrinter(std::vector<ColumnSpec> columns) : _columns(std::move(columns)) {}

  /// @brief Add a row to the table.
  template <typename... Args>
  void AddRow(Args&&... args) {
    _rows.push_back({toString(std::forward<Args>(args))...});
  }

  /// @brief Render the complete table, with widths etc. defined by constructor
  [[nodiscard]]
  std::vector<std::string> Render() const
  {
    if (_columns.empty()) return {};

    const auto widths = CalculateWidths();

    std::vector<std::string> result;
    result.reserve(_rows.size() + 4);

    // ┌──────┬────────┐
    result.push_back(MakeBorder(widths, "┌", "┬", "┐", "─"));

    // │ Name │ Value  │
    result.push_back(MakeRow(widths, MakeHeader()));

    // ├──────┼────────┤
    result.push_back(MakeBorder(widths, "├", "┼", "┤", "─"));

    // data rows
    for (const auto& row : _rows)
      result.push_back(MakeRow(widths, row));

    // └──────┴────────┘
    result.push_back(MakeBorder(widths, "└", "┴", "┘", "─"));

    return result;
  }

 private:
  std::vector<ColumnSpec> _columns;
  std::vector<std::vector<std::string>> _rows;
  /// @brief Convert a value to a string using stream formatting.
  template <typename T>
  static std::string toString(T&& value)
  {
    std::ostringstream os;
    os << std::forward<T>(value);
    return os.str();
  }
  /// @brief If name is longer than desired name we truncate for example instead of SuperLongBlarb we have SuperLong...
  static std::string Truncate(const std::string_view value, const std::size_t width)
  {
    if (value.size() <= width) return std::string(value);

    if (width == 0) return {};

    if (width == 1) return "…";

    return std::string(value.substr(0, width - 1)) + "…";
  }
  /// @brief Alling to values to table deepening on values
  static std::string Align(
    const std::string_view value,
    const std::size_t width,
    const Alignment alignment)
  {
    if (value.size() >= width){
      return std::string(value);
    }
    const std::size_t padding = width - value.size();

    switch (alignment) {
      case Alignment::Left:
        return std::string(value) + std::string(padding, ' ');

      case Alignment::Right:
        return std::string(padding, ' ') + std::string(value);

      case Alignment::Center: {
        const std::size_t left = padding / 2;
        const std::size_t right = padding - left;

        return std::string(left, ' ') + std::string(value) + std::string(right, ' ');
      }
    }
    return std::string(value);
  }

  std::vector<std::size_t> CalculateWidths() const
  {
    std::vector<std::size_t> widths;
    widths.reserve(_columns.size());

    for (std::size_t col = 0; col < _columns.size(); ++col) {
      auto width = _columns[col].name.size();

      for (const auto& row : _rows) {
        if (col < row.size()) {
          width = std::max(width, row[col].size());
        }
      }

      width = std::max(width, _columns[col].minWidth);

      if (_columns[col].maxWidth != 0){
        width = std::min(width, _columns[col].maxWidth);
      }
      widths.push_back(width);
    }

    return widths;
  }

  std::vector<std::string> MakeHeader() const
  {
    std::vector<std::string> header(_columns.size());
    for (std::size_t i = 0; i < _columns.size(); i++) {
      header[i] = _columns[i].name;
    }
    return header;
  }

  std::string MakeBorder(const std::vector<std::size_t>& widths,
                         const std::string_view left, const std::string_view middle,
                         const std::string_view right, const std::string_view horizontal) const {
    std::string result;
    result += left;

    for (std::size_t i = 0; i < widths.size(); ++i) {
      for (std::size_t j = 0; j < widths[i] + 2; ++j) {
        result += horizontal;
      }
      if (i + 1 < widths.size()) result += middle;
    }

    result += right;
    return result;
  }

  std::string MakeRow(const std::vector<std::size_t>& widths,
                      const std::vector<std::string>& row) const {
    std::string result = "│";

    for (std::size_t col = 0; col < widths.size(); ++col) {
      std::string value;

      if (col < row.size()){
        value = Truncate(row[col], widths[col]);
      }
      value = Align(value, widths[col], _columns[col].alignment);

      result += " ";
      result += value;
      result += " │";
    }
    return result;
  }
};
} // namespace M3::Utils
