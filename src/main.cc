// Copyright (C) 2021 Thomas Rodgers.  Use, modification, and distribution are
// subject to the Boost Software License, Version 1.0.  (See accompanying file
// LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

#include <algorithm>
#include <array>
#include <cstdio>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

template<size_t fwidth>
  std::string
  read_fixed_field(std::istream& stm)
  {
    std::array<char, fwidth> buf;
    buf[0] = 0;
    if (stm.read(buf.data(), buf.max_size()))
      {
         std::string_view s{ buf.data(), buf.max_size() };
         s.remove_prefix(std::min(s.find_first_not_of(" "), s.size()));
         return std::string{ s };
      }
    return std::string{ };
  }

class DEMInfo
{
public:
  DEMInfo() = default;
  DEMInfo(const DEMInfo&) = default;

  friend std::istream&
  operator>>(std::istream& stm, DEMInfo& inf)
  {
    // read in file name and description
    inf.fname_ = read_fixed_field<40>(stm);
    inf.description_ = read_fixed_field<40>(stm);

    // skip to the dimensions in rows and columns
    stm.seekg(853);
    stm >> inf.rows_;
    stm >> inf.cols_;

    // skip the rest for now
    stm.seekg(1024);

    return stm;
  }

  friend std::ostream&
  operator<<(std::ostream& stm, const DEMInfo& inf)
  {
    return stm << "FileName   : " << inf.fname_
               << "\nDescription: " << (!inf.description_.empty() ? inf.description_ : "<none>")
               << "\nRows       : " << inf.rows_
               << "\nColumns    : " << inf.cols_;
  }

  const std::string&
  file_name() const noexcept
  { return fname_; }

  const std::string&
  description() const noexcept
  { return description_; }

  std::pair<int, int>
  dimensions() const noexcept
  { return std::make_pair(rows_, cols_); }

private:
  std::string fname_;
  std::string description_;
  int rows_, cols_;
};

class DEMRecord
{
    using vector_t = std::vector<int>;

public:
  DEMRecord()
    : col_{ -1 }
  { }

  DEMRecord(DEMRecord&&) = default;

  using iterator = vector_t::iterator;
  using const_iterator = vector_t::const_iterator;

  iterator
  begin() noexcept
  { return std::begin(elevations_); }

  const_iterator
  begin() const noexcept
  { return std::begin(elevations_); }

  iterator
  end() noexcept
  { return std::end(elevations_); }

  const_iterator
  end() const noexcept
  { return std::end(elevations_); }

  friend std::istream&
  operator>>(std::istream& stm, DEMRecord& rec)
  {
    // there are 146 samples in the first block, 170 in each subsequent block
    constexpr const auto samples_block_0 = 146;
    constexpr const auto samples_block_n = 170;

    constexpr const auto header_width = 146;
    constexpr const auto int_field_width = 6;
    constexpr const auto block_width = 1024;

    // read row (and discard) and column
    auto start = stm.tellg(); // record current position in the stream
    int row;
    stm >> row >> rec.col_;

    int m, n;
    stm >> m >> n; // n is always 1
    stm.seekg(start, std::ios_base::beg); // seek back to starting point for this record

    // compute how many block_width blocks in this record
    auto l = (m * int_field_width) / block_width;
    if ((l % block_width) != 0)
      ++l;

    // allocate enough space to hold all samples in this record
    rec.elevations_.reserve(m);

    // create a buffer to hold the entire record
    std::vector<char> buf;
    buf.resize(l * block_width);

    // and slurp it in
    if (stm.read(buf.data(), buf.size()))
    {
      std::string_view s{ buf.data(), buf.max_size() };

      s.remove_prefix(header_width); // discard record header
      for (auto block = 0; block < l; ++block)
      {

        auto points = block ? samples_block_n : samples_block_0;
        for (auto i = 0; i < points; ++i, --m)
        {
          if (!m)
            break;

          std::string p{ s.begin(), s.begin()+int_field_width };
          s.remove_prefix(int_field_width);

          int z;
          if (!std::sscanf(p.c_str(), "%d", &z))
            break;
          rec.elevations_.emplace_back(z);
        }
        // the first block has 2 padding spaces before the subsequent block starts
        // subsequent blocks have 4 spaces.
        s.remove_prefix(block ? 4 : 2);
      }
    }

    return stm;
  }

  friend std::ostream&
  operator<<(std::ostream& stm, const DEMRecord& rec)
  {
    stm << "\nColumn: " << rec.col_
        << "\nPoints: " << rec.elevations_.size();

    int i = 0;
    const char* sep = "";
    for (auto z : rec.elevations_)
      {
        if ((i++ % 10) == 0)
          stm << '\n' << z;
        else
          stm << sep << z;
        sep = ", ";
      }
    return stm;
  }

private:
  int col_;
  std::vector<int> elevations_;
};

class DEMFile
{
  using vector_t = std::vector<DEMRecord>;

public:
  DEMFile() = default;
  DEMFile(const DEMFile&) = default;

  const DEMInfo&
  info() const noexcept
  { return info_; }

  using iterator = vector_t::iterator;
  using const_iterator = vector_t::const_iterator;

  iterator
  begin() noexcept
  { return std::begin(records_); }

  const_iterator
  begin() const noexcept
  { return std::begin(records_); }

  iterator
  end() noexcept
  { return std::end(records_); }

  const_iterator
  end() const noexcept
  { return std::end(records_); }

  friend std::istream&
  operator>>(std::istream& stm, DEMFile& f)
  {
    stm >> f.info_;

    const auto [_, cols] = f.info_.dimensions();
    f.records_.reserve(cols);

    for (auto i = 0; i < cols; ++i)
      {
        DEMRecord rec;
        stm >> rec;
        f.records_.emplace_back(std::move(rec));
      }

    return stm;
  }

  friend std::ostream&
  operator<<(std::ostream& stm, const DEMFile& f)
  {
    return stm << f.info_
               << "\nRecords    : " << f.records_.size();
  }

private:
  DEMInfo info_;
  vector_t records_;
};

int main(int argc, char const **argv)
{
  if (argc != 2)
    {
      std::cerr << "No input file specified." << std::endl;
      return -1;
    }

  fs::path fpath{ argv[1] };
  if (!fs::exists(fpath))
    {
      std::cerr << "No input file named " << fpath << std::endl;
      return -1;
    }

  std::cout << "Reading " << fpath << "...";
  std::cout.flush();
  std::ifstream stm{ fpath };
  DEMFile f;
  stm >> f;
  std::cout << "Done.\n";

  std::cout << f << std::endl;

  return 0;
}
