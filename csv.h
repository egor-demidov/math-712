//
// Created by egor on 9/8/24.
//

#ifndef CSV_H
#define CSV_H

#include <string>
#include <filesystem>
#include <fstream>
#include <array>
#include <iostream>

template <unsigned long n_cols>
class CSVWriter {
public:
    explicit CSVWriter(std::filesystem::path init_file_path, std::array<std::string, n_cols> init_col_names)
        : file_path{std::move(init_file_path)}
    , col_names{init_col_names}
    , ofs{file_path} {

        if (!ofs.good()) {
            std::cerr << "Unable to write to " << file_path << std::endl;
            exit(EXIT_FAILURE);
        }

        for (unsigned long i = 0; i < n_cols; i ++) {
            ofs << col_names[i];
            if (i < n_cols -1)
                ofs << ", ";
            else
                ofs << "\n";
        }
    }

    ~CSVWriter() {
        ofs << std::endl;
    }

    void append_line(std::array<double, n_cols> const & data) {
        for (unsigned long i = 0; i < n_cols; i ++) {
            ofs << data[i];
            if (i < n_cols -1)
                ofs << ", ";
            else
                ofs << "\n";
        }
    }

private:
    const std::filesystem::path file_path;
    const std::array<std::string, n_cols> col_names;
    std::ofstream ofs;
};

#endif //CSV_H
