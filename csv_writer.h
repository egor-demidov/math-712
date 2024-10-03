//
// Created by egor on 10/3/24.
//

#ifndef CSV_WRITER_H
#define CSV_WRITER_H

#include <stdio.h>
#include <stdlib.h>

struct csv_writer {
    FILE * file;
    long n_cols;
};

typedef struct csv_writer csv_writer_t;

csv_writer_t open_csv_writer(const char * file_name, const char * col_names[], long n_cols);

void close_csv_writer(csv_writer_t writer);

void append_csv_line(csv_writer_t writer, double values[]);

#endif //CSV_WRITER_H
