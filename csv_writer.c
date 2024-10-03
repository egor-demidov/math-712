//
// Created by egor on 10/3/24.
//

#include "csv_writer.h"

csv_writer_t open_csv_writer(const char * file_name, const char * col_names[], long n_cols) {
    csv_writer_t writer;
    writer.file = fopen(file_name, "w");

    if (writer.file == NULL) {
        fprintf(stderr, "Unable to create a CSV writer to %s\n", file_name);
        exit(EXIT_FAILURE);
    }

    writer.n_cols = n_cols;

    for (long n = 0; n < writer.n_cols - 1; n ++) {
        fprintf(writer.file, "%s, ", col_names[n]);
    }
    fprintf(writer.file, "%s\n", col_names[writer.n_cols - 1]);

    return writer;
}

void close_csv_writer(csv_writer_t writer) {
    fclose(writer.file);
}

void append_csv_line(csv_writer_t writer, double values[]) {
    for (long n = 0; n < writer.n_cols - 1; n ++) {
        fprintf(writer.file, "%lf, ", values[n]);
    }
    fprintf(writer.file, "%lf\n", values[writer.n_cols - 1]);
}
