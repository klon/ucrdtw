
int ucrdtw(double* data, long long data_size, double* query, long query_size, double warp_width, int verbose, long long* location, double* distance);

int ucrdtwf(FILE* data_stream, FILE* query_stream, long query_length, double warp_width, int verbose, long long* location, double* distance);
