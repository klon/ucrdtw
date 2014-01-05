/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected (c) 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "ucrdtw.h"

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) ((x-y)*(x-y))

#define INF 1e20 // Pseudo-infinite number for this code

/// Data structure for sorting the query
typedef struct index {
    double value;
    int index;
} index_t;

/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int index_comp(const void* a, const void* b) {
    index_t* x = (index_t*) a;
    index_t* y = (index_t*) b;
    return fabs(y->value) - fabs(x->value);   // high to low
}

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
typedef struct deque {
    int* dq;
    int size, capacity;
    int f, r;
} deque_t;

/// Initial the queue at the beginning step of envelope calculation
void deque_init(deque_t* d, int capacity) {
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int*) calloc(d->capacity, sizeof(int));
    d->f = 0;
    d->r = d->capacity - 1;
}

/// Destroy the queue
void deque_free(deque_t* d) {
    free(d->dq);
}

/// Insert to the queue at the back
void deque_push_back(deque_t* d, int v) {
    d->dq[d->r] = v;
    d->r--;
    if (d->r < 0)
        d->r = d->capacity - 1;
    d->size++;
}

/// Delete the current (front) element from queue
void deque_pop_front(deque_t* d) {
    d->f--;
    if (d->f < 0)
        d->f = d->capacity - 1;
    d->size--;
}

/// Delete the last element from queue
void deque_pop_back(deque_t* d) {
    d->r = (d->r + 1) % d->capacity;
    d->size--;
}

/// Get the value at the current position of the circular queue
int deque_front(deque_t* d) {
    int aux = d->f - 1;

    if (aux < 0)
        aux = d->capacity - 1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queue
int deque_back(deque_t* d) {
    int aux = (d->r + 1) % d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int deque_empty(deque_t* d) {
    return d->size == 0;
}

/// Finding the envelope of min and max value for LB_Keogh
/// Implementation idea is introduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void lower_upper_lemire(double *t, int len, int r, double *l, double *u) {
    struct deque du, dl;
    int i;

    deque_init(&du, 2 * r + 2);
    deque_init(&dl, 2 * r + 2);

    deque_push_back(&du, 0);
    deque_push_back(&dl, 0);

    for (i = 1; i < len; i++) {
        if (i > r) {
            u[i - r - 1] = t[deque_front(&du)];
            l[i - r - 1] = t[deque_front(&dl)];
        }
        if (t[i] > t[i - 1]) {
            deque_pop_back(&du);
            while (!deque_empty(&du) && t[i] > t[deque_back(&du)]) {
                deque_pop_back(&du);
            }
        } else {
            deque_pop_back(&dl);
            while (!deque_empty(&dl) && t[i] < t[deque_back(&dl)]) {
                deque_pop_back(&dl);
            }
        }
        deque_push_back(&du, i);
        deque_push_back(&dl, i);
        if (i == 2 * r + 1 + deque_front(&du)) {
            deque_pop_front(&du);
        } else if (i == 2 * r + 1 + deque_front(&dl)) {
            deque_pop_front(&dl);
        }
    }
    for (i = len; i < len + r + 1; i++) {
        u[i - r - 1] = t[deque_front(&du)];
        l[i - r - 1] = t[deque_front(&dl)];
        if (i - deque_front(&du) >= 2 * r + 1) {
            deque_pop_front(&du);
        }
        if (i - deque_front(&dl) >= 2 * r + 1) {
            deque_pop_front(&dl);
        }
    }
    deque_free(&du);
    deque_free(&dl);
}

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give significant benefits.
/// And using the first and last points can be computed in constant time.
/// The pruning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
double lb_kim_hierarchy(double *t, double *q, int j, int len, double mean, double std, double best_so_far) {
    /// 1 point at front and back
    double d, lb;
    double x0 = (t[j] - mean) / std;
    double y0 = (t[(len - 1 + j)] - mean) / std;
    lb = dist(x0,q[0]) + dist(y0, q[len - 1]);
    if (lb >= best_so_far)
        return lb;

    /// 2 points at front
    double x1 = (t[(j + 1)] - mean) / std;
    d = min(dist(x1,q[0]), dist(x0,q[1]));
    d = min(d, dist(x1,q[1]));
    lb += d;
    if (lb >= best_so_far)
        return lb;

    /// 2 points at back
    double y1 = (t[(len - 2 + j)] - mean) / std;
    d = min(dist(y1,q[len-1]), dist(y0, q[len-2]));
    d = min(d, dist(y1,q[len-2]));
    lb += d;
    if (lb >= best_so_far)
        return lb;

    /// 3 points at front
    double x2 = (t[(j + 2)] - mean) / std;
    d = min(dist(x0,q[2]), dist(x1, q[2]));
    d = min(d, dist(x2,q[2]));
    d = min(d, dist(x2,q[1]));
    d = min(d, dist(x2,q[0]));
    lb += d;
    if (lb >= best_so_far)
        return lb;

    /// 3 points at back
    double y2 = (t[(len - 3 + j)] - mean) / std;
    d = min(dist(y0,q[len-3]), dist(y1, q[len-3]));
    d = min(d, dist(y2,q[len-3]));
    d = min(d, dist(y2,q[len-2]));
    d = min(d, dist(y2,q[len-1]));
    lb += d;

    return lb;
}

/// LB_Keogh 1: Create Envelope for the query
/// Note that because the query is known, envelope can be created once at the beginning.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
double lb_keogh_cumulative(int* order, double *t, double *uo, double *lo, double *cb, int j, int len, double mean, double std, double best_so_far) {
    double lb = 0;
    double x, d;
    int i;

    for (i = 0; i < len && lb < best_so_far; i++) {
        x = (t[(order[i] + j)] - mean) / std;
        d = 0;
        if (x > uo[i]) {
            d = dist(x, uo[i]);
        } else if (x < lo[i]) {
            d = dist(x, lo[i]);
        }
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
/// l,u: lower and upper envelope of the current data
double lb_keogh_data_cumulative(int* order, double *tz, double *qo, double *cb, double *l, double *u, int len, double mean, double std, double best_so_far) {
    double lb = 0;
    double uu, ll, d;
    int i;

    for (i = 0; i < len && lb < best_so_far; i++) {
        uu = (u[order[i]] - mean) / std;
        ll = (l[order[i]] - mean) / std;
        d = 0;
        if (qo[i] > uu) {
            d = dist(qo[i], uu);
        } else {
            if (qo[i] < ll) {
                d = dist(qo[i], ll);
            }
        }
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double* A, double* B, double *cb, int m, int r, double best_so_far) {
    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i, j, k;
    double x, y, z, min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two arrays of size O(r).
    cost = (double*) calloc(2 * r + 1, sizeof(double));
    cost_prev = (double*) calloc(2 * r + 1, sizeof(double));
    for (k = 0; k < 2 * r + 1; k++) {
        cost[k] = INF;
        cost_prev[k] = INF;
    }

    for (i = 0; i < m; i++) {
        k = max(0, r - i);
        min_cost = INF;

        for (j = max(0, i - r); j <= min(m - 1, i + r); j++, k++) {
            /// Initialize all row and column
            if ((i == 0) && (j == 0)) {
                cost[k] = dist(A[0], B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j - 1 < 0) || (k - 1 < 0)) {
                y = INF;
            } else {
                y = cost[k - 1];
            }
            if ((i - 1 < 0) || (k + 1 > 2 * r)) {
                x = INF;
            } else {
                x = cost_prev[k + 1];
            }
            if ((i - 1 < 0) || (j - 1 < 0)) {
                z = INF;
            } else {
                z = cost_prev[k];
            }

            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + dist(A[i], B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost) {
                min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distance with lower bound together are larger than best_so_far
        if (i + r < m - 1 && min_cost + cb[i + r + 1] >= best_so_far) {
            free(cost);
            free(cost_prev);
            return min_cost + cb[i + r + 1];
        }

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}


/// Calculate the nearest neighbor of a times series in a larger time series expressed as location and distance,
/// using the UCR suite optimizations.
int ucrdtw(double* data, long long data_size, double* query, long query_size, double warp_width, int verbose, long long* location, double* distance) {
    long m = query_size;
    int r = warp_width <= 1 ? floor(warp_width * m) : floor(warp_width);

    double best_so_far; /// best-so-far
    double *q, *t; /// data array
    int *order; ///new order of the query
    double *u, *l, *qo, *uo, *lo, *tz, *cb, *cb1, *cb2, *u_d, *l_d;

    double d = 0.0;
    long long i, j;
    double ex, ex2, mean, std;

    long long loc = 0;
    double t1, t2;
    int kim = 0, keogh = 0, keogh2 = 0;
    double dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
    double *buffer, *u_buff, *l_buff;
    index_t *q_tmp;

    /// For every EPOCH points, all cumulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;

    if (verbose) {
        t1 = clock();
    }

    /// calloc everything here
    q = (double*) calloc(m, sizeof(double));
    if (q == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }
    memcpy((void*)q, (void*)query, m * sizeof(double));

    qo = (double*) calloc(m, sizeof(double));
    if (qo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    uo = (double*) calloc(m, sizeof(double));
    if (uo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    lo = (double*) calloc(m, sizeof(double));
    if (lo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    order = (int *) calloc(m, sizeof(int));
    if (order == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    q_tmp = (index_t *) calloc(m, sizeof(index_t));
    if (q_tmp == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    u = (double*) calloc(m, sizeof(double));
    if (u == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    l = (double*) calloc(m, sizeof(double));
    if (l == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    cb = (double*) calloc(m, sizeof(double));
    if (cb == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
    }

    cb1 = (double*) calloc(m, sizeof(double));
    if (cb1 == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    cb2 = (double*) calloc(m, sizeof(double));
    if (cb2 == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    u_d = (double*) calloc(m, sizeof(double));
    if (u == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    l_d = (double*) calloc(m, sizeof(double));
    if (l == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    t = (double*) calloc(m, sizeof(double) * 2);
    if (t == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    tz = (double*) calloc(m, sizeof(double));
    if (tz == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    buffer = (double*) calloc(EPOCH, sizeof(double));
    if (buffer == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    u_buff = (double*) calloc(EPOCH, sizeof(double));
    if (u_buff == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    l_buff = (double*) calloc(EPOCH, sizeof(double));
    if (l_buff == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    /// Read query
    best_so_far = INF;
    ex = ex2 = 0;
    for (i = 0; i < m; i++) {
        d = q[i];
        ex += d;
        ex2 += d * d;
    }
    /// Do z-normalize the query, keep in same array, q
    mean = ex / m;
    std = ex2 / m;
    std = sqrt(std - mean * mean);
    for (i = 0; i < m; i++)
        q[i] = (q[i] - mean) / std;

    /// Create envelope of the query: lower envelope, l, and upper envelope, u
    lower_upper_lemire(q, m, r, l, u);

    /// Sort the query one time by abs(z-norm(q[i]))
    for (i = 0; i < m; i++) {
        q_tmp[i].value = q[i];
        q_tmp[i].index = i;
    }
    qsort(q_tmp, m, sizeof(index_t), index_comp);

    /// also create another arrays for keeping sorted envelope
    for (i = 0; i < m; i++) {
        int o = q_tmp[i].index;
        order[i] = o;
        qo[i] = q[o];
        uo[i] = u[o];
        lo[i] = l[o];
    }

    /// Initial the cummulative lower bound
    for (i = 0; i < m; i++) {
        cb[i] = 0;
        cb1[i] = 0;
        cb2[i] = 0;
    }

    i = 0;          /// current index of the data in current chunk of size EPOCH
    j = 0;          /// the starting index of the data in the circular array, t
    ex = ex2 = 0;
    int done = 0;
    int it = 0, ep = 0, k = 0;
    long long I; /// the starting index of the data in current chunk of size EPOCH
    long long data_index = 0;
    while (!done) {
        /// Read first m-1 points
        ep = 0;
        if (it == 0) {
            for (k = 0; k < m - 1 && data_index < data_size; k++) {
                buffer[k] = data[data_index++];
            }
        } else {
            for (k = 0; k < m - 1; k++)
                buffer[k] = buffer[EPOCH - m + 1 + k];
        }

        /// Read buffer of size EPOCH or when all data has been read.
        ep = m - 1;
        while (ep < EPOCH && data_index < data_size) {
            buffer[ep] = data[data_index++];
            ep++;
        }

        /// Data are read in chunk of size EPOCH.
        /// When there is nothing to read, the loop is end.
        if (ep <= m - 1) {
            done = 1;
        } else {
            lower_upper_lemire(buffer, ep, r, l_buff, u_buff);
            /// Do main task here..
            ex = 0;
            ex2 = 0;
            for (i = 0; i < ep; i++) {
                /// A bunch of data has been read and pick one of them at a time to use
                d = buffer[i];

                /// Calcualte sum and sum square
                ex += d;
                ex2 += d * d;

                /// t is a circular array for keeping current data
                t[i % m] = d;

                /// Double the size for avoiding using modulo "%" operator
                t[(i % m) + m] = d;

                /// Start the task when there are more than m-1 points in the current chunk
                if (i >= m - 1) {
                    mean = ex / m;
                    std = ex2 / m;
                    std = sqrt(std - mean * mean);

                    /// compute the start location of the data in the current circular array, t
                    j = (i + 1) % m;
                    /// the start location of the data in the current chunk
                    I = i - (m - 1);

                    /// Use a constant lower bound to prune the obvious subsequence
                    lb_kim = lb_kim_hierarchy(t, q, j, m, mean, std, best_so_far);

                    if (lb_kim < best_so_far) {
                        /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                        /// uo, lo are envelope of the query.
                        lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, best_so_far);
                        if (lb_k < best_so_far) {
                            /// Take another linear time to compute z_normalization of t.
                            /// Note that for better optimization, this can merge to the previous function.
                            for (k = 0; k < m; k++) {
                                tz[k] = (t[(k + j)] - mean) / std;
                            }

                            /// Use another lb_keogh to prune
                            /// qo is the sorted query. tz is unsorted z_normalized data.
                            /// l_buff, u_buff are big envelope for all data in this chunk
                            lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff + I, u_buff + I, m, mean, std, best_so_far);
                            if (lb_k2 < best_so_far) {
                                /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                /// Note that cb and cb2 will be cumulative summed here.
                                if (lb_k > lb_k2) {
                                    cb[m - 1] = cb1[m - 1];
                                    for (k = m - 2; k >= 0; k--)
                                        cb[k] = cb[k + 1] + cb1[k];
                                } else {
                                    cb[m - 1] = cb2[m - 1];
                                    for (k = m - 2; k >= 0; k--)
                                        cb[k] = cb[k + 1] + cb2[k];
                                }

                                /// Compute DTW and early abandoning if possible
                                dist = dtw(tz, q, cb, m, r, best_so_far);

                                if (dist < best_so_far) {   /// Update best_so_far
                                                    /// loc is the real starting location of the nearest neighbor in the file
                                    best_so_far = dist;
                                    loc = (it) * (EPOCH - m + 1) + i - m + 1;
                                }
                            } else
                                keogh2++;
                        } else
                            keogh++;
                    } else
                        kim++;

                    /// Reduce absolute points from sum and sum square
                    ex -= t[j];
                    ex2 -= t[j] * t[j];
                }
            }

            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
            if (ep < EPOCH)
                done = 1;
            else
                it++;
        }
    }

    i = (it) * (EPOCH - m + 1) + ep;

    free(q);
    free(qo);
    free(uo);
    free(lo);
    free(order);
    free(q_tmp);
    free(u);
    free(l);
    free(cb);
    free(cb1);
    free(cb2);
    free(u_d);
    free(l_d);
    free(t);
    free(tz);
    free(buffer);
    free(u_buff);
    free(l_buff);

    if (verbose) {
        t2 = clock();
        printf("Location : %lld\n", loc);
        printf("Distance : %.6f\n", sqrt(best_so_far));
        printf("Data Scanned : %lld\n", i);
        printf("Total Execution Time : %.4f secs\n", (t2 - t1) / CLOCKS_PER_SEC);
        printf("\n");
        printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i) * 100);
        printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i) * 100);
        printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i) * 100);
        printf("DTW Calculation     : %6.2f%%\n", 100 - (((double) kim + keogh + keogh2) / i * 100));
    }
    *location = loc;
    *distance = sqrt(best_so_far);
    return 0;
}

/// Calculate the nearest neighbor of a times series in a larger time series expressed as location and distance,
/// using the UCR suite optimizations. This function supports streaming the data and the query to search.
int ucrdtwf(FILE* data_file, FILE* query_file, long query_length, double warp_width, int verbose, long long* location, double* distance) {
    FILE* fp = data_file;
    FILE* qp = query_file;
    long m = query_length;
    int r = warp_width <= 1 ? floor(warp_width * m) : floor(warp_width);

    double best_so_far;          /// best-so-far
    double *t, *q;       /// data array and query array
    int *order;          ///new order of the query
    double *u, *l, *qo, *uo, *lo, *tz, *cb, *cb1, *cb2, *u_d, *l_d;

    double d;
    long long i, j;
    double ex, ex2, mean, std;

    long long loc = 0;
    double t1, t2;
    int kim = 0, keogh = 0, keogh2 = 0;
    double dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
    double *buffer, *u_buff, *l_buff;
    index_t *q_tmp;

    /// For every EPOCH points, all cumulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;

    if (verbose) {
        t1 = clock();
    }

    /// calloc everything here
    q = (double*) calloc(m, sizeof(double));
    if (q == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    qo = (double*) calloc(m, sizeof(double));
    if (qo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    uo = (double*) calloc(m, sizeof(double));
    if (uo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    lo = (double*) calloc(m, sizeof(double));
    if (lo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    order = (int*) calloc(m, sizeof(int));
    if (order == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    q_tmp = (index_t*) calloc(m, sizeof(index_t));
    if (q_tmp == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    u = (double*) calloc(m, sizeof(double));
    if (u == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    l = (double*) calloc(m, sizeof(double));
    if (l == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    cb = (double*) calloc(m, sizeof(double));
    if (cb == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
    }

    cb1 = (double*) calloc(m, sizeof(double));
    if (cb1 == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    cb2 = (double*) calloc(m, sizeof(double));
    if (cb2 == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    u_d = (double*) calloc(m, sizeof(double));
    if (u == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    l_d = (double*) calloc(m, sizeof(double));
    if (l == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    t = (double*) calloc(m, sizeof(double) * 2);
    if (t == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    tz = (double*) calloc(m, sizeof(double));
    if (tz == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    buffer = (double*) calloc(EPOCH, sizeof(double));
    if (buffer == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    u_buff = (double*) calloc(EPOCH, sizeof(double));
    if (u_buff == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    l_buff = (double*) calloc(EPOCH, sizeof(double));
    if (l_buff == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    /// Read query file
    best_so_far = INF;
    i = 0;
    j = 0;
    ex = ex2 = 0;

    while (fscanf(qp, "%lf", &d) != EOF && i < m) {
        ex += d;
        ex2 += d * d;
        q[i] = d;
        i++;
    }

    /// Do z-normalize the query, keep in same array, q
    mean = ex / m;
    std = ex2 / m;
    std = sqrt(std - mean * mean);
    for (i = 0; i < m; i++)
        q[i] = (q[i] - mean) / std;

    /// Create envelope of the query: lower envelope, l, and upper envelope, u
    lower_upper_lemire(q, m, r, l, u);

    /// Sort the query one time by abs(z-norm(q[i]))
    for (i = 0; i < m; i++) {
        q_tmp[i].value = q[i];
        q_tmp[i].index = i;
    }
    qsort(q_tmp, m, sizeof(index_t), index_comp);

    /// also create another arrays for keeping sorted envelope
    for (i = 0; i < m; i++) {
        int o = q_tmp[i].index;
        order[i] = o;
        qo[i] = q[o];
        uo[i] = u[o];
        lo[i] = l[o];
    }

    /// Initial the cummulative lower bound
    for (i = 0; i < m; i++) {
        cb[i] = 0;
        cb1[i] = 0;
        cb2[i] = 0;
    }

    i = 0;          /// current index of the data in current chunk of size EPOCH
    j = 0;          /// the starting index of the data in the circular array, t
    ex = ex2 = 0;
    int done = 0;
    int it = 0, ep = 0, k = 0;
    long long I; /// the starting index of the data in current chunk of size EPOCH

    while (!done) {
        /// Read first m-1 points
        ep = 0;
        if (it == 0) {
            for (k = 0; k < m - 1; k++)
                if (fscanf(fp, "%lf", &d) != EOF)
                    buffer[k] = d;
        } else {
            for (k = 0; k < m - 1; k++)
                buffer[k] = buffer[EPOCH - m + 1 + k];
        }

        /// Read buffer of size EPOCH or when all data has been read.
        ep = m - 1;
        while (ep < EPOCH) {
            if (fscanf(fp, "%lf", &d) == EOF)
                break;
            buffer[ep] = d;
            ep++;
        }

        /// Data are read in chunk of size EPOCH.
        /// When there is nothing to read, the loop is end.
        if (ep <= m - 1) {
            done = 1;
        } else {
            lower_upper_lemire(buffer, ep, r, l_buff, u_buff);
            /// Do main task here..
            ex = 0;
            ex2 = 0;
            for (i = 0; i < ep; i++) {
                /// A bunch of data has been read and pick one of them at a time to use
                d = buffer[i];

                /// Calcualte sum and sum square
                ex += d;
                ex2 += d * d;

                /// t is a circular array for keeping current data
                t[i % m] = d;

                /// Double the size for avoiding using modulo "%" operator
                t[(i % m) + m] = d;

                /// Start the task when there are more than m-1 points in the current chunk
                if (i >= m - 1) {
                    mean = ex / m;
                    std = ex2 / m;
                    std = sqrt(std - mean * mean);

                    /// compute the start location of the data in the current circular array, t
                    j = (i + 1) % m;
                    /// the start location of the data in the current chunk
                    I = i - (m - 1);

                    /// Use a constant lower bound to prune the obvious subsequence
                    lb_kim = lb_kim_hierarchy(t, q, j, m, mean, std, best_so_far);

                    if (lb_kim < best_so_far) {
                        /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                        /// uo, lo are envelope of the query.
                        lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, best_so_far);
                        if (lb_k < best_so_far) {
                            /// Take another linear time to compute z_normalization of t.
                            /// Note that for better optimization, this can merge to the previous function.
                            for (k = 0; k < m; k++) {
                                tz[k] = (t[(k + j)] - mean) / std;
                            }

                            /// Use another lb_keogh to prune
                            /// qo is the sorted query. tz is unsorted z_normalized data.
                            /// l_buff, u_buff are big envelope for all data in this chunk
                            lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff + I, u_buff + I, m, mean, std, best_so_far);
                            if (lb_k2 < best_so_far) {
                                /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                /// Note that cb and cb2 will be cumulative summed here.
                                if (lb_k > lb_k2) {
                                    cb[m - 1] = cb1[m - 1];
                                    for (k = m - 2; k >= 0; k--)
                                        cb[k] = cb[k + 1] + cb1[k];
                                } else {
                                    cb[m - 1] = cb2[m - 1];
                                    for (k = m - 2; k >= 0; k--)
                                        cb[k] = cb[k + 1] + cb2[k];
                                }

                                /// Compute DTW and early abandoning if possible
                                dist = dtw(tz, q, cb, m, r, best_so_far);

                                if (dist < best_so_far) {   /// Update best_so_far
                                                    /// loc is the real starting location of the nearest neighbor in the file
                                    best_so_far = dist;
                                    loc = (it) * (EPOCH - m + 1) + i - m + 1;
                                }
                            } else
                                keogh2++;
                        } else
                            keogh++;
                    } else
                        kim++;

                    /// Reduce absolute points from sum and sum square
                    ex -= t[j];
                    ex2 -= t[j] * t[j];
                }
            }

            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
            if (ep < EPOCH)
                done = 1;
            else
                it++;
        }
    }

    i = (it) * (EPOCH - m + 1) + ep;

    free(q);
    free(qo);
    free(uo);
    free(lo);
    free(order);
    free(q_tmp);
    free(u);
    free(l);
    free(cb);
    free(cb1);
    free(cb2);
    free(u_d);
    free(l_d);
    free(t);
    free(tz);
    free(buffer);
    free(u_buff);
    free(l_buff);

    if (verbose) {
        t2 = clock();
        printf("Data Scanned : %lld\n", i);
        printf("Total Execution Time : %.4f secs\n", (t2 - t1) / CLOCKS_PER_SEC);
        printf("\n");
        printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i) * 100);
        printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i) * 100);
        printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i) * 100);
        printf("DTW Calculation     : %6.2f%%\n", 100 - (((double) kim + keogh + keogh2) / i * 100));
    }
    *location = loc;
    *distance = sqrt(best_so_far);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        printf("Error: Invalid arguments\n");
        printf("Command usage: dtwfind <data-file>  <query-file> <query-length> <warp-width>\n\n");
        printf("For example: dtwfind data.txt query.txt 128 0.05\n");
        return -2;
    }

    FILE* data_file = fopen(argv[1], "r");
    if (data_file == NULL) {
        printf("ERROR: File not found: %s\n", argv[1]);
    }

    FILE* query_file = fopen(argv[2], "r");
    if (query_file == NULL) {
        printf("ERROR: File not found: %s\n", argv[2]);
    }

    long long location = -1;
    double distance = -1;
    int ret = ucrdtwf(data_file, query_file, atoll(argv[3]), atof(argv[4]), 1, &location, &distance);
    printf("Location: %lld\nDistance: %.6f\n", location, distance);
    return ret;
}
