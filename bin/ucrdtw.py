import _ucrdtw
import numpy
import sys
import argparse
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate best DTW location and distance')
    parser.add_argument('data', metavar='DATA_FILE', type=str, help='Path to data file')
    parser.add_argument('query', metavar='QUERY_FILE', type=str, help='Path to query file')
    parser.add_argument('query_size', metavar='QUERY_SIZE', type=int, default=0, help='Max size of query')
    parser.add_argument('warp_width', metavar='WARP_WIDTH', type=float, default=0.05, help='Width of allowed warp as fraction of query size')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Print verbose info')
    args = parser.parse_args()

    data = numpy.fromfile(args.data, sep=' ')
    query = numpy.fromfile(args.query, sep=' ', count=args.query_size if args.query_size > 0 else -1)
    print _ucrdtw.ucrdtw(data, query, args.warp_width, True)

