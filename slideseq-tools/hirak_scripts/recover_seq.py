import numpy as np
import edlib
import sys
from collections import defaultdict
import timeit
import logging
from joblib import Parallel, delayed
import multiprocessing as mp

def setup_custom_logger(name):
    formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    #handler = logging.FileHandler('log.txt', mode='w')
    #handler.setFormatter(formatter)
    screen_handler = logging.StreamHandler(stream=sys.stdout)
    screen_handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    #logger.addHandler(handler)
    logger.addHandler(screen_handler)
    return logger


logger = setup_custom_logger('recover')
#bcode_filename = '/home/hsarkar/code/slide-seq-data/P25255/barcode_all.file'
#fp = open(bcode_filename)
#lines = [line.strip() for line in fp.readlines()]
logger.info('Strart reading from stream')
lines = [line.strip() for line in sys.stdin]
logger.info('Lines read from stream {}'.format(len(lines)))

barcode_matrix = '8C18X6C1X8M1X'
bi1 = 0
be1 = 0 + 8
si = 0 + 8 
se = 0 + 8 + 18
bi2 = 0 + 8 + 18
be2 = 0 + 8 + 18 + 6
ui = 0 + 8 + 18 + 6 + 1
ue = 0 + 8 + 18 + 6 + 1 + 8

# spacer_seq = {}
# for line in lines:
#     if line[si:se] in spacer_seq:
#         spacer_seq[line[si:se]] += 1
#     else:
#         spacer_seq[line[si:se]] = 1

spacer_seq = defaultdict(int)
for line in lines:
    spacer_seq[line[si:se]] += 1
freq_spacer = sorted(spacer_seq.items(), key=lambda item: item[1], reverse=True)[0][0]
logger.info('most frequent spacer {}'.format(freq_spacer))

logger.info('Start edlib aligner')
align_results = [
    edlib.align(
        freq_spacer, line, mode = "HW", task = 'locations',k=3
    ) for line in lines
]
logger.info('Without multiprocessing..')
# lines_seq = [(line, freq_spacer) for line in lines]
# # align_results = Parallel(n_jobs=4)(
# #     delayed(edlib.align)(
# #         line[1], 
# #         line[0], 
# #         mode = "HW", 
# #         task = 'locations',
# #         k=3
# #     ) for line in lines_seq
# # )
# def edlib_align(str_pair):
#     return(
#         edlib.align(
#             str_pair[1],
#             str_pair[0],
#             mode = 'HW',
#             task = 'locations',
#             k = 3
#         )
#     )

# pool = mp.Pool(processes=10)
# align_results = pool.map(edlib_align, lines_seq)
# logger.info('With multiprocessing..')

logger.info('Filtering list')
num_recovered = len(
    list(
    filter(lambda a: (
        a[1]['editDistance'] != -1 and 
        a[1]['locations'][0][0] >= 5 and 
        a[1]['locations'][0][0] <= 10), 
            enumerate(align_results)))
)

logger.info('Recovered {:.0%} lines'.format(num_recovered/len(lines)))