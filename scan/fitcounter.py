#!/usr/bin/python3 -u

from scan.biotools import *


class FitCounter:
    def __init__(self, name):
        self.name = name
        self.counters_reset()
        self.timer_reset()

    def timer_reset(self):
        self.t0 = None
        self.t1 = None
        self.total_time = 0

    def timer_start(self):
        self.t0 = time.time()
        self.t1 = None
        self.total_time = 0

    def timer_finish(self):
        self.t1 = time.time()
        self.total_time = self.total_time + (self.t1 - self.t0)

    def counters_reset(self):
        self.count_files = 0
        self.count_skipped = 0
        self.count_errors = 0
        self.count_struct = 0
        self.count_res_tuple = 0
        self.count_hits = 0
        self.timer_start()

    def new_file(self):
        self.count_files = self.count_files + 1

    def new_skipped(self):
        self.count_skipped = self.count_skipped + 1

    def new_struct(self):
        self.count_struct = self.count_struct + 1

    def new_res_tuple(self):
        self.count_res_tuple = self.count_res_tuple + 1

    def new_hit(self):
        self.count_hits = self.count_hits + 1

    def new_error(self):
        self.count_errors = self.count_errors + 1

    def __str__(self):
        out = []
        out.append("fitscan statistics after {:>4d} sec:".format(int(self.total_time)))
        out.append("  {:>7} files with structutes".format(self.count_files))
        out.append("  {:>7} files skipped".format(self.count_skipped))
        out.append("  {:>7} structures in all models/chains".format(self.count_struct))
        out.append("  {:>7} hits".format(self.count_hits))
        out.append("  {:>7} {}".format(self.count_res_tuple, self.name))
        out.append("  {:>7} errors".format(self.count_errors))
        return "\n".join(out)

    def print_stats(self):
        self.timer_finish()
        print(str(self))

    def __add__(self, counter):
        ret: FitCounter = FitCounter(self.name)
        ret.total_time = self.total_time + counter.total_time
        ret.count_files = self.count_files + counter.count_files
        ret.count_skipped = self.count_skipped + counter.count_skipped
        ret.count_struct = self.count_struct + counter.count_struct
        ret.count_res_tuple = self.count_res_tuple + counter.count_res_tuple
        ret.count_hits = self.count_hits + counter.count_hits
        ret.count_errors = self.count_errors + counter.count_errors
        return ret
