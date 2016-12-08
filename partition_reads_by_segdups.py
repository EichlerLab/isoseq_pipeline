from __future__ import print_function
from __future__ import division

import tabix
from intervaltree import Interval, IntervalTree

import attr
from attr.validators import instance_of
import bisect
import argparse

@attr.s
class BedInterval(object):
    "Class for a genomic region with chr, start, and end"

    chr = attr.ib(validator=instance_of(str))
    start = attr.ib(validator=instance_of(int))
    end = attr.ib(validator=instance_of(int))
    name = attr.ib(default="", validator=instance_of(str), cmp=False)
    reads = attr.ib(default=0, validator=instance_of(int), cmp=False)
    
    def to_region(self):
        "Returns list in format for pytabix queries"
        return (self.chr, self.start, self.end)

    def to_interval(self):
        "Returns intervaltree Interval with start and end"
        return Interval(start, end)

@attr.s
class Partition(object):
    """Class for creating read partitions based on segdups"""

    conflicts = attr.ib(default=attr.Factory(dict), cmp=False)
    intervals = attr.ib(default=attr.Factory(list), cmp=False)
    reads = attr.ib(default=0, validator=instance_of(int))

    def __len__(self):
        return self.reads

    def add_segdups(self, segdups):
        for segdup in segdups:
            if segdup.chr not in self.conflicts:
                self.conflicts[segdup.chr] = IntervalTree()
            self.conflicts[segdup.chr].add(segdup.to_interval())

    def merge_conflicts(self, conflicts):
        for chr, tree in conflicts.items():
            if chr not in self.conflicts:
                self.conflicts[chr] = tree
            else:
                self.conflicts[chr] |= tree # IntervalTree Union

    def add_interval(self, interval, segdups):
        """Takes a BedInterval.
        Adds new intersecting segdups to conflicts,
        adds reads to read counter."""
        self.add_segdups(segdups)
        self.reads += interval.reads
        self.intervals.append(interval)

    def partition_conflicts(self, partition):
        """Returns True if partitions have overlapping segdup conflicts.
        Returns False otherwise."""

        for chr, tree in partition.conflicts.items():
            if chr not in self.conflicts:
                continue
            for iv in partition.conflicts[chr]:
                if self.conflicts[chr].overlaps(iv):
                    return True
        return False

    def merge(self, partition):
        """Add conflicts, intervals, and reads of partition to self"""
        self.merge_conflicts(partition.conflicts)
        self.reads += partition.reads
        self.intervals.extend(partition.intervals)

@attr.s
class PartitionSet(object):
    """Maintains a set of partitions that don't have internal segdup conflicts.
    Adds regions to partition until it contains the minimum number of reads."""
    segdups = attr.ib(validator=instance_of(tabix.open))
    read_threshold = attr.ib(default=200, validator=instance_of(int))
    small_partitions = attr.ib(default=attr.Factory(list))
    full_partitions = attr.ib(default=attr.Factory(list))
    
    def add_interval(self, interval):
        """Selects a partition to add the interval to based on size and lack of segdup conflicts.
        If no suitable interval exists, create a new one."""

        segdups = self.segdup_matches(interval)
        index = self.find_partition_index_for_interval(segdups, interval)
        if index is not None:
            partition = self.small_partitions.pop(index)
        else:
            partition = Partition()

        partition.add_interval(interval, segdups)
        if len(partition) < self.read_threshold:
            bisect.insort_left(self.small_partitions, partition)
        else:
            bisect.insort_left(self.full_partitions, partition)

    def segdup_matches(self, interval):
        "Returns a list of segdup BedIntervals for a given interval"
        segdups = []
        for record in self.segdups.query(*interval.to_region()):
            chr, se = record[3].split(":")
            start, end = map(int, se.split("-"))
            segdups.append(BedInterval(chr, start, end, name=se))
        return segdups

    def segdup_conflicts(self, segdups, partition):
        "Returns whether an interval of segdups overlaps conflict regions in given partition"
        for segdup in segdups:
            if segdup.chr not in partition.conflicts:
                continue
            if partition.conflicts[segdup.chr].overlaps(segdup.to_interval()):
                return True
        return False

    def find_partition_index_for_interval(self, segdups, interval):
        """Return index of partition to add interval to.
        Returns None if no suitable partition exists."""
        for i, partition in enumerate(self.small_partitions):
            if not self.segdup_conflicts(segdups, partition):
                return i
        return None

    def find_partition_index_for_partition(self, partition):
        """Return index of partition query partition can be added to.
        Returns None if no suitable partition exists.
        """

        for i, part in enumerate(self.small_partitions):
            if not part.partition_conflicts(partition):
                return i

        for i, part in enumerate(self.full_partitions):
            if not part.partition_conflicts(partition):
                return i + len(self.small_partitions)
        return None

    def finish_small_partitions(self):
        """Method to be called after all intervals have been added.
        Combines small partitions with each other and with full partitions
        to reduce number of small partitions.
        """

        while len(self.small_partitions) > 0:
            partition = self.small_partitions.pop()
            index = self.find_partition_index_for_partition(partition)

            if index < len(self.small_partitions):
                part = self.small_partitions.pop(index)
            else:
                part = self.full_partitions.pop(index - len(self.small_partitions))

            if part is None:
                bisect.insort_left(self.full_partitions, partition)
                continue
            part.merge(partition)
            if partition.reads < self.read_threshold:
                bisect.insort_left(self.small_partitions, part)
            else:
                bisect.insort_left(self.full_partitions, part)

    def write(self, outfile):
        assert len(self.small_partitions) == 0
        print("chr", "start", "end", "group", "reads", sep="\t", file=outfile)
        for i, part in enumerate(self.full_partitions):
            for iv in part.intervals:
                print(iv.chr, iv.start, iv.end, i, part.reads, sep="\t", file=outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("segdups", help="Path to bgzipped and tabixed bed file of segdups \
                        with chr, start, end of query and chr:start-end of match in name field")
    parser.add_argument("read_locations", help="Bed file with chr, start, end, and nreads of initial groups")
    parser.add_argument("outfile", help="Path to tab-delimited output file in headered bed format")
    parser.add_argument("--count_threshold", default=200, type=int, help="Target minimum partition size")

    args = parser.parse_args()

    segdups = tabix.open(args.segdups)

    partitions = PartitionSet(segdups, args.count_threshold)

    with open(args.read_locations, "r") as regions:
        for region in regions:
            chr, start, end, nreads = region.split()
            interval = BedInterval(chr, int(start), int(end), reads=int(nreads))
            partitions.add_interval(interval)
    partitions.finish_small_partitions()
    partitions.write(open(args.outfile, "w"))
