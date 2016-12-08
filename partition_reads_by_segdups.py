from __future__ import print_function
from __future__ import division

import tabix
from intervaltree import IntervalTree

import attr
from attr.validators import instance_of
import bisect
import argparse

@attr.s
class Interval(object):
    "Class for a genomic region with chr, start, and end"

    chr = attr.ib(validator=instance_of(str))
    start = attr.ib(validator=instance_of(int))
    end = attr.ib(validator=instance_of(int))
    reads = attr.ib(default=None, validator=instance_of(int), cmp=False)
    name = attr.ib(default=None, validator=instance_of(str), cmp=False)

    def toregion(self):
        "Returns list in format for pytabix queries"
        return (self.chr, self.start, self.end)

@attr.s
class Partition(object):
    """Class for creating read partitions based on segdups"""

    conflicts = attr.ib(default=attr.Factory(dict), cmp=False)
    intervals = attr.ib(default=attr.Factory(list), cmp=False)
    reads = attr.ib(default=0, validator=instance_of(int))

    def add_interval(self, interval, segdups):
        """Takes an Interval.
        Adds new intersecting segdups to conflicts,
        adds reads to read counter."""
        for segdup in segdups:
            if segdup.chr not in self.conflicts:
                self.conflicts[segdup.chr] = intervaltree.IntervalTree()
            self.conflicts[segdup.chr].add(segdup.start, segdup.end, segdup.name)
        self.reads += interval.reads
        bisect.insort_left(interval, self.intervals)

@attr.s
class PartitionSet(object):
    """Maintains a set of partitions that don't have internal segdup conflicts.
    Adds regions to partition until it contains the minimum number of reads."""
    small_partitions = attr.ib(default=attr.Factory(list))
    full_partitions = attr.ib(default=attr.Factory(list))
    segdups = attr.ib(default=attr.Factory(tabix.open))
    read_threshold = attr.ib(default=200, validator=instance_of(int))

    def add_interval(self, interval):
        """Selects a partition to add the interval to based on size and lack of segdup conflicts.
        If no suitable interval exists, create a new one."""

        segdups = self.segdup_matches(interval)
        index = self.find_partition_index(interval)
        if index is not None:
            partition = self.small_partitions.pop(index)
        else:
            partition = Partition()
        partition.add_interval()
        if partition.reads < self.read_threshold:
            # Insert in sorted order
            i = bisect.insort_left(partition, self.small_partitions)
        else:
            # Insert in order of filling
            self.full_partitions.append(partition)

    def segdup_matches(self, interval):
        "Returns a list of segdup Intervals for a given interval"
        segdups = []
        for record in self.segdups.query(*interval.toregion()):
            chr, se = record[3].split(":")
            start, end = map(int, se.split("-"))
            segdups.append(Interval(chr, start, end, name=se))
        return segdups

    def segdup_conflicts(self, segdups, partition):
        "Returns whether an interval of segdups overlaps conflict regions in given partition"
        for segdup in segdups:
            if segdup.chr not in partition.conflicts:
                continue
            if partition.conflicts.overlaps(segdup.start, segdup.end):
                return True
        return False

    def find_partition_index(self, interval):
        """Return index of partition to add interval to.
        Returns None if no suitable partition exists."""
        for i, partition in self.small_partitions:
            if not self.segdup_conflict(segdups, partition):
                return i
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("segdups")
    parser.add_argument("read_locations")
    parser.add_argument("--count_threshold", default=200, type=int)

    args = parser.parse_args()

