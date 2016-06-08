import re
import pysam

def natural_key(s):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', s)]


b37_Chrom_Names = []


def get_reference_names_from_bam(filename):
    res = []
    with pysam.AlignmentFile(filename, "rb") as f:
        res = [x['SN'] for x in f.header['SQ']]

    return res

class GenomePosition(object):

    """docstring for GenomePosition"""

    def __init__(self, ref_name, pos):
        self.ref_name = ref_name
        self.pos = pos

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.ref_name, self.pos) == (other.ref_name, other.pos)
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return (natural_key(self.ref_name), self.pos) < \
                (natural_key(other.ref_name), other.pos)
        return NotImplemented

    def __hash__(self):
        return hash((self.ref_name, self.pos))

    def __str__(self):
        return '{}:{}'.format(self.ref_name, self.pos)

    def reference_name(self):
        return self.ref_name

    def position(self):
        return self.pos

    def genome_position_with_ci(self, slop):
        return GenomePositionWithCi(self, Interval(-slop, slop))


class Interval(object):

    """docstring for Interval"""

    def __init__(self, a, b):
        assert a <= b
        self.a = a
        self.b = b

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.a, self.b) == (other.a, other.b)
        return NotImplemented

    def __len__(self):
        return self.b - self.a

    def __str__(self):
        return '[{}, {}]'.format(self.a, self.b)

    def before(self, pos):
        return pos < self.a

    def covers(self, pos):
        return pos >= self.a and pos <= self.b

    def overlaps(self, other):
        return (self.a >= other.a and self.a < other.b) or \
            (other.a >= self.a and other.a < self.b)

    def reciprocal_overlaps(self, other, threshold):
        if not self.overlaps(other):
            return False

        o = min(self.b, other.b) - max(self.a, other.a)
        n = len(self)
        m = len(other)

        return float(o)/max(n, m) >= threshold


class GenomePositionWithCi(object):

    """docstring for GenomePositionWithCi"""

    def __init__(self, genome_pos, confidence_interval):
        self.genome_pos = genome_pos
        self.confidence_interval = confidence_interval

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.genome_pos, self.confidence_interval) == (other.genome_pos, other.confidence_interval)
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self.genome_pos < other.genome_pos
        return NotImplemented

    def __str__(self):
        return '{}, {}'.format(self.genome_pos, self.confidence_interval)

    def __hash__(self):
        return hash(self.genome_pos)

    def reference_name(self):
        return self.genome_pos.ref_name

    def position(self):
        return self.genome_pos.pos

    def matches(self, g_pos):
        return self.to_genome_region().covers(g_pos)

    def to_genome_region(self):
        return GenomeRegion(
            GenomePosition(
                self.reference_name(), self.position() + self.confidence_interval.a),
            GenomePosition(
                self.reference_name(), self.position() + self.confidence_interval.b)
        )


class GenomeRegion(object):

    """docstring for GenomeRegion"""

    def __init__(self, start, end):
        assert start.ref_name == end.ref_name
        self.start = start
        self.end = end

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.start, self.end) == (other.start, other.end)
        return NotImplemented

    def reference_name(self):
        return self.start.ref_name

    def start_pos(self):
        return self.start.pos

    def end_pos(self):
        return self.end.pos

    def before(self, g_pos):
        if (self.start.reference_name() == g_pos.reference_name()):
            return Interval(self.start.position(), self.end.position()).before(g_pos.position())

    def covers(self, g_pos):
        return self.start.reference_name() == g_pos.reference_name() \
        and Interval(self.start.position(), self.end.position()).covers(g_pos.position())

    def overlaps(self, other):
        if self.start.ref_name != other.start.ref_name or self.end.ref_name != other.end.ref_name:
            return False
        return Interval(self.start.pos, self.end.pos).overlaps(Interval(other.start.pos, other.end.pos))

    def reciprocal_overlaps(self, other, threshold):
        if self.start.ref_name != other.start.ref_name or self.end.ref_name != other.end.ref_name:
            return False
        return Interval(self.start.pos, self.end.pos).reciprocal_overlaps(Interval(other.start.pos, other.end.pos), threshold)
