class ChromNameIdConverter(object):
    """docstring for ChromNameIdConverter"""

    Chrom_names = [
        '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
        '21', '22', 'X', 'Y', 'MT', 'GL000207.1', 'GL000226.1',
        'GL000229.1', 'GL000231.1', 'GL000210.1', 'GL000239.1',
        'GL000235.1', 'GL000201.1', 'GL000247.1', 'GL000245.1',
        'GL000197.1', 'GL000203.1', 'GL000246.1', 'GL000249.1',
        'GL000196.1', 'GL000248.1', 'GL000244.1', 'GL000238.1',
        'GL000202.1', 'GL000234.1', 'GL000232.1', 'GL000206.1',
        'GL000240.1', 'GL000236.1', 'GL000241.1', 'GL000243.1',
        'GL000242.1', 'GL000230.1', 'GL000237.1', 'GL000233.1',
        'GL000204.1', 'GL000198.1', 'GL000208.1', 'GL000191.1',
        'GL000227.1', 'GL000228.1', 'GL000214.1', 'GL000221.1',
        'GL000209.1', 'GL000218.1', 'GL000220.1', 'GL000213.1',
        'GL000211.1', 'GL000199.1', 'GL000217.1', 'GL000216.1',
        'GL000215.1', 'GL000205.1', 'GL000219.1', 'GL000224.1',
        'GL000223.1', 'GL000195.1', 'GL000212.1', 'GL000222.1',
        'GL000200.1', 'GL000193.1', 'GL000194.1', 'GL000225.1',
        'GL000192.1'
    ]

    @classmethod
    def name_to_id(cls, chrom_name):
        return cls.Chrom_names.index(chrom_name)

    @classmethod
    def is_valid(cls, chrom_name):
        return chrom_name in cls.Chrom_names


class GenomePosition(object):

    """docstring for GenomePosition"""

    def __init__(self, ref_id, ref_name, pos):
        self.ref_id = ref_id
        self.ref_name = ref_name
        self.pos = pos

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.ref_id, self.pos) == (other.ref_id, other.pos)
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return (self.ref_id, self.pos) < \
                (other.ref_id, other.pos)
        return NotImplemented

    def __hash__(self):
        return hash((self.ref_id, self.pos))

    def __str__(self):
        return '{}:{}'.format(self.ref_name, self.pos)

    def reference_id(self):
        return self.ref_id

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

    def after(self, pos):
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

    def reference_id(self):
        return self.genome_pos.ref_id

    def reference_name(self):
        return self.genome_pos.ref_name

    def position(self):
        return self.genome_pos.pos

    def after(self, g_pos):
        return self.to_genome_region().after(g_pos)

    def matches(self, g_pos):
        return self.to_genome_region().covers(g_pos)

    def to_genome_region(self):
        return GenomeRegion(
            GenomePosition(
                self.reference_id(),
                self.reference_name(),
                self.position() + self.confidence_interval.a
            ),
            GenomePosition(
                self.reference_id(),
                self.reference_name(),
                self.position() + self.confidence_interval.b
            )
        )


class GenomeRegion(object):

    """docstring for GenomeRegion"""

    def __init__(self, start, end):
        """
        Args:
            start: its class must be GenomePosition.
            end: its class must be GenomePosition.
        """
        assert isinstance(start, GenomePosition) and isinstance(end, GenomePosition)
        assert start.ref_id == end.ref_id
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

    def after(self, g_pos):
        """
        Args:
            g_pos: its class can be GenomePositionWithCi other than GenomePosition.
        """
        if self.start.ref_id == g_pos.reference_id():
            return Interval(self.start.pos, self.end.pos).after(g_pos.position())
        return self.start.ref_id > g_pos.reference_id()

    def covers(self, g_pos):
        """
        Args:
            g_pos: its class can be GenomePositionWithCi other than GenomePosition.
        """
        return self.start.ref_id == g_pos.reference_id() \
            and Interval(self.start.pos, self.end.pos).covers(g_pos.position())

    def overlaps(self, other):
        if self.start.ref_id != other.start.ref_id or self.end.ref_id != other.end.ref_id:
            return False
        return Interval(self.start.pos, self.end.pos).overlaps(Interval(other.start.pos, other.end.pos))

    def reciprocal_overlaps(self, other, threshold):
        if self.start.ref_id != other.start.ref_id or self.end.ref_id != other.end.ref_id:
            return False
        return Interval(self.start.pos, self.end.pos).reciprocal_overlaps(Interval(other.start.pos, other.end.pos), threshold)
