from abc import ABCMeta, abstractmethod
import re
from generic import GenomePosition, GenomeRegion, GenomePositionWithCi, Interval


class VariantTypes(object):
    """docstring for VariantType"""
    DEL = 0
    INS = 1
    INV = 2
    DUP = 3


class Variant(object):
    """docstring for Variant"""

    __metaclass__ = ABCMeta

    def __init__(self, name, pos1, pos2, micro_hom_seq=None, micro_ins_seq=None):
        self.name = name
        self.pos1 = pos1
        self.pos2 = pos2
        self.micro_hom_seq = micro_hom_seq
        self.micro_ins_seq = micro_ins_seq

    def __eq__(self, other):
        return (self.variant_type(), self.pos1, self.pos2) == \
            (other.variant_type(), other.pos1, other.pos2)

    def __lt__(self, other):
        return (self.pos1, self.pos2) < (other.pos1, other.pos2)

    def __hash__(self):
        return hash((self.variant_type(), self.pos1, self.pos2))

    def __str__(self):
        return '{}, {}, {}, L:{}'.format(self.pos1, self.pos2, self.name, len(self))

    def __len__(self):
        return self.pos2.position() - self.pos1.position() - 1

    def _genome_region(self):
        return GenomeRegion(self.pos1.genome_pos, self.pos2.genome_pos)

    def matched_with_both_overlaps(self, other):
        if self.variant_type() != other.variant_type():
            return False
        return self.pos1.to_genome_region().overlaps(other.pos1.to_genome_region()) and \
            self.pos2.to_genome_region().overlaps(other.pos2.to_genome_region())

    def matched_with_reciprocal_overlap(self, other, threshold):
        if self.variant_type() != other.variant_type():
            return False
        return self._genome_region().reciprocal_overlaps(other._genome_region(), threshold)

    @abstractmethod
    def variant_type(self):
        pass


class Deletion(Variant):
    """docstring for Deletion"""

    def variant_type(self):
        return VariantTypes.DEL


class VariantFile(object):
    """docstring for VariantFile"""

    __metaclass__ = ABCMeta

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self._file_obj = self._open_file()
        return self

    def __exit__(self, *args):
        self._close_file()

    def __iter__(self):
        for line in self._file_obj:
            if self.is_valid_line(line):
                data = self.line_to_data(line)
                variant_type = self.get_variant_type(data)
                if variant_type == VariantTypes.DEL:
                    yield self.data_to_deletion(data)
                if variant_type == VariantTypes.INS:
                    yield self.data_to_insertion(data)
                if variant_type == VariantTypes.INV:
                    yield self.data_to_inversion(data)
                if variant_type == VariantTypes.DUP:
                    yield self.data_to_duplication(data)

    def _open_file(self):
        return open(self.filename, 'r')

    def _close_file(self):
        self._file_obj.close()

    @abstractmethod
    def get_variant_type(self, data):
        pass

    def is_valid_line(self, line):
        return True

    def line_to_data(self, line):
        return line.split()

    @abstractmethod
    def data_to_deletion(self, data):
        pass

    @abstractmethod
    def data_to_insertion(self, data):
        pass

    @abstractmethod
    def data_to_inversion(self, data):
        pass

    @abstractmethod
    def data_to_duplication(self, data):
        pass


class SvsimBedpeFile(VariantFile):
    """docstring for SvsimBedpeFile"""

    def __init__(self, filename, slop):
        super(SvsimBedpeFile, self).__init__(filename)
        self.slop = slop

    def get_variant_type(self, data):
        variant_types = ['DEL', 'INS', 'INV', 'DUP']
        pattern = '(?P<variant_type>{})'.format('|'.join(variant_types))
        m = re.search(pattern, data[6])
        return variant_types.index(m.group('variant_type'))

    def data_to_deletion(self, data):
        return Deletion(
            data[6],
            GenomePosition(data[0], int(data[2])).genome_position_with_ci(self.slop),
            GenomePosition(data[3], int(data[5])).genome_position_with_ci(self.slop)
        )

    def data_to_insertion(self, data):
        pass

    def data_to_inversion(self, data):
        self._file_obj.next()
        return Deletion(
            data[6],
            GenomePosition(data[0], int(data[2])),
            GenomePosition(data[3], int(data[5]))
        )

    def data_to_duplication(self, data):
        pass


class SpritesBedpeFile(VariantFile):
    """docstring for SpritesBedpeFile"""

    def get_variant_type(self, data):
        variant_types = ['DEL', 'INS', 'INV', 'DUP']
        pattern = '(?P<variant_type>{})'.format('|'.join(variant_types))
        m = re.search(pattern, data[6])
        return variant_types.index(m.group('variant_type'))

    def data_to_deletion(self, data):
        delta = int(data[2]) - int(data[1]) - 1
        if data[6].endswith('5F'):
            return Deletion(
                data[6],
                GenomePositionWithCi(GenomePosition(data[0], int(data[2])), Interval(-delta, 0)),
                GenomePositionWithCi(GenomePosition(data[3], int(data[5])), Interval(-delta, 0))
            )
        if data[6].endswith('5R'):
            return Deletion(
                data[6],
                GenomePositionWithCi(GenomePosition(data[0], int(data[1])+1), Interval(0, delta)),
                GenomePositionWithCi(GenomePosition(data[3], int(data[4])+1), Interval(0, delta))
            )
        return Deletion(
            data[6],
            GenomePositionWithCi(GenomePosition(data[0], int(data[2])), Interval(-delta, 0)),
            GenomePositionWithCi(GenomePosition(data[3], int(data[5])), Interval(-delta, 0))
        )

    def data_to_insertion(self, data):
        pass

    def data_to_inversion(self, data):
        pass

    def data_to_duplication(self, data):
        pass
