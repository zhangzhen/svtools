from abc import ABCMeta, abstractmethod
from generic import GenomePosition, GenomeRegion


class TargetRegion(object):

    """docstring for TargetRegion"""

    __metaclass__ = ABCMeta

    def __init__(self, g_region, anchor):
        self.g_region = g_region
        self.anchor = anchor

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(self.g_region.reference_name(),
            self.g_region.start_pos(), self.g_region.end_pos(), self.anchor.pos)

    def reference_name(self):
        return self.g_region.reference_name()

    @abstractmethod
    def anchor_matches(self, variant):
        pass

    @abstractmethod
    def matches(self, variant):
        pass

    @classmethod
    def make_taregion(cls, ref_name, start_pos, end_pos, anchor_pos):
        return cls(GenomeRegion(GenomePosition(ref_name, start_pos),
                                GenomePosition(ref_name, end_pos)),
                   GenomePosition(ref_name, anchor_pos))


class TargetRegionForLeftAnchor(TargetRegion):

    """docstring for TargetRegionForLeftAnchor"""

    def __init__(self, g_region, anchor):
        assert g_region.reference_name() == anchor.reference_name(
        ) and g_region.end_pos() > anchor.position()
        super(TargetRegionForLeftAnchor, self).__init__(g_region, anchor)

    def anchor_matches(self, variant):
        return variant.pos1.matches(self.anchor)

    def matches(self, variant):
        return self.g_region.covers(variant.pos2)


class TargetRegionForRightAnchor(TargetRegion):

    """docstring for TargetRegionForRightAnchor"""

    def __init__(self, g_region, anchor):
        assert g_region.reference_name() == anchor.reference_name(
        ) and g_region.start_pos() < anchor.position()
        super(TargetRegionForRightAnchor, self).__init__(g_region, anchor)

    def anchor_matches(self, variant):
        return variant.pos2.matches(self.anchor)

    def matches(self, variant):
        return self.g_region.covers(variant.pos1)


def get_target_regions_from_file(filename):
    with open(filename, 'r') as f:
        for line in f:
            data = line.split()
            ref_name = data[0]
            start_pos = int(data[1])
            end_pos = int(data[2])
            anchor_pos = int(data[3])
            if end_pos > anchor_pos:
                yield TargetRegionForLeftAnchor.make_taregion(ref_name, start_pos, end_pos, anchor_pos)
            if start_pos < anchor_pos:
                yield TargetRegionForRightAnchor.make_taregion(ref_name, start_pos, end_pos, anchor_pos)


def eval_target_regions(taregions, truth_variants):
    taregions = list(taregions)
    num_taregions = len(taregions)
    num_anchor_matched = 0
    num_matched = 0

    unmatched_taregions = []

    i = j = 0

    while i < len(taregions) and j < len(truth_variants):
        if taregions[i].anchor_matches(truth_variants[j]):

            i += 1
            j += 1

    for ta in taregions:
        for v in truth_variants:
            if ta.anchor_matches(v):
                num_anchor_matched += 1
                if ta.matches(v):
                    num_matched += 1
                else:
                    unmatched_taregions.append(ta)

    precision = float(num_matched) / num_anchor_matched
    print "{}\t{}\t{}".format(num_taregions, num_anchor_matched, num_matched)
    print "Precision(%): {}".format(round(precision*100, 2))
    print "\n".join([str(x) for x in unmatched_taregions])