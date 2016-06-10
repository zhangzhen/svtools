from abc import ABCMeta, abstractmethod
from generic import ChromNameIdConverter, GenomePosition, GenomeRegion


class TargetRegion(object):

    """docstring for TargetRegion"""

    __metaclass__ = ABCMeta

    def __init__(self, g_region, anchor):
        self.g_region = g_region
        self.anchor = anchor

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(
            self.g_region.reference_name(),
            self.g_region.start_pos(),
            self.g_region.end_pos(),
            self.anchor.pos
        )

    def reference_name(self):
        return self.g_region.reference_name()

    @abstractmethod
    def anchor_matches(self, variant):
        pass

    @abstractmethod
    def anchor_before(self, variant):
        pass

    @abstractmethod
    def after(self, variant):
        pass

    @abstractmethod
    def covers(self, variant):
        pass

    @classmethod
    def make_taregion(cls, ref_name, start_pos, end_pos, anchor_pos):
        return cls(
            GenomeRegion(
                GenomePosition(ChromNameIdConverter.name_to_id(ref_name), ref_name, start_pos),
                GenomePosition(ChromNameIdConverter.name_to_id(ref_name), ref_name, end_pos)
            ),
            GenomePosition(ChromNameIdConverter.name_to_id(ref_name), ref_name, anchor_pos)
        )


class TargetRegionForLeftAnchor(TargetRegion):

    """docstring for TargetRegionForLeftAnchor"""

    def __init__(self, g_region, anchor):
        assert g_region.reference_name() == anchor.reference_name(
        ) and g_region.end_pos() > anchor.position()
        super(TargetRegionForLeftAnchor, self).__init__(g_region, anchor)

    def anchor_before(self, variant):
        return variant.pos1.after(self.anchor)

    def anchor_matches(self, variant):
        return variant.pos1.matches(self.anchor)

    def after(self, variant):
        return self.g_region.after(variant.pos2)

    def covers(self, variant):
        return self.g_region.covers(variant.pos2)


class TargetRegionForRightAnchor(TargetRegion):

    """docstring for TargetRegionForRightAnchor"""

    def __init__(self, g_region, anchor):
        assert g_region.reference_name() == anchor.reference_name(
        ) and g_region.start_pos() < anchor.position()
        super(TargetRegionForRightAnchor, self).__init__(g_region, anchor)

    def anchor_before(self, variant):
        return variant.pos2.after(self.anchor)

    def anchor_matches(self, variant):
        return variant.pos2.matches(self.anchor)

    def after(self, variant):
        return self.g_region.after(variant.pos1)

    def covers(self, variant):
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

    num_taregions = len(taregions)
    num_truth_varirants = len(truth_variants)

    num_anchor_matched = 0
    num_covered = 0

    uncovered_taregions = []
    matched_truth_variant_names = set()

    for ta in taregions:
        for v in truth_variants:
            if ta.anchor_matches(v):
                num_anchor_matched += 1
                if ta.covers(v):
                    matched_truth_variant_names.add(v.name)
                    num_covered += 1
                else:
                    uncovered_taregions.append(ta)
                break

    num_matched_truth_variants = len(matched_truth_variant_names)

    recall = float(num_matched_truth_variants) / num_truth_varirants
    precision = float(num_covered) / num_anchor_matched

    # print "#truth\t#taregions\t#matched\t#covered\trecall\tprecision"
    print "{}\t{}\t{}%\t{}\t{}\t{}\t{}%".format(
        num_truth_varirants,
        num_matched_truth_variants,
        round(recall*100, 2),
        num_taregions,
        num_anchor_matched,
        num_covered,
        round(precision*100, 2)
    )

    print "\nUncovered target regions whose anchor is matched:"
    print "\n".join([str(x) for x in uncovered_taregions])
