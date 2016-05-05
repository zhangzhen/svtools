from abc import ABCMeta, abstractmethod
from generic import GenomePosition, GenomeRegion


class TargetRegion(object):
    """docstring for TargetRegion"""

    __metaclass__ = ABCMeta

    def __init__(self, g_region, anchor):
        self.g_region = g_region
        self.anchor = anchor

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
        assert g_region.ref_name() == anchor.ref_name and g_region.end_pos() > anchor.pos
        super(TargetRegionForLeftAnchor, self).__init__(g_region, anchor)

    def anchor_matches(self, variant):
        return variant.pos1.matches(self.anchor)

    def matches(self, variant):
        return self.anchor_matches(variant) and self.g_region.covers(variant.pos2)


class TargetRegionForRightAnchor(TargetRegion):
    """docstring for TargetRegionForRightAnchor"""
    def __init__(self, g_region, anchor):
        assert g_region.ref_name() == anchor.ref_name and g_region.start_pos() < anchor.pos
        super(TargetRegionForRightAnchor, self).__init__(g_region, anchor)

    def anchor_matches(self, variant):
        return variant.pos2.matches(self.anchor)

    def matches(self, variant):
        return self.anchor_matches(variant) and self.g_region.covers(variant.pos1)
