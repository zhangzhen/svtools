"""
Tests for `svtools` module.
"""
import mock
from svtools.variants import Deletion
from StringIO import StringIO
from svtools.targetregions import TargetRegionForLeftAnchor, TargetRegionForRightAnchor


class TestVariantFile(object):

    @classmethod
    def setup_class(cls):
        pass

    @mock.patch('svtools.variants.open')
    def test_svsim_bedpe_file(self, mock_open):
        slop = 3
        expected_d1 = variants.Deletion(
            'DEL0718::1::1',
            variants.GenomePosition('1', 576974).genome_position_with_ci(slop),
            variants.GenomePosition('1', 577575).genome_position_with_ci(slop)
        )
        mock_open.return_value = StringIO('1\t576973\t576974\t1\t577574\t577575\tDEL0718::1::1\t255\t+\t+')
        with variants.SvsimBedpeFile('any name', slop) as f:
            mock_open.assert_called_once_with('any name', 'r')
            d1 = iter(f).next()
            assert d1 == expected_d1

    @mock.patch('svtools.variants.open')
    def test_sprites_bedpe_file(self, mock_open):
        expected_d1 = variants.Deletion(
            'DEL.74.5R',
            variants.GenomePositionWithCi(variants.GenomePosition('11', 9637171), variants.Interval(0, 39)),
            variants.GenomePositionWithCi(variants.GenomePosition('11', 9637497), variants.Interval(0, 39))
        )
        expected_d2 = variants.Deletion(
            'DEL.88.5F',
            variants.GenomePositionWithCi(variants.GenomePosition('11', 67881097), variants.Interval(-67, 0)),
            variants.GenomePositionWithCi(variants.GenomePosition('11', 67881201), variants.Interval(-67, 0))
        )
        mock_open.return_value = StringIO(
            '\n'.join(
                [
                    '11\t9637170\t9637210\t11\t9637496\t9637536\tDEL.74.5R',
                    '11\t67881029\t67881097\t11\t67881133\t67881201\tDEL.88.5F'
                ]
            )
        )
        with variants.SpritesBedpeFile('any name') as f:
            mock_open.assert_called_once_with('any name', 'r')
            d1 = iter(f).next()
            assert d1 == expected_d1
            d2 = iter(f).next()
            assert d2 == expected_d2

    @classmethod
    def teardown_class(cls):
        pass


class TestVariant(object):

    @classmethod
    def setup_class(cls):
        cls.slop = 50

    def test_not_matched(self):
        d1 = variants.Deletion(
            'DELL:1877:11:9413301_9415177::11::3224',
            variants.GenomePosition('11', 9637171).genome_position_with_ci(self.slop),
            variants.GenomePosition('11', 9639049).genome_position_with_ci(self.slop)
        )
        d2 = variants.Deletion(
            'DEL.74.5R',
            variants.GenomePositionWithCi(variants.GenomePosition('11', 9637171), variants.Interval(0, 39)),
            variants.GenomePositionWithCi(variants.GenomePosition('11', 9637497), variants.Interval(0, 39))
        )

        assert not d1.matched_with_both_overlaps(d2)

    def test_matched(self):
        d1 = variants.Deletion(
            'DELL:10518:1:9907630_9918147::1::27',
            variants.GenomePosition('1', 10052764).genome_position_with_ci(self.slop),
            variants.GenomePosition('1', 10063283).genome_position_with_ci(self.slop)
        )
        d2 = variants.Deletion(
            'DEL.1.5F',
            variants.GenomePositionWithCi(variants.GenomePosition('1', 10052714), variants.Interval(0, 26)),
            variants.GenomePositionWithCi(variants.GenomePosition('1', 10063233), variants.Interval(0, 26))
        )

        assert d1.matched_with_both_overlaps(d2)

    @classmethod
    def teardown_class(cls):
        pass


class TestTargetRegion(object):

    @classmethod
    def setup_class(cls):
        cls.slop = 50

    def test_matches(self):
        taregion1 = TargetRegionForLeftAnchor.make_taregion(11, 4910142, 4910792, 4801135)
        variant = Deletion('', )
