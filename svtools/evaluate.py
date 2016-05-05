def _find_matches(truth_variants, predicted_variants, matcher):
    return ((truth, predicted) for truth in truth_variants for predicted in predicted_variants
            if matcher(truth, predicted))


def _find_no_matches(truth_variants, predicted_variants, matcher):
    a = set(predicted_variants)
    b = set((v2 for v1, v2 in _find_matches(truth_variants, predicted_variants, matcher)))
    return list(a - b)


def get_variants(cls, filename, *args):
    with cls(filename, *args) as f:
        return [v for v in iter(f)]


def _show_difference(variants1, variants2):
    a = set(variants1)
    b = set(variants2)
    increased = b - a
    reduced = a - b
    print '{} vs. {} ({} increased, {} removed)'.format(len(a), len(b), len(increased), len(reduced))
    print '\n'.join(('+ {}'.format(v) for v in sorted(list(increased))))
    print '\n'.join(('- {}'.format(v) for v in sorted(list(reduced))))


def compare_variants(truth, pred1, pred2, matcher):
    matches1 = _find_matches(truth, pred1, matcher)
    matches2 = _find_matches(truth, pred2, matcher)
    _show_difference(
        (a for a, b in matches1),
        (a for a, b in matches2)
    )
    print '\n'
    no_matches1 = _find_no_matches(truth, pred1, matcher)
    no_matches2 = _find_no_matches(truth, pred2, matcher)
    _show_difference(no_matches1, no_matches2)
