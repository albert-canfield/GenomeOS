from genomeos.lib import LAYERS, LIBRARIES, by_layer


def test_every_layer_has_libraries():
    for layer in LAYERS:
        assert by_layer(layer), layer


def test_library_ids_match_layers_and_have_sources():
    for lib in LIBRARIES.values():
        assert lib.id.startswith(lib.layer + "."), lib.id
        assert lib.source, lib.id
        assert lib.purpose


def test_gene_symbols_look_like_hgnc():
    import re

    for lib in LIBRARIES.values():
        for g in lib.genes:
            assert re.fullmatch(r"[A-Z0-9][A-Z0-9-]*", g), (lib.id, g)
