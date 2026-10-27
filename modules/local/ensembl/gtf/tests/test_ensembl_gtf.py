import importlib.util
from pathlib import Path


_SCRIPT = Path(__file__).resolve().parents[5] / "bin" / "ensembl_gtf.py"
spec = importlib.util.spec_from_file_location("ensembl_gtf", _SCRIPT)
eg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eg)


def test_release_path_current_uses_main_ensembl_current_gtf_layout():
    assert eg.release_path_for_value("current") == "current/gtf"


def test_release_path_latest_uses_main_ensembl_current_gtf_layout():
    assert eg.release_path_for_value("latest") == "current/gtf"


def test_release_path_numbered_release_uses_release_gtf_layout():
    assert eg.release_path_for_value("115") == "release-115/gtf"


def test_release_path_prefixed_release_uses_release_gtf_layout():
    assert eg.release_path_for_value("release-115") == "release-115/gtf"


def test_find_species_searches_ensemblgenomes_collections():
    from unittest.mock import patch
    fungi_root = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/current/gtf/"
    species_dir = f"{fungi_root}fungi_ascomycota1_collection/debaryomyces_hansenii_cbs767_gca_000006445/"
    gtf_file = "Debaryomyces_hansenii_cbs767_gca_000006445.ASM644v2.63.gtf.gz"
    listing = lambda *names: "".join(f'<a href="{n}">{n}</a>\n' for n in names)

    def fake_fetch(url, timeout=60):
        if url == fungi_root:
            return listing("fungi_ascomycota1_collection/")
        if url == f"{fungi_root}fungi_ascomycota1_collection/":
            return listing("debaryomyces_hansenii_cbs767_gca_000006445/", "candida_albicans/")
        if url == species_dir:
            return listing(gtf_file)
        raise eg.EnsemblSpeciesNotFound(f"404 {url}")

    with patch.object(eg, "fetch_text", side_effect=fake_fetch):
        root, name, division = eg._find_species("https://ftp.ensembl.org/pub", "current", "debaryomyces_hansenii")

    assert (root, name, division) == (species_dir, gtf_file, "EnsemblGenomes/fungi")
