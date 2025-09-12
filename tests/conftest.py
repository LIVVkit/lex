from livvkit import __main__
import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def generate_livv_output():
    outdir = "simple_extn_output"
    __main__.main(["-V", "tests/simple_test.yml", "-o", outdir])
    return Path(outdir, "index.json")