from pathlib import Path

import pytest
from livvkit import __main__


@pytest.fixture(scope="session")
def generate_livv_output():
    outdir = "simple_extn_output"
    __main__.main(["-V", "tests/simple_test.yml", "-o", outdir])
    return Path(outdir, "index.json")
