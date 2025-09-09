import lex.convert_cfg as lconvert
import livvkit.util.functions as fcn
from pathlib import Path
import os


def test_convert_to_yaml():
    in_file = "json_to_convert.json"
    out_file = "yml_to_convert.yml"
    ref_file = "yml_reference.yml"

    assert not Path(out_file).exists()
    lconvert.json_to_yaml(Path(in_file))
    assert Path(out_file).exists()
    ref_yml = fcn.read_yaml(Path(ref_file))
    test_yml = fcn.read_yaml(Path(out_file))
    assert ref_yml == test_yml
    os.remove(out_file)
    assert not Path(out_file).exists()
