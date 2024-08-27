# File: pyoz_input.py

import sys
import argparse
import toml
from typing import Dict, Any, Tuple, Optional


class InputParser:
    def __init__(self):
        self.parser = self._create_parser()

    def _create_parser(self) -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(description="pyOZ - Ornstein-Zernike equation solver")
        parser.add_argument("-i", "--input", required=True, help="Input configuration file")
        parser.add_argument("-o", "--output", help="Output file for redirecting stdout")
        parser.add_argument("-g", "--gamma", help="Initial gamma function file")
        parser.add_argument(
            "-b", "--binarygamma", action="store_true", help="Gamma file is in binary format"
        )
        return parser

    def parse_cmdline(self, argv: Optional[list] = None) -> Dict[str, Any]:
        if argv is None:
            argv = sys.argv[1:]
        args = self.parser.parse_args(argv)
        return vars(args)

    def parse_input(
        self, cmdline: Dict[str, Any]
    ) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any], Any]:
        with open(cmdline["input"], "r") as f:
            config = toml.load(f)

        ctrl = config.get("control", {})
        syst = config.get("system", {})
        parm = config.get("parameters", {})
        outp = config.get("output", {})
        const = self._create_constants(config.get("constants", {}))

        return ctrl, syst, parm, outp, const

    def _create_constants(self, const_config: Dict[str, Any]) -> Any:
        class Constants:
            def __init__(self, **kwargs):
                for key, value in kwargs.items():
                    setattr(self, key, value)

        return Constants(**const_config)


def parse_cmdline(argv: Optional[list] = None) -> Dict[str, Any]:
    parser = InputParser()
    return parser.parse_cmdline(argv)


def parse_input(
    cmdline: Dict[str, Any]
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any], Any]:
    parser = InputParser()
    return parser.parse_input(cmdline)


# Usage example:
# cmdline = parse_cmdline()
# ctrl, syst, parm, outp, const = parse_input(cmdline)
