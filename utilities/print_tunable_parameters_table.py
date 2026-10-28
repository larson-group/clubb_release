#!/usr/bin/env python3
"""Compare model parameter JSON definitions as a Markdown table.

Usage:
    python utilities/print_tunable_parameters_table.py input/clubb_params.jsonc
    python utilities/print_tunable_parameters_table.py first.json second.json -bare -output_file /tmp/parameters.md

main() reads each catalog's selected -config (default: default), or a complete
model/SILHS snapshot. -bare shows only differing parameters. -format writes a
normalized resolved snapshot beside each source with a -std.json suffix. Column
arrays remain arrays; absent parameters are shown as missing.
"""
import argparse
import json
from pathlib import Path
import sys

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from input.scripts.parameters_json_to_namelist import parameter_settings
from input.scripts.namelist import load_json


def main():
    """Render the selected JSON parameter values and optionally write normalized snapshots."""
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False, add_help=False)
    parser.add_argument("-h", "-help", action="help")
    parser.add_argument("infiles", nargs="+", type=Path)
    parser.add_argument("-config", default="default")
    parser.add_argument("-output_file", type=Path)
    parser.add_argument("-bare", action="store_true")
    parser.add_argument("-format", action="store_true")
    args = parser.parse_args()
    try:
        settings = [parameter_settings(args.config, path) if "presets" in load_json(path) else parameter_settings(path, config=args.config) for path in args.infiles]
        models = [values["model"] for values in settings]
        names = dict.fromkeys(name for values in models for name in values)
        rows = ["| Parameter | " + " | ".join(str(path) for path in args.infiles) + " |",
                "| --- | " + " | ".join("---" for _ in args.infiles) + " |"]
        for name in names:
            values = [model.get(name) for model in models]
            if args.bare and all(value == values[0] for value in values):
                continue
            rows.append("| " + name + " | " + " | ".join("missing" if value is None else json.dumps(value) for value in values) + " |")
        text = "\n".join(rows) + "\n"
        if args.output_file:
            args.output_file.write_text(text, encoding="utf-8")
        else:
            print(text, end="")
        if args.format:
            for path, values in zip(args.infiles, settings):
                path.with_name(path.stem + "-std.json").write_text(json.dumps(values, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError) as exc:
        parser.exit(2, f"{parser.prog}: error: {exc}\n")


if __name__ == "__main__":
    main()
