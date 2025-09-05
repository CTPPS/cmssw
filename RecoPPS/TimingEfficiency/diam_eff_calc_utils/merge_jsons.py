import json
import argparse
import os


def merge_json(input_path, output_path):
    if os.path.exists(output_path):
        with open(output_path, "r", encoding="utf-8") as f:
            output_data = json.load(f)
    else:
        output_data = {}

    with open(input_path, "r", encoding="utf-8") as f:
        input_data = json.load(f)

    if isinstance(input_data, dict) and isinstance(output_data, dict):
        output_data.update(input_data)
    else:
        raise ValueError("Both JSON files must be JSON objects (dictionaries)")

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=4, ensure_ascii=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge input JSON into the output JSON"
    )
    parser.add_argument(
        "--output_json",
        dest="output_json",
        help="Path to the output JSON file (will be updated)",
    )
    parser.add_argument(
        "--input_json",
        dest="input_json",
        help="Path to the input JSON file to copy content from",
    )
    args = parser.parse_args()

    merge_json(args.input_json, args.output_json)
