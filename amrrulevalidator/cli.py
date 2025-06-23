import argparse
from pathlib import Path
from amrrulevalidator.utils.resources import ResourceManager
from amrrulevalidator.validate import run_validate
from amrrulevalidator.clean import run_clean
from amrrulevalidator.convert_rules import run_convert_to_latest_spec


def main():
    parser = argparse.ArgumentParser(prog="amrrule", description="AMRrulevalidator")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # validate subcommand
    validate_parser = subparsers.add_parser("validate", help="Validate draft AMRrules file")
    validate_parser.add_argument("--input", required=True, help="AMRrules tsv file to validate")
    validate_parser.add_argument("--output", required=True, help="Validated AMRrules file, in tsv format, annotated with values to check")
    
    # clean subcommand
    clean_parser = subparsers.add_parser("clean", help="Produces a clean AMRrules file for import into the interpretation engine")
    clean_parser.add_argument("--input", required=True, help="AMRrules tsv file to clean (should have been validated first)")
    clean_parser.add_argument("--output", required=True, help="Cleaned AMRrules file, in tsv format")
    
    # update-resources subcommand
    update_parser = subparsers.add_parser("update-resources", help="Update CARD and AMRFinderPlus resources used for validation")

    # convert to latest spec subcommand
    convert_parser = subparsers.add_parser("convert-to-latest-spec", help="Convert AMRrules file to the latest specification")
    convert_parser.add_argument("--input", required=True, help="AMRrules tsv file to convert")
    convert_parser.add_argument("--output", required=True, help="Converted AMRrules file, in tsv format")
    
    args = parser.parse_args()
    
    if args.command == "validate":
        print(f"Validating {args.input}...")
        input_path = Path(args.input)
        output_path = Path(args.output)
        resource_manager = ResourceManager()
        success = run_validate(input_path, output_path, resource_manager)
        if success:
            print(f"Validation complete. Output written to {args.output}")
            return 0
        else:
            print("Validation failed.")
            return 1
    elif args.command == "clean":
        print(f"Cleaning {args.input}...")
        input_path = Path(args.input)
        output_path = Path(args.output)
        success = run_clean(input_path, output_path)
        if success:
            print(f"Cleaning complete. Output written to {args.output}")
            return 0
        else:
            print("Cleaning failed.")
            return 1
    elif args.command == "update-resources":
        print("Updating resources...")
        resource_manager = ResourceManager()
        success = resource_manager.setup_all_resources()
        if success:
            print("Resources updated successfully.")
            return 0
        else:
            print("Failed to update some resources.")
            return 1
    elif args.command == "convert-to-latest-spec":
        input_path = Path(args.input)
        output_path = Path(args.output)
        resource_manager = ResourceManager()
        success = run_convert_to_latest_spec(input_path, output_path, resource_manager)
        if success:
            return 0
        else:
            print("Conversion failed.")
            return 1
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    main()
