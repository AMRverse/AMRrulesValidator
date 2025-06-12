import argparse


def main():
    parser = argparse.ArgumentParser(prog="amrrules", description="AMR rules validator")
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
    
    args = parser.parse_args()
    
    if args.command == "validate":
        print("TODO validate")
        return 0
    elif args.command == "clean":
        print("TODO clean")
        return 0
    elif args.command == "update-resources":
        print("TODO update-resources")
        return 0
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    main()
