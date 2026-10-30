"""Stage commands for the native scPerturb workflow."""
import runpy
import sys
COMMANDS = {"preflight": "stages.preflight", "dragen": "stages.native_crispr", "guide-qc": "stages.guide_qc", "reanalysis": "stages.reanalysis", "report": "stages.report"}
def main():
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        print("Usage: python -m scperturb {" + ",".join(COMMANDS) + "} [stage arguments]")
        return 0
    command = sys.argv.pop(1)
    if command not in COMMANDS:
        raise SystemExit("Unknown stage: " + command)
    runpy.run_module("scperturb." + COMMANDS[command], run_name="__main__")
    return 0
if __name__ == "__main__":
    raise SystemExit(main())
