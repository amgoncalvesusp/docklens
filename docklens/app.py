"""
app.py — entry point for the DockLens desktop application.

Run with:  python -m docklens.app
"""

from __future__ import annotations

import logging
import sys


LOGGER = logging.getLogger(__name__)


def main():
    if "--self-check" in sys.argv:
        from .self_check import run_self_check

        return run_self_check()
    if "--check-manifest" in sys.argv:
        try:
            index = sys.argv.index("--check-manifest")
            manifest_path = sys.argv[index + 1]
        except (ValueError, IndexError):
            return 2
        import json

        from . import batch_runner
        from .integration_manifest import ManifestError, load_manifest
        from .integration_result import write_integration_result

        try:
            manifest = load_manifest(manifest_path)
            result = batch_runner.run_paired(
                manifest.receptor_path,
                manifest.poses_path,
                key_residues=manifest.key_residues,
                hbond_preset=manifest.hbond_preset,
            )
            if getattr(manifest, "result_path", None) is not None:
                write_integration_result(manifest, result)
        except Exception as error:  # noqa: BLE001 - process boundary contains parser failures
            code = error.code if isinstance(error, ManifestError) else "analysis_failed"
            if not isinstance(error, (ManifestError, OSError, ValueError)):
                LOGGER.exception("Unexpected DockingHub manifest analysis failure")
            print("DockLens manifest error: %s" % code, file=sys.stderr)
            return 2
        print(json.dumps({
            "schema": "docklens-manifest-check-v1",
            "poses": len(result.summaries),
            "interactions": len(result.details),
            "errors": sum(record.status == "error" for record in result.input_qc),
        }, sort_keys=True))
        return 0
    if "--manifest" in sys.argv:
        try:
            index = sys.argv.index("--manifest")
            manifest_path = sys.argv[index + 1]
        except (ValueError, IndexError):
            return 2
        from .integration_manifest import ManifestError, load_manifest

        try:
            manifest = load_manifest(manifest_path)
        except Exception as error:  # noqa: BLE001 - process boundary contains manifest failures
            code = error.code if isinstance(error, ManifestError) else "manifest_failed"
            if not isinstance(error, ManifestError):
                LOGGER.exception("Unexpected DockingHub manifest load failure")
            print("DockLens manifest error: %s" % code, file=sys.stderr)
            return 2
        from .main_window import launch

        return launch(launch_manifest=manifest)
    from .main_window import launch

    return launch()


if __name__ == "__main__":
    sys.exit(main())
