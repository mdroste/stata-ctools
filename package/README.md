# Editable Stata package sources

Edit `ado/`, `help/`, and `metadata/` here. `build/` remains the published flat
Stata net-install directory so existing installation URLs continue to work.
Its ado/help/manifest files are generated copies; never edit those copies.

```sh
python3 scripts/sync_release.py
python3 scripts/sync_package.py
python3 scripts/sync_release.py --check
python3 scripts/sync_package.py --check
```

`release.json` remains authoritative for version, date and minimum Stata version.
`make` synchronizes both sets of generated files before building. Plugin objects
and depfiles live in `build/obj/<platform>/<compiler-and-flags-digest>/`.

Stage a complete native package with `scripts/stage_package.py --revision REV`.
For all platforms, use `--platform all --binaries DIR`, where DIR contains flat
plugins or downloaded CI artifact directories. Local staging and CI call this
same CLI and validate the same manifest, licenses, platform scope and hashes.
The staged `BUILD_INFO.json` records the requested build revision and SHA-256 of
every payload. Publication still requires the licensed Stata correctness gate.
