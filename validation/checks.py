"""Checks: executable checks and their negative controls."""

import fnmatch
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from scripts import stage_package as packaging
from scripts.package_contract import check_package
from scripts.stage_package import stage
from validation import ROOT, stata
from validation.native import unit
from validation.stata import stata_command


def release_workflow_contract():
    source = (ROOT / ".github/workflows/build.yml").read_text()
    for job in ("commit-plugins", "release"):
        block = source.split("\n  " + job + ":\n")[-1]
        condition = re.search(r"    if: (.*?)\n    needs:", block, re.S).group(1)
        condition = condition.removeprefix(">-").strip()
        assert "needs.package.result" in condition
        for event in ("push", "pull_request", "workflow_dispatch", "release"):
            for ref in ("refs/heads/main", "refs/heads/dev", "refs/tags/v1.0.0"):
                for requested in (False, True):
                    for status in ("success", "failure", "skipped", "cancelled"):
                        for cancelled in (False, True):
                            expression = condition.replace("always()", "True")
                            expression = expression.replace("!cancelled()", repr(not cancelled))
                            for key, value in {
                                "needs.package.result": status,
                                "github.event_name": event,
                                "github.ref": ref,
                                "inputs.publish_plugins": requested,
                            }.items():
                                expression = expression.replace(key, repr(value))
                            expression = " ".join(
                                expression.replace("&&", "and").replace("||", "or").split()
                            )
                            publish = (
                                event == "release"
                                if job == "release"
                                else event == "workflow_dispatch"
                                and requested
                                and ref == "refs/heads/main"
                            )
                            allowed = publish and status == "success" and not cancelled
                            assert eval(expression, {"__builtins__": {}}, {}) == allowed, (
                                job,
                                event,
                                ref,
                                requested,
                                status,
                                cancelled,
                            )
    dispatch = source.split("  workflow_dispatch:\n")[1].split("  release:\n")[0]
    assert "publish_plugins:" in dispatch
    assert "type: boolean" in dispatch and "default: false" in dispatch
    paths = re.findall(r"    paths:\n((?:      - .*\n)+)", source)
    assert len(paths) == 2
    for block in paths:
        patterns = re.findall(r"'([^']+)'", block)
        for name in ("build/new.ado", "build/cimport.sthlp", "build/ctools.pkg"):
            assert any(fnmatch.fnmatchcase(name, pattern) for pattern in patterns)
    assert "stata-validation:" not in source
    assert "CTOOLS_STATA_RUNNERS" not in source
    for platform in ("linux", "mac_arm", "mac_x86", "windows"):
        assert f"name: ctools_{platform}" in source
    assert "needs.sign-windows.result == 'success'" in source
    # A release uploads only under a tag that names config/release.json's version.
    from scripts.sync_release import release_tag_error

    version = json.loads((ROOT / "config/release.json").read_text())["version"]
    assert release_tag_error("v" + version) is None
    for tag in (version, "v0.0.0", f"v{version}-rc1", f"V{version}", ""):
        assert release_tag_error(tag), tag
    release_block = source.split("\n  release:\n")[-1]
    tag_step = release_block.index('sync_release.py --check --tag "$RELEASE_TAG"')
    assert tag_step < release_block.index("archive_package.py")
    assert "RELEASE_TAG: ${{ github.event.release.tag_name }}" in release_block
    # The Linux plugin is compiled and load-tested on the manylinux2014 floor.
    linux_block = source.split("\n  build-linux:\n")[1].split("\n  build-macos-arm:\n")[0]
    assert "make linux-docker BUILD_REVISION=" in linux_block
    assert "make linux-docker-load" in linux_block
    assert "run: make linux BUILD_REVISION" not in linux_block
    with patch.dict(os.environ, {"CTOOLS_STATA_EXECUTABLE": "/test path/Stata"}):
        assert stata_command("driver.do")[0] == "/test path/Stata"
        assert stata_command("driver.do")[-1] == str(Path("driver.do").resolve())
    with patch(
        "validation.stata.os",
        SimpleNamespace(
            name="nt", environ={"CTOOLS_STATA_EXECUTABLE": "C:/Program Files/Stata/StataMP.exe"}
        ),
    ):
        command = stata_command("driver.do")
        assert command[0] == "C:/Program Files/Stata/StataMP.exe"
        assert command[1:3] == ["/e", "do"]
    print("PASS publication requires explicit action and successful checks; pushes only build")


def release_packages():
    with tempfile.TemporaryDirectory(prefix="ctools-release-integrity-") as tmp:
        directory = stage(ROOT / "build", Path(tmp) / "package", ["mac_arm"], "test-only")
        assert not check_package(directory, integrity=True)
        info_path = directory / "BUILD_INFO.json"
        original = info_path.read_bytes()
        for name in ("ctools_mac_arm.plugin", "cppmlhdfe_p.ado"):
            payload_path = directory / name
            payload = payload_path.read_bytes()
            payload_path.write_bytes(b"corrupt or truncated payload")
            assert any("checksum mismatch" in e for e in check_package(directory, integrity=True))
            payload_path.write_bytes(payload)
        for mutation in ("missing", "scope", "inventory", "digest", "manifest", "invalid"):
            info = json.loads(original)
            if mutation == "missing":
                info_path.unlink()
            elif mutation == "invalid":
                info_path.write_text("null")
            else:
                if mutation == "scope":
                    info["platforms"] = ["linux"]
                if mutation == "inventory":
                    info["sha256"].pop("cppmlhdfe_p.ado")
                if mutation == "digest":
                    info["sha256"]["ctools.pkg"] = "0" * 64
                if mutation == "manifest":
                    (directory / "ctools.pkg").write_text(
                        (directory / "ctools.pkg").read_text() + "d changed\n"
                    )
                info_path.write_text(json.dumps(info))
            assert check_package(directory, integrity=True), mutation
            if mutation == "manifest":
                p = directory / "ctools.pkg"
                p.write_text(p.read_text().removesuffix("d changed\n"))
            info_path.write_bytes(original)
        assert not check_package(directory, integrity=True)
    print("PASS checksum corruption, inventory, scope, manifest, and missing build-info checks")


def release_cleanup():
    with tempfile.TemporaryDirectory(prefix="ctools-clean-") as tmp:
        directory = Path(tmp)
        build, source = directory / "build", directory / "src"
        build.mkdir()
        source.mkdir()
        (build / "obj").mkdir()
        artifacts = [
            build / "ctools.plugin",
            build / "ctools_windows.plugin",
            build / "x.o",
            source / "x.o",
            build / "obj/x.o",
        ]
        for artifact in artifacts:
            artifact.touch()
        result = subprocess.run(
            [
                "make",
                "-s",
                "-f",
                str(ROOT / "Makefile"),
                "clean",
                "DETECTED_OS=Windows",
                f"BUILD_DIR={build}",
                f"SRC_DIR={source}",
                "LIBOMP_PREFIX=/nonexistent",
            ],
            cwd=ROOT,
            check=True,
        )
        assert all(not p.exists() for p in artifacts)
        # A matching directory cannot be unlinked: genuine cleanup errors propagate.
        (build / "bad.o").mkdir()
        result = subprocess.run(
            ["python3", "scripts/clean_build.py", str(build), str(source)],
            cwd=ROOT,
            capture_output=True,
        )
        assert result.returncode != 0
    print("PASS MSYS/POSIX recipe cleanup, generic plugin, and error propagation")


def package_must_reject(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError("invalid release input was accepted")


def package_packages():
    with tempfile.TemporaryDirectory(prefix="ctools-package-test-") as tmp:
        directory = Path(tmp)
        source = directory / "source"
        source.mkdir()
        # Real source manifests and helpers; dummy binaries test staging only.
        for entry in (ROOT / "build").iterdir():
            if entry.is_file() and (
                entry.suffix in (".ado", ".sthlp", ".pkg")
                or entry.name in ("stata.toc", "LICENSE", "THIRD_PARTY_NOTICES")
            ):
                shutil.copy2(entry, source / entry.name)
        # Every other non-plugin file the manifest lists (e.g. codec licences)
        for line in (ROOT / "build" / "ctools.pkg").read_text().splitlines():
            name = line[2:]
            if line.startswith("f ") and not name.endswith(".plugin"):
                shutil.copy2(ROOT / "build" / name, source / name)
        for target in packaging.PLATFORMS:
            (source / f"ctools_{target}.plugin").write_bytes(b"test fixture, not a plugin")
        for targets in (["mac_arm"], list(packaging.PLATFORMS)):
            staged = packaging.stage(
                source,
                directory / ("stage-" + targets[0] + str(len(targets))),
                targets,
                "test-only",
            )
            assert not packaging.check_package(staged)
            assert (staged / "cppmlhdfe_p.ado").is_file()
            assert {p.name for p in staged.glob("*.plugin")} == {
                f"ctools_{p}.plugin" for p in targets
            }
            (staged / "cppmlhdfe_p.ado").unlink()
            assert packaging.check_package(staged)
        (source / "ctools_linux.plugin").unlink()
        package_must_reject(
            lambda: packaging.stage(
                source, directory / "bad", list(packaging.PLATFORMS), "test-only"
            )
        )
        assert not (directory / "bad").exists()
        (source / "cppmlhdfe_p.ado").unlink()
        package_must_reject(
            lambda: packaging.stage(source, directory / "bad-helper", ["mac_arm"], "test-only")
        )
    print("PASS single/all-platform manifests, helper completeness, and atomic staging")


def package_inventory():
    public = {p.stem for p in (ROOT / "build").glob("*.sthlp")} - {"ctools", "ctools_licenses"}
    registered = {p.stem for p in (ROOT / "validation/stata").glob("c*.do")} - {"contracts"}
    available, _ = stata.cases()
    assert public <= available.keys()
    assert public == registered, (public - registered, registered - public)
    manifest = (ROOT / "build/ctools.pkg").read_text().splitlines()
    for name in public:
        assert (ROOT / "src" / name).is_dir(), name
        for suffix in (".ado", ".sthlp"):
            assert (ROOT / "build" / (name + suffix)).is_file(), name
            assert f"f {name}{suffix}" in manifest, name
    print(f"PASS {len(public)} public-command help/source/package/test inventory")


# Mata and Stata functions added after Stata 14.1, the package minimum.
package_STATA16 = (
    "pathresolve", "pathgetparent", "issamefile", "isascii", "expm1", "ln1p", "ln1m",
    "ustrsplit", "st_frame",
)
package_OBS_LIMIT = "2147483647"


def package_stata_contract():
    """Ado code stays within Stata 14.1, guards the plugin's observation limit,
    and carries the same ado/plugin protocol as the plugin source."""
    newer = re.compile(r"\b(?:" + "|".join(package_STATA16) + r")\w*\s*\(")
    assert newer.search('mata: st_local("d", pathgetparent(st_local("using")))')
    assert newer.search("mata: st_framecurrent()") and not newer.search("pathbasename(p)")
    errors = []
    for path in sorted((ROOT / "build").glob("*.ado")):
        for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            code = line.split("//")[0]
            if not code.lstrip().startswith("*") and newer.search(code):
                errors.append(f"{path.name}:{number}: needs Stata 16: {code.strip()}")
    # Every public command stops at Stata's plugin limit of 2^31-1 observations;
    # cimport checks the row count it scanned.
    public = {p.stem for p in (ROOT / "build").glob("*.sthlp")} - {"ctools", "ctools_licenses"}
    for name in sorted(public):
        source = (ROOT / "build" / f"{name}.ado").read_text()
        guard = "`nobs' > " if name == "cimport" else "_N > "
        if guard + package_OBS_LIMIT not in source:
            errors.append(f"{name}.ado: no {guard}{package_OBS_LIMIT} observation guard")
    release = json.loads((ROOT / "config/release.json").read_text())
    protocol = str(release["protocol"])
    stamped = {
        "src/ctools_version.h": re.findall(
            r'#define CTOOLS_PROTOCOL "(\d+)"', (ROOT / "src/ctools_version.h").read_text()
        ),
        "build/_ctools_load.ado": re.findall(
            r'local __ado_protocol "(\d+)"', (ROOT / "build/_ctools_load.ado").read_text()
        ),
    }
    for name, values in stamped.items():
        if values != [protocol]:
            errors.append(f"{name}: protocol {values}, config/release.json has {protocol}")
    plugin = (ROOT / "src/ctools_plugin.c").read_text()
    if 'SF_macro_save("_ctools_plugin_protocol", CTOOLS_PROTOCOL)' not in plugin:
        errors.append("src/ctools_plugin.c: the version call does not report CTOOLS_PROTOCOL")
    loader = (ROOT / "build/_ctools_load.ado").read_text()
    if '"`ctools_plugin_protocol\'" != "`__ado_protocol\'"' not in loader:
        errors.append("build/_ctools_load.ado: the plugin protocol is not compared")
    assert not errors, "\n".join(errors)
    print(f"PASS Stata 14.1 functions, {len(public)} observation guards, protocol {protocol}")


def package_licenses():
    """Every bundled library's license is in the notices, and net install
    installs both texts as ctools_licenses.sthlp."""

    def norm(text):
        return " ".join(text.split())

    notices = norm((ROOT / "docs/licenses/THIRD_PARTY_NOTICES").read_text())
    vendor = sorted(
        p for p in (ROOT / "src/vendor").glob("*/*") if p.name in ("LICENSE", "COPYING")
    )
    assert vendor
    missing = [str(p.relative_to(ROOT)) for p in vendor if norm(p.read_text()) not in notices]
    assert not missing, f"THIRD_PARTY_NOTICES lacks {missing}"
    installed = norm((ROOT / "build/ctools_licenses.sthlp").read_text())
    for name in ("LICENSE", "docs/licenses/THIRD_PARTY_NOTICES"):
        assert norm((ROOT / name).read_text()) in installed, name
    assert "{asis} " + notices[:200] in installed
    manifest = (ROOT / "build/ctools.pkg").read_text().splitlines()
    assert "f ctools_licenses.sthlp" in manifest
    print(f"PASS {len(vendor)} vendor licenses in the notices and the installed help file")


def package_description():
    """The package description lists every public command."""
    listed = set()
    for line in (ROOT / "build/ctools.pkg").read_text().splitlines():
        match = re.match(r"d   (\S+)\s+- ", line)
        if match:
            listed.update(match[1].split("/"))
    public = {p.stem for p in (ROOT / "build").glob("*.sthlp")} - {"ctools", "ctools_licenses"}
    missing = sorted(public - listed)
    assert not missing, f"ctools.pkg description omits {missing}"
    print(f"PASS package description lists {len(public)} commands")


docs_BENCHMARKS = ("docs/benchmarks/hdfe_datasets_", "docs/benchmarks/readme_profile/")
docs_DATED = re.compile(r"(?<!\d)(19|20)\d{2}-?[01]\d-?[0-3]\d(?!\d)")
docs_ARTIFACTS = {".do", ".py", ".plugin", ".log", ".patch", ".dta", ".o"}


def docs_main(argv=None):
    files = subprocess.run(
        (
            ["/bin/zsh", "-lic", "git ls-files --cached --others --exclude-standard docs"]
            if sys.platform == "darwin"
            else ["git", "ls-files", "--cached", "--others", "--exclude-standard", "docs"]
        ),
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    files = [name for name in files if (ROOT / name).exists()]
    errors = []
    for name in files:
        if docs_DATED.search(name):
            errors.append(f"{name}: dated path; move snapshots to reports/")
        elif Path(name).suffix in docs_ARTIFACTS:
            errors.append(f"{name}: run artifact; keep drivers in validation/ and output in temp/")
        elif name.startswith("docs/benchmarks/") and not name.startswith(docs_BENCHMARKS):
            errors.append(f"{name}: not regenerated evidence; see docs/development/README.md")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"docs tree checked ({len(files)} files)")


# The runner executes many named cases in one Stata session, so a case that
# draws random numbers before its own `set seed` gets whatever RNG state the
# previous case left and its data change with the selection and order of cases.
seeds_PREFIX = (
    r"(?:(?:qui(?:e|et|etl|etly)?|n(?:oi|ois|oisi|oisil|oisily)|cap(?:t|tu|tur|ture)?"
    r"|(?:by|bys|byso|bysor|bysort)\s[^:]*:|version\s+[\d.]+\s*:|frame\s+\S+\s*:)\s+)*"
)
seeds_FUNCTION = re.compile(
    r"(?<![\w.])(?:runiform|runiformint|rnormal|rbeta|rbinomial|rcauchy|rchi2|rexponential"
    r"|rgamma|rhypergeometric|rigaussian|rlaplace|rlogistic|rnbinomial|rpoisson|rt|rweibull"
    r"|rweibullph|uniform|rdiscrete|jumble|_jumble)\s*\("
)
seeds_COMMAND = re.compile(
    r"^\s*" + seeds_PREFIX + r"(?:sample|bsample|csample|cbsample|splitsample|drawnorm"
    r"|bootstrap|permute|simulate|corr2data)\b(?!\s*=)"
)
seeds_SEED = re.compile(
    r"^\s*" + seeds_PREFIX + r"(?:set\s+(?:seed|rngstate)\b|(?:drawnorm|bootstrap|permute"
    r"|simulate|corr2data)\b.*\bseed\(\s*[^)\s])|\brseed\(\s*[^)\s]"
)
# The bootstrap options of the ctools estimators draw from Stata's generator
# (docs/development/BOOTSTRAP_PLAN.md, 5.8); a seed() inside them seeds it.
# Quoted text is ignored, so messages that name the options do not count.
seeds_BOOT_COMMAND = r"^\s*" + seeds_PREFIX + r"(?:creghdfe|cppmlhdfe|civreghdfe)\b.*"
seeds_BOOT_OPTION = r"(?:\b(?:vce|fevce)\(\s*boot|\bwildbootstrap\b)"
seeds_BOOT = re.compile(seeds_BOOT_COMMAND + seeds_BOOT_OPTION)
seeds_BOOT_SEED = re.compile(seeds_BOOT_COMMAND + seeds_BOOT_OPTION + r".*\br?seed\(\s*[^)\s]")


def seeds_unquoted(statement):
    return re.sub(r'"[^"]*"', '""', statement)


seeds_CALL = re.compile(r"^\s*" + seeds_PREFIX + r"([A-Za-z_]\w*)")
seeds_DO = re.compile(
    r"^\s*" + seeds_PREFIX + r"(?:do|run)\s+\"?(validation/[\w./]+\.do)\"?(?:\s+\"?(\w+))?"
)
seeds_PROGRAM = re.compile(
    r"^\s*(?:cap(?:t|tu|tur|ture)?\s+)?(?:qui(?:etly)?\s+)?program\s+(?:define\s+|def\s+)?"
    r"(?!drop\b|dir\b|list\b)([A-Za-z_]\w*)(?!.*,\s*plugin\b)"
)
seeds_MATA = re.compile(r"^\s*" + seeds_PREFIX + r"mata\s*:?\s*\{?\s*$")


def seeds_statements(text, first=1):
    """(line, statement) pairs without comments, with /// continuations joined."""
    text = re.sub(r"/\*.*?\*/", lambda m: "\n" * m[0].count("\n"), text, flags=re.S)
    found, pending, start = [], "", None
    for number, line in enumerate(text.split("\n"), first):
        if not pending and re.match(r"\s*\*", line):
            continue
        joined = re.search(r"(?:^|\s)///", line)
        line = re.sub(r"(?:^|\s)//.*$", "", line)
        start = number if start is None else start
        pending += " " + line
        if not joined:
            found.append((start, pending.strip()))
            pending, start = "", None
    return found


def seeds_split(statements):
    """Separate executed statements from `program` bodies (Mata blocks run in place)."""
    top, programs, stack = [], {}, []
    for number, statement in statements:
        if stack and re.fullmatch(r"end", statement):
            stack.pop()
            continue
        if stack and stack[-1] == "mata":
            (programs[stack[-2]][-1] if len(stack) > 1 else top).append((number, statement))
            continue
        program = seeds_PROGRAM.match(statement)
        if program and not stack:
            stack.append(program[1])
            programs.setdefault(program[1], []).append([(number, statement)])
            top.append((number, ("program", program[1])))
            continue
        if seeds_MATA.match(statement):
            stack.append("mata")
            continue
        if stack:
            programs[stack[-1]][-1].append((number, statement))
        else:
            top.append((number, statement))
    return top, {name: [body[1:] for body in bodies] for name, bodies in programs.items()}


def seeds_blocks(text, marker):
    """Yield (name, statements) for each named block, the shared preamble first."""
    starts = list(re.finditer(r'^if "`' + marker + r"'\" == \"([a-z0-9_]+)\" \{\n", text, re.M))
    preamble = seeds_statements(text[: starts[0].start()]) if starts else []
    for i, match in enumerate(starts):
        end = starts[i + 1].start() if i + 1 < len(starts) else len(text)
        first = text[: match.end()].count("\n") + 1
        yield match[1], preamble + seeds_statements(text[match.end() : end], first)


def seeds_findings(root):
    """Return 'path:line [case] statement' for each case that draws before seeding."""
    root = Path(root)
    cases = {}
    for folder, marker in (
        ("validation/stata", "__case"),
        ("validation/performance", "__workload"),
    ):
        for path in sorted((root / folder).glob("*.do")):
            for name, statements in seeds_blocks(path.read_text(errors="replace"), marker):
                cases[(path.relative_to(root).as_posix(), name)] = seeds_split(statements)
    files = {}
    for path in sorted((root / "validation/support").glob("*.do")):
        files[path.relative_to(root).as_posix()] = seeds_split(
            seeds_statements(path.read_text(errors="replace"))
        )
    programs = {}
    for _, defined in [*cases.values(), *files.values()]:
        for name, bodies in defined.items():
            programs.setdefault(name, []).extend(bodies)
    memo = {}

    def first_event(statements, active):
        """('seed'|'draw', line, statement) for the first RNG event, or None."""
        for number, statement in statements:
            if isinstance(statement, tuple):
                continue
            unquoted = seeds_unquoted(statement)
            if seeds_SEED.search(statement) or seeds_BOOT_SEED.search(unquoted):
                return ("seed", number, statement)
            if (
                seeds_FUNCTION.search(statement)
                or seeds_COMMAND.match(statement)
                or seeds_BOOT.search(unquoted)
            ):
                return ("draw", number, statement)
            called = seeds_DO.match(statement)
            if called:
                event = target_event(("do",) + called.groups(), active)
            else:
                called = seeds_CALL.match(statement)
                event = target_event(("program", called[1]), active) if called else None
            if event:
                return (event[0], number, f"{statement}  [{event[0]} at line {event[1]}]")
        return None

    def target_event(key, active):
        if key in active:
            return None
        if key not in memo:
            if key[0] == "program":
                bodies = programs.get(key[1], [])
            elif key[2] and (key[1], key[2]) in cases:
                bodies = [cases[(key[1], key[2])][0]]
            else:
                bodies = [files[key[1]][0]] if key[1] in files else []
            # A name defined twice is judged by its worst definition.
            events = [first_event(body, active | {key}) for body in bodies]
            memo[key] = next((e for e in events if e and e[0] == "draw"), None) or next(
                (e for e in events if e), None
            )
        return memo[key]

    findings = []
    for (path, name), (top, defined) in cases.items():
        seeded = False
        for number, statement in top:
            if isinstance(statement, tuple):
                # A program defined here may be called under any name, so its
                # body must not draw unless the case has already been seeded.
                events = [first_event(body, frozenset()) for body in defined[statement[1]]]
                event = next((e for e in events if e and e[0] == "draw"), None)
            else:
                event = first_event([(number, statement)], frozenset())
            if event and event[0] == "seed":
                seeded = True
            elif event and not seeded:
                findings.append(f"{path}:{event[1]} [{name}] {event[2][:200]}")
                break
    return findings


def seeds_selftest():
    with tempfile.TemporaryDirectory(prefix="ctools-seeds-") as tmp:
        root = Path(tmp)
        for folder in ("stata", "performance", "support"):
            (root / "validation" / folder).mkdir(parents=True)
        (root / "validation/support/helpers.do").write_text(
            "program make_data\n    gen x = rnormal()\nend\n"
            "program seeded_data\n    args s\n    set seed `s'\n    gen x = runiform()\nend\n"
        )
        (root / "validation/support/setup.do").write_text("program quiet_helper\nend\n")
        case = 'if "`__case\'" == "{}" {{\n{}\n}}\n'
        clean = {
            "seeded": "do validation/support/setup.do\nset seed 1\ngen x = runiform()",
            "helper_seeds": "seeded_data 4\nmake_data",
            "comments": "* gen x = runiform()\n// csample 10\n/* rnormal() */\n"
            "set seed 2\ncsample 5",
            "mata_seed": "mata: rseed(3)\ngen x = rnormal()",
            "drawnorm": "drawnorm x, n(10) seed(5)\ngen y = runiform()",
            "subcase": "do validation/stata/other.do seeded\ngen y = runiform()",
            "program_after_seed": "set seed 6\nprogram local_draw\n    gen x = runiform()\nend",
            "strings": 'di "runiform values"\nforeach sample in a b {\n    di "`sample\'"\n}',
            "boot_seeded": "creghdfe y x, absorb(g) vce(bootstrap, reps(5) seed(3))\ngen x = rnormal()",
            "boot_message": 'di "vce(bootstrap) is new"\nset seed 1\ngen x = runiform()',
        }
        dirty = {
            "unseeded": "gen x = runiform()",
            "late_seed": "gen x = rpoisson(2)\nset seed 1",
            "continued": "gen x = 1 + ///\n    rnormal()\nset seed 1",
            "helper_draws": "make_data\nset seed 1",
            "sampling": "quietly capture noisily bysort g: csample 10",
            "mata_block": "mata:\nx = runiform(2, 2)\nend",
            "program_before_seed": "program p_draw\n    bsample\nend\nset seed 1\np_`which'",
            "unseeded_subcase": "do validation/stata/other.do unseeded",
            "plugin_declaration": 'program p_io, plugin using("x.plugin")\ngen x = runiform()',
            "boot_unseeded": "creghdfe y x, absorb(g) vce(bootstrap, reps(5))\nset seed 1",
            "wild_unseeded": "creghdfe y x, absorb(g) vce(cluster g) wildbootstrap(reps(99))",
        }
        blocks = {**clean, **dirty}
        (root / "validation/stata/other.do").write_text(
            "gettoken __case 0 : 0\n"
            + "".join(case.format(name, body) for name, body in blocks.items())
        )
        (root / "validation/performance/work.do").write_text(
            'if "`__workload\'" == "bench" {\n    set obs 10\n    gen x = rnormal()\n}\n'
        )
        flagged = {re.search(r"\[(\w+)\]", line)[1] for line in seeds_findings(root)}
        assert flagged == set(dirty) | {"bench"}, sorted(flagged ^ (set(dirty) | {"bench"}))
    print(f"PASS seed checker flags {len(dirty) + 1} unseeded and passes {len(clean)} seeded cases")


def seeds_main(argv=None):
    seeds_selftest()
    findings = seeds_findings(ROOT)
    if findings:
        raise SystemExit("Cases draw random numbers before `set seed`:\n" + "\n".join(findings))
    print("PASS every named validation case and workload seeds before its first random draw")


def windows_main(argv=None):
    if os.name != "nt":
        raise SystemExit("This regression requires a real Windows worker")
    with tempfile.TemporaryDirectory(prefix="ctools-win-paths-") as tmp:
        binary = Path(tmp) / "paths.exe"
        subprocess.run(
            [
                os.environ.get("CC", "gcc"),
                "-std=c11",
                "-O1",
                "-DSYSTEM=STWIN32",
                "-DSD_FASTMODE",
                "-I",
                str(ROOT / "src"),
                str(unit(Path(tmp), "oct21_windows")),
                str(ROOT / "src/formats/ctools_output.c"),
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run([str(binary), tmp], check=True)
        build, source = Path(tmp) / "build", Path(tmp) / "src"
        build.mkdir()
        source.mkdir()
        for name in ("ctools.plugin", "ctools_windows.plugin", "stale.o"):
            (build / name).touch()
        (source / "stale.o").touch()
        subprocess.run(
            [
                "make",
                "clean",
                "DETECTED_OS=Windows",
                f"BUILD_DIR={build.as_posix()}",
                f"SRC_DIR={source.as_posix()}",
            ],
            cwd=ROOT,
            check=True,
        )
        assert not list(build.iterdir()) and not list(source.iterdir())


def run(out):
    """Run independently so deliberate mocks cannot contaminate later checks."""
    import subprocess
    import sys

    targets = [
        ("build.compiler_probes", "validation.check_build", "build_compiler_probes"),
        ("build.dependencies", "validation.check_build", "build_dependency_metadata"),
        ("build.fail_fast", "validation.check_build", "build_fail_fast_builds"),
        ("build.incremental", "validation.check_build", "build_incremental_builds"),
        ("build.depfile_paths", "validation.check_build", "build_depfile_paths"),
        ("release.archive", "validation.check_build", "archive_main"),
        ("release.workflow", "validation.checks", "release_workflow_contract"),
        ("release.macos_signing", "validation.check_macos_signing", "main"),
        ("release.integrity", "validation.checks", "release_packages"),
        ("release.cleanup", "validation.checks", "release_cleanup"),
        ("package.staging", "validation.checks", "package_packages"),
        ("package.inventory", "validation.checks", "package_inventory"),
        ("package.stata", "validation.checks", "package_stata_contract"),
        ("package.description", "validation.checks", "package_description"),
        ("package.licenses", "validation.checks", "package_licenses"),
        ("runner", "validation.check_runner", "test_runner"),
        ("options", "validation.check_runner", "options_main"),
        ("benchmark.evidence", "validation.check_runner", "test_performance"),
        ("docs", "validation.checks", "docs_main"),
        ("seeds", "validation.checks", "seeds_main"),
    ]
    failures = []
    for name, module, function in targets:
        log = out / (name + ".log")
        with log.open("w") as handle:
            result = subprocess.run(
                [sys.executable, "-c", f"from {module} import {function}; {function}()"],
                cwd=ROOT,
                stdout=handle,
                stderr=subprocess.STDOUT,
            )
        if result.returncode:
            failures.append(name)
            print(f"FAIL {name}: {log}\n" + log.read_text(errors="replace")[-6000:], flush=True)
        else:
            print("PASS " + name, flush=True)
    return failures


def main(argv):
    import argparse

    from validation import check_build

    if argv[0] == "dependencies":
        return check_build.dependencies_main(argv[1:])
    if argv[0] == "windows":
        return windows_main(argv[1:])
    if argv[0] == "docs":
        return docs_main(argv[1:])
    if argv[0] == "seeds":
        return seeds_main(argv[1:])
    if argv[0] == "options":
        from validation import options

        parser = argparse.ArgumentParser(description="Check the public option inventory")
        parser.add_argument("--write", action="store_true")
        args = parser.parse_args(argv[1:])
        rows = options.records()
        gaps = [f"{r['command']}.{r['option']}" for r in rows if r["status"] == "gap"]
        if gaps:
            raise SystemExit("Options without coverage: " + ", ".join(gaps))
        document = ROOT / "docs/OPTION_COVERAGE.md"
        expected = options.report(rows)
        if args.write:
            document.write_text(expected)
        elif document.read_text() != expected:
            raise SystemExit("Option inventory is stale; run --checks options --write")
        print(f"{len(rows)} public options checked; no coverage gaps")
        return
    if argv[0] == "package":
        from scripts.package_contract import check_package

        parser = argparse.ArgumentParser()
        parser.add_argument("directory", type=Path)
        parser.add_argument("--metadata-only", action="store_true")
        parser.add_argument("--integrity", action="store_true")
        args = parser.parse_args(argv[1:])
        errors = check_package(args.directory, args.metadata_only, args.integrity)
        if errors:
            raise SystemExit("\n".join(errors))
        print("Package checked")
        return
    raise SystemExit("Unknown repository check: " + argv[0])
