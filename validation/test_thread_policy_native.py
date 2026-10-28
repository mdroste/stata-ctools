"""Thread policy regressions: default cap, precedence, CPU limits, saved setting.

The harness includes ctools_threads.c. Each case runs in a fresh process with a
controlled environment and a temporary home, so the developer's own saved
setting and OMP_* variables never leak into the results. On Linux the affinity
case also runs under taskset. CC may name a multi-word compiler command.
"""
from native.support.compiler import compile_native
from native.support.runner import harness
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SYSTEM = platform.system()
CAP = 16
SCRUB = ("CTOOLS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "XDG_CONFIG_HOME")


def build(tmp):
    source = tmp / "thread_policy.c"
    source.write_text(harness("thread_policy.c"))
    binary = tmp / ("thread_policy.exe" if SYSTEM == "Windows" else "thread_policy")
    default = "/usr/bin/clang" if SYSTEM == "Darwin" else "cc"
    command = [*shlex.split(os.environ.get("CC", default)), "-std=gnu11", "-O1", "-g",
               "-I", str(ROOT / "src")]
    sanitizers = os.environ.get("CTOOLS_SANITIZERS", "undefined")
    if sanitizers:
        command += ["-fno-omit-frame-pointer", "-fsanitize=" + sanitizers]
    command += {"Darwin": ["-DSYSTEM=APPLEMAC"], "Windows": ["-DSYSTEM=STWIN32"]}.get(
        SYSTEM, ["-DSYSTEM=STUNIX"])
    command += [str(source), "-o", str(binary)]
    if SYSTEM != "Windows":
        command += ["-pthread", "-lm"]
    compile_native(command, check=True, capture_output=True, text=True)
    return binary


def environment(home, **values):
    env = {k: v for k, v in os.environ.items() if k not in SCRUB}
    env.update(HOME=str(home), APPDATA=str(home / "AppData"),
               ASAN_OPTIONS="detect_leaks=0:abort_on_error=1", UBSAN_OPTIONS="halt_on_error=1")
    env.update(values)
    return env


def settings_path(home, env):
    if SYSTEM == "Windows":
        return str(home / "AppData" / "ctools" / "ctools.conf")
    if SYSTEM == "Darwin":
        return f"{home}/Library/Application Support/ctools/ctools.conf"
    base = env.get("XDG_CONFIG_HOME") or f"{home}/.config"
    return f"{base}/ctools/ctools.conf"


def report(binary, env, prefix=()):
    out = subprocess.run([*prefix, str(binary), "report"], env=env, check=True,
                         capture_output=True, text=True, timeout=30).stdout
    # Values may contain spaces ("saved setting"); each ends at the next key.
    fields = dict(re.findall(r"(\w+)=(.*?)(?=\s+\w+=|\s*$)", out))
    return {k: v if k == "source" else int(v) for k, v in fields.items()}


def check(name, got, **expected):
    wrong = {k: (got[k], v) for k, v in expected.items() if got[k] != v}
    assert not wrong, f"{name}: {wrong} in {got}"
    assert got["max"] == got["threads"] and got["clamp"] == got["limit"], (name, got)
    print("PASS", name, flush=True)


def main():
    with tempfile.TemporaryDirectory(prefix="ctools-thread-policy-") as directory:
        tmp = Path(directory)
        binary = build(tmp)

        homes = [tmp / "home"]
        if SYSTEM == "Linux":
            homes.append(tmp / "xdg-home")
        for home in homes:
            home.mkdir()
            (home / "AppData").mkdir()
            env = environment(home)
            if home.name == "xdg-home":
                env["XDG_CONFIG_HOME"] = str(home / "xdg")
            work = home / "work"
            work.mkdir()
            subprocess.run([str(binary), "unit", str(work), settings_path(home, env)], env=env,
                           check=True, capture_output=True, text=True, timeout=60)
            print("PASS unit", home.name, flush=True)

        home = tmp / "cases"
        home.mkdir()
        (home / "AppData").mkdir()
        base = report(binary, environment(home))
        usable = base["available"]
        limit = base["limit"]
        assert 1 <= usable <= base["detected"] and limit == min(usable, 256), base
        check("automatic", base, threads=min(usable, CAP), source="automatic", saved=0)
        check("omp list", report(binary, environment(home, OMP_NUM_THREADS="3,2")),
              threads=min(3, limit), source="OMP_NUM_THREADS", omp=3)
        check("ctools over omp",
              report(binary, environment(home, OMP_NUM_THREADS="3", CTOOLS_NUM_THREADS="2")),
              threads=min(2, limit), source="CTOOLS_NUM_THREADS")
        check("explicit above usable", report(binary, environment(home, CTOOLS_NUM_THREADS="100000")),
              threads=limit, source="CTOOLS_NUM_THREADS")
        check("thread limit", report(binary, environment(home, OMP_THREAD_LIMIT="1",
                                                         CTOOLS_NUM_THREADS="8")),
              threads=1, limit=1, omp_limit=1)
        check("invalid ignored", report(binary, environment(home, CTOOLS_NUM_THREADS="abc",
                                                            OMP_NUM_THREADS="0")),
              threads=min(usable, CAP), source="automatic", env=-1, omp=-1)

        conf = Path(settings_path(home, environment(home)))
        conf.parent.mkdir(parents=True)
        conf.write_text("threads=5\n")
        check("saved over omp", report(binary, environment(home, OMP_NUM_THREADS="3")),
              threads=min(5, limit), source="saved setting", saved=5)
        check("ctools over saved", report(binary, environment(home, CTOOLS_NUM_THREADS="2")),
              threads=min(2, limit), source="CTOOLS_NUM_THREADS", saved=5)
        conf.write_text("threads=lots\n")
        check("invalid saved ignored", report(binary, environment(home)),
              threads=min(usable, CAP), source="automatic", saved=-1)
        conf.unlink()

        taskset = shutil.which("taskset") if SYSTEM == "Linux" else None
        if taskset:
            check("affinity", report(binary, environment(home), [taskset, "-c", "0"]),
                  available=1, limit=1, threads=1)
            check("affinity explicit", report(binary, environment(home, CTOOLS_NUM_THREADS="8"),
                                              [taskset, "-c", "0"]),
                  available=1, threads=1, source="CTOOLS_NUM_THREADS")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        print(exc.stdout or "", exc.stderr or "", flush=True)
        raise
