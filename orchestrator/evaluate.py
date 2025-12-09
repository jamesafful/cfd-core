
import json, subprocess, time, pathlib, hashlib, statistics, os, math
from schemas import Metrics

ROOT = pathlib.Path(__file__).resolve().parents[1]

def build_dir_for(flags: str) -> pathlib.Path:
    h = hashlib.sha1(flags.encode()).hexdigest()[:10]
    p = ROOT / "build" / h
    p.mkdir(parents=True, exist_ok=True)
    return p

def run_build(flags: str, bd: pathlib.Path):
    cfg = f'cmake -S {ROOT} -B {bd} -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="{flags}"'
    subprocess.check_call(cfg, shell=True)
    subprocess.check_call(f"cmake --build {bd} -j", shell=True)

def run_case(bd: pathlib.Path, nx: int = 800, cfl: float = 0.5, rk: int = 3, deterministic: bool=True):
    exe = bd / "agentic_cfd"
    t0 = time.time()
    out = subprocess.check_output(
        f"{exe} --problem sod1d --nx {nx} --final_time 0.2 --cfl {cfl} --rk {rk} {'--deterministic' if deterministic else ''}",
        shell=True, text=True)
    dt = (time.time()-t0)*1000.0
    m = json.loads(out)
    runtime = float(m.get("kernel_ms", dt))
    return runtime, m

_BASELINE_FILE = ROOT / "build" / "baseline.json"

def get_baseline_ms() -> float:
    if _BASELINE_FILE.exists():
        return json.loads(_BASELINE_FILE.read_text())['baseline_ms']
    bd = build_dir_for("")
    run_build("", bd)
    times = []
    for _ in range(3):
        t, _ = run_case(bd, nx=800, cfl=0.5, rk=3, deterministic=True)
        times.append(t)
    base = statistics.mean(times)
    _BASELINE_FILE.parent.mkdir(parents=True, exist_ok=True)
    _BASELINE_FILE.write_text(json.dumps({"baseline_ms": base}))
    return base

def evaluate(cfg: dict) -> Metrics:
    bd = build_dir_for(cfg.get("compiler_flags",""))
    run_build(cfg.get("compiler_flags",""), bd)
    kruns = int(os.environ.get("KRUNS", 3))
    times=[]; last={}
    for _ in range(kruns):
        t, m = run_case(bd, nx=800, cfl=cfg["cfl"], rk=cfg["rk_stages"], deterministic=True)
        times.append(t); last=m
    runtime = statistics.mean(times)
    l2 = float(last.get("l2", 1.0))
    speedup = get_baseline_ms() / max(runtime, 1e-6)
    score = ( -math.log(max(l2,1e-12)) * 0.4 + speedup * 0.4 + 0.2 )
    return Metrics(l2=l2, runtime_ms=runtime, determinism=1.0, score=score)
