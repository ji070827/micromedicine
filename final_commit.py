import subprocess
import py_compile
import glob

log = open("final_commit.log", "w", encoding="utf-8", buffering=1)

def w(s):
    log.write(s + "\n")
    log.flush()

# 1. 语法验证
w("=== 语法验证 ===")
files = ['app.py', 'scripts/build_molecule_report.py', 'scripts/generate_large_library.py']
for f in files:
    try:
        py_compile.compile(f, doraise=True)
        w(f"OK   {f}")
    except Exception as e:
        w(f"FAIL {f}: {e}")

# 2. git 提交推送
w("\n=== git 提交 ===")
for cmd in [["git","add","-A"], ["git","commit","-m","前端改造：逐分子成药性详情+2D结构SVG+实时进度轮询+超算结果文件检测"], ["git","push","origin","main"]]:
    w(f"\n>>> {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=120)
    w(f"rc={r.returncode}")
    if r.stdout: w(r.stdout[-1000:])
    if r.stderr: w("STDERR: " + r.stderr[-800:])

log.close()
import os
os.remove("final_commit.py") if os.path.exists("final_commit.py") else None
print("done")