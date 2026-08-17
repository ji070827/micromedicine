import subprocess

log = open("commit_whu.log", "w", encoding="utf-8", buffering=1)

def run(cmd):
    log.write(f"\n>>> {' '.join(cmd)}\n")
    log.flush()
    r = subprocess.run(cmd, capture_output=True, timeout=120,
                       encoding='utf-8', errors='replace')
    log.write(f"--- rc={r.returncode}\n")
    if r.stdout:
        log.write(r.stdout + "\n")
    if r.stderr:
        log.write("STDERR:\n" + r.stderr + "\n")
    log.flush()
    return r.returncode

run(["git", "add", "-A"])
run(["git", "status", "-sb"])
run(["git", "commit", "-m", "新增武大超算部署方案：setup_whu.sh + slurm作业脚本 + AF3支持Singularity + WHU部署文档"])
run(["git", "push", "origin", "main"])
run(["git", "log", "--oneline", "-3"])

log.close()
import os
os.remove("commit_whu.py")
print("done")