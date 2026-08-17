import subprocess
import os

log = open("commit_fda.log", "w", encoding="utf-8", buffering=1)

def run(cmd, timeout=120):
    log.write(f"\n>>> {' '.join(cmd)}\n")
    log.flush()
    r = subprocess.run(cmd, capture_output=True, timeout=timeout,
                       encoding='utf-8', errors='replace')
    log.write(f"--- rc={r.returncode}\n")
    if r.stdout:
        log.write(r.stdout + "\n")
    if r.stderr:
        log.write("STDERR:\n" + r.stderr + "\n")
    log.flush()
    return r.returncode

run(["git", "add", "-A"])
run(["git", "commit", "-m", "路线一改用真实FDA批准成药库（47个真实药物），新增通用库加载器/generate_fda_drug_library；同步文档"])
run(["git", "push", "origin", "main"])
run(["git", "log", "--oneline", "-3"])

log.close()

# 清理自身
for f in ["commit_fda.py"]:
    if os.path.exists(f):
        os.remove(f)
print("done")