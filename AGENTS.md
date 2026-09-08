# AGENTS.md — Build Concurrency Policy (RAM-gated)

## Rule: choose `-j` from available RAM before compiling

1. Check available RAM first:
   - `free -m | awk '/^Mem:/ {print $7}'` (MB available), or
   - `awk '/MemAvailable/ {print int($2/1024)}' /proc/meminfo` (MB available)
2. Mapping:
   - `< 1000 MB` available → DO NOT compile. Wait: `sleep 60`, then re-check in a loop until `>= 1000 MB`.
   - `1000–1999 MB` available → `-j1`
   - `2000–2999 MB` available → `-j2`
   - `>= 3000 MB` available → `-j3` (max; never exceed `-j3`)
3. Example wait loop:
   ```bash
   while true; do
     AVAIL=$(awk '/MemAvailable/ {print int($2/1024)}' /proc/meminfo)
     echo "Available RAM: ${AVAIL} MB"
     if [ "$AVAIL" -ge 1000 ]; then break; fi
     echo "Less than 1GB available, waiting 60s..."
     sleep 60
   done
   if [ "$AVAIL" -ge 3000 ]; then JOBS=3;
   elif [ "$AVAIL" -ge 2000 ]; then JOBS=2;
   else JOBS=1; fi
   make -j"$JOBS"
   ```
4. Re-check RAM right before each compile/link step; if it dropped below 1GB, stop and go back to waiting.
