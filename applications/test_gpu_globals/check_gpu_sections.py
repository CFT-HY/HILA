#!/usr/bin/env python3
# Check GPU data sections of a CUDA binary: section size vs. extent of the variables in it.
# Usage: python3 check_gpu_sections.py <executable>   (needs cuobjdump; CUOBJDUMP env overrides)
import os, re, subprocess, sys
# For GPU data sections (.nv.constant3, .nv.global, .nv.global.init): section size vs. the
# extent of the variables placed in it.  Uses cuobjdump -elf (CUOBJDUMP env to override).
out = subprocess.run([os.environ.get('CUOBJDUMP', 'cuobjdump'), '-elf', sys.argv[1]],
                     capture_output=True, text=True).stdout
secs = {}
for line in out.splitlines():
    m = re.match(r'\s*([0-9a-f]+)\s+[0-9a-f]+\s+([0-9a-f]+)\s.*\s(\.nv\.(constant3|global|global\.init))$', line)
    if m:
        secs[int(m.group(1), 16)] = [m.group(3), int(m.group(2), 16), 0]
for line in out.splitlines():
    m = re.match(r'\s*0x[0-9a-f]+\s+(0x[0-9a-f]+)\s+(0x[0-9a-f]+)\s+0x11\s+\S+\s+(0x[0-9a-f]+)\s+(\S+)$', line)
    if m and int(m.group(3), 16) in secs:
        s = secs[int(m.group(3), 16)]
        s[2] = max(s[2], int(m.group(1), 16) + int(m.group(2), 16))
bad = False
for name, size, ext in secs.values():
    ok = size >= ext
    bad |= not ok
    print(f"  {name:16s} size {size:#6x}  symbols up to {ext:#6x}  {'OK' if ok else 'BROKEN'}")
sys.exit(1 if bad else 0)
