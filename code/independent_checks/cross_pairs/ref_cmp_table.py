#!/usr/bin/env python3
# Compare the manuscript's Section 8 KR table (theorem_yy.txt) with the referee's computation (summary file).
import re, sys
ms = open(sys.argv[1], encoding="utf-8").read()
sec = ms.split("8. TABLE")[1].split("Observations")[0]
man = {int(a): (int(b), int(c)) for a, b, c in re.findall(r"(\d+):(\d+)\((\d+)\)", sec)}
mine_txt = open(sys.argv[2]).read()
mine = {int(a): (int(b), int(c)) for a, b, c in re.findall(r"(\d+):(\d+)\((\d+)\)", mine_txt)}
diff = [(d, man[d], mine.get(d)) for d in sorted(man) if mine.get(d) != man[d]]
print("manuscript entries:", len(man), "range", min(man), max(man), " mismatches:", diff)
