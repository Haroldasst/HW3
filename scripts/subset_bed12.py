import argparse
import random
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--n", type=int, default=10000)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

if args.n <= 0:
    raise ValueError("n turi būti teigiamas skaičius")

rng = random.Random(args.seed)
reservoir = []
total = 0

with open(args.input) as infile:
    for line in infile:
        if not line.strip() or line.startswith("#"):
            continue

        total += 1

        if len(reservoir) < args.n:
            reservoir.append((total, line))
        else:
            j = rng.randrange(total)
            if j < args.n:
                reservoir[j] = (total, line)

if total == 0:
    raise ValueError("Įvesties BED12 failas tuščias")

selected = sorted(reservoir, key=lambda x: x[0])

output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w") as outfile:
    for _, line in selected:
        outfile.write(line)

print(f"Viso transkriptų modelių: {total}")
print(f"Atrinkta transkriptų modelių: {len(selected)}")
print(f"Atsitiktinės atrankos seed: {args.seed}")
