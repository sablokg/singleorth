#!/usr/bin/env bash

# Gaurav Sablok
# gsablok@proton.me
#
# Pangenome orthogroup alignment + phylogeny pipeline.
# Supports macse, prank, muscle, or all three.

set -euo pipefail

echo "Gaurav Sablok gsablok@proton.me"
echo "Pangenome orthology alignment + phylogeny pipeline"
echo

# ---------------------------------------------------------------------------
# Input collection
# ---------------------------------------------------------------------------

read -r -p "Enter the directory path to the sequences: " dirpath
if [[ ! -d "${dirpath}" ]]; then
	echo "Error: '${dirpath}' is not a directory." >&2
	exit 1
fi
directorypath="$(cd "${dirpath}" && pwd)"

if ! compgen -G "${directorypath}"/*.fa >/dev/null; then
	echo "Error: no .fa files found in ${directorypath}" >&2
	exit 1
fi

read -r -p "Enter the number of species involved in the orthology search: " species
if ! [[ "${species}" =~ ^[0-9]+$ ]]; then
	echo "Error: species count must be a positive integer." >&2
	exit 1
fi

read -r -p "Enter the number of threads: " threads
if ! [[ "${threads}" =~ ^[0-9]+$ ]]; then
	echo "Error: thread count must be a positive integer." >&2
	exit 1
fi

# Requires the jar path only when macse is (or might be) used; asked
# unconditionally is simplest and harmless if unused.
read -r -p "Enter the path to macse.jar (leave blank if not using macse): " macse_jar

echo
echo "Select the alignment tool:"
PS3="Enter choice [1-4]: "
options=("macse" "prank" "muscle" "all")
select tool in "${options[@]}"; do
	case "${tool}" in
	macse | prank | muscle | all)
		break
		;;
	*)
		echo "Invalid selection, choose 1-4."
		;;
	esac
done
echo "Selected tool: ${tool}"
echo

# ---------------------------------------------------------------------------
# Validate that every .fa file has the expected number of sequences
# ---------------------------------------------------------------------------

echo "Checking that all files contain ${species} sequences..."
bad_files=()
for f in "${directorypath}"/*.fa; do
	count=$(grep -c ">" "${f}")
	if [[ "${count}" != "${species}" ]]; then
		bad_files+=("${f} (found ${count})")
	fi
done

if ((${#bad_files[@]} > 0)); then
	echo "Error: the following files do not have ${species} sequences:" >&2
	printf '  %s\n' "${bad_files[@]}" >&2
	exit 1
fi
echo "All files verified: each contains ${species} sequences."
echo

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

AMAS_URL="https://raw.githubusercontent.com/marekborowiec/AMAS/master/amas/AMAS.py"

fetch_amas() {
	if [[ ! -f AMAS.py ]]; then
		wget -q "${AMAS_URL}" -O AMAS.py
		chmod 755 AMAS.py
	fi
}

# concatenate + partition + tree-build, run inside the tool-specific dir
run_downstream() {
	local prefix="$1" # e.g. macse, prank, muscle
	fetch_amas
	python3 AMAS.py -in ./*.trimmed.fasta -f fasta -d dna -c "${threads}" \
		-p "${prefix}alignmentpartitions.txt" -t "${prefix}alignmentconcatenated.fasta"

	iqtree --seqtype DNA -s "${prefix}alignmentconcatenated.fasta" \
		--alrt 1000 -b 1000 -T "${threads}"

	raxmlHPC-PTHREADS -s "${prefix}alignmentconcatenated.fasta" --no-seq-check -O \
		-m GTRGAMMA -p 12345 -n "${prefix}phylogeny_GAMMA" -T "${threads}" -N 50

	raxmlHPC-PTHREADS -s "${prefix}alignmentconcatenated.fasta" --no-seq-check -O \
		-m GTRGAMMA -p 12345 -n "${prefix}phylogeny_GTRCAT" -T "${threads}" -N 50 -b 1000
}

run_macse() {
	if [[ -z "${macse_jar}" || ! -f "${macse_jar}" ]]; then
		echo "Error: macse selected but no valid --macse jar path was given." >&2
		exit 1
	fi
	echo "== Aligning with MACSE =="
	mkdir -p "${directorypath}/macse_run"
	cd "${directorypath}/macse_run"

	for j in "${directorypath}"/*.fa; do
		base="$(basename "${j%.*}")"
		java -Xmx100g -jar "${macse_jar}" -prog alignSequences \
			-gc_def 12 -seq "${j}" \
			-out_AA "${base}.AA" -out_NT "${base}.NT" \
			>"${base}.macse.run.log.txt"
	done

	for i in *.NT; do
		mv "${i}" "${i%.*}.ntaligned.fasta"
	done
	echo "Renamed aligned files with .ntaligned.fasta suffix"

	for i in *.ntaligned.fasta; do
		trimal -in "${i}" -out "${i%.*}.trimmed.fasta" -nogaps
	done

	run_downstream "macse"
	echo "MACSE analysis complete. Output in ${directorypath}/macse_run"
}

run_prank() {
	echo "== Aligning with PRANK =="
	command -v prank >/dev/null 2>&1 || sudo apt-get install -y prank
	mkdir -p "${directorypath}/prank_run"
	cd "${directorypath}/prank_run"

	for i in "${directorypath}"/*.fa; do
		base="$(basename "${i%.*}")"
		prank -d="${i}" -o="${base}.prankaligned" -codon
	done

	for i in *.prankaligned.best.fas; do
		[[ -e "${i}" ]] || continue
		trimal -in "${i}" -out "${i%.*}.trimmed.fasta" -nogaps
	done

	run_downstream "prank"
	echo "PRANK analysis complete. Output in ${directorypath}/prank_run"
}

run_muscle() {
	echo "== Aligning with MUSCLE =="
	command -v muscle >/dev/null 2>&1 || sudo apt-get install -y muscle
	mkdir -p "${directorypath}/muscle_run"
	cd "${directorypath}/muscle_run"

	for i in "${directorypath}"/*.fa; do
		base="$(basename "${i%.*}")"
		muscle -in "${i}" -out "${base}.musclealigned.fasta"
	done

	for i in *.musclealigned.fasta; do
		trimal -in "${i}" -out "${i%.*}.trimmed.fasta" -nogaps
	done

	run_downstream "muscle"
	echo "MUSCLE analysis complete. Output in ${directorypath}/muscle_run"
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

case "${tool}" in
macse) run_macse ;;
prank) run_prank ;;
muscle) run_muscle ;;
all)
	run_macse
	run_prank
	run_muscle
	;;
esac

echo
echo "Pipeline finished."
