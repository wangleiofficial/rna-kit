from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from Bio.Align import PairwiseAligner

from .exceptions import MetricCalculationError
from .structures import PDBStructure, ResidueRecord


@dataclass(frozen=True)
class ChainAlignment:
    native_chain: str
    prediction_chain: str
    matched_residues: int
    native_indices: tuple[int, ...]
    prediction_indices: tuple[int, ...]
    native_sequence: str
    prediction_sequence: str


@dataclass(frozen=True)
class StructureAlignment:
    native_indices: tuple[int, ...]
    prediction_indices: tuple[int, ...]
    native_index_spec: str
    prediction_index_spec: str
    matched_residues: int
    chain_alignments: tuple[ChainAlignment, ...]


class StructureMatcher:
    def __init__(self) -> None:
        self.aligner = PairwiseAligner(mode="global")
        self.aligner.match_score = 2.0
        self.aligner.mismatch_score = -1.0
        self.aligner.open_gap_score = -2.0
        self.aligner.extend_gap_score = -0.5

    def align(self, native: PDBStructure, prediction: PDBStructure) -> StructureAlignment:
        native_chains = native.chain_records()
        prediction_chains = prediction.chain_records()
        if not native_chains or not prediction_chains:
            raise MetricCalculationError("At least one residue is required in both structures.")

        native_chain_names = list(native_chains)
        prediction_chain_names = list(prediction_chains)
        candidates = {
            (native_chain, prediction_chain): self._align_chain(
                native_chain,
                native_chains[native_chain],
                prediction_chain,
                prediction_chains[prediction_chain],
            )
            for native_chain in native_chain_names
            for prediction_chain in prediction_chain_names
        }

        chosen_pairs = self._choose_chain_pairs(native_chain_names, prediction_chain_names, candidates)
        chain_alignments = tuple(
            candidates[(native_chain, prediction_chain)]
            for native_chain, prediction_chain in chosen_pairs
            if candidates[(native_chain, prediction_chain)].matched_residues > 0
        )
        matched_native_indices = tuple(
            index for alignment in chain_alignments for index in alignment.native_indices
        )
        matched_prediction_indices = tuple(
            index for alignment in chain_alignments for index in alignment.prediction_indices
        )

        if not matched_native_indices or not matched_prediction_indices:
            raise MetricCalculationError(
                "Unable to infer a common residue mapping between the reference and prediction structures."
            )

        native_view = native.with_selected_indices(list(matched_native_indices))
        prediction_view = prediction.with_selected_indices(list(matched_prediction_indices))
        return StructureAlignment(
            native_indices=matched_native_indices,
            prediction_indices=matched_prediction_indices,
            native_index_spec=native_view.index_spec(),
            prediction_index_spec=prediction_view.index_spec(),
            matched_residues=len(matched_native_indices),
            chain_alignments=chain_alignments,
        )

    def _choose_chain_pairs(
        self,
        native_chain_names: list[str],
        prediction_chain_names: list[str],
        candidates: dict[tuple[str, str], ChainAlignment],
    ) -> tuple[tuple[str, str], ...]:
        @lru_cache(maxsize=None)
        def solve(native_index: int, used_prediction_mask: int) -> tuple[int, tuple[tuple[str, str], ...]]:
            if native_index >= len(native_chain_names):
                return 0, ()

            best_score, best_pairs = solve(native_index + 1, used_prediction_mask)
            native_chain = native_chain_names[native_index]

            for prediction_index, prediction_chain in enumerate(prediction_chain_names):
                if used_prediction_mask & (1 << prediction_index):
                    continue

                candidate = candidates[(native_chain, prediction_chain)]
                if candidate.matched_residues == 0:
                    continue

                next_score, next_pairs = solve(
                    native_index + 1,
                    used_prediction_mask | (1 << prediction_index),
                )
                total_score = candidate.matched_residues + next_score
                if total_score > best_score:
                    best_score = total_score
                    best_pairs = ((native_chain, prediction_chain),) + next_pairs

            return best_score, best_pairs

        return solve(0, 0)[1]

    def _align_chain(
        self,
        native_chain: str,
        native_records: list[tuple[int, ResidueRecord]],
        prediction_chain: str,
        prediction_records: list[tuple[int, ResidueRecord]],
    ) -> ChainAlignment:
        native_sequence = "".join(record.nt for _, record in native_records)
        prediction_sequence = "".join(record.nt for _, record in prediction_records)
        alignment = self.aligner.align(native_sequence, prediction_sequence)[0]

        matched_native_indices: list[int] = []
        matched_prediction_indices: list[int] = []
        for (native_start, native_end), (prediction_start, prediction_end) in zip(
            alignment.aligned[0], alignment.aligned[1]
        ):
            block_length = min(native_end - native_start, prediction_end - prediction_start)
            for offset in range(block_length):
                native_position = native_start + offset
                prediction_position = prediction_start + offset
                if native_sequence[native_position] != prediction_sequence[prediction_position]:
                    continue
                matched_native_indices.append(native_records[native_position][0])
                matched_prediction_indices.append(prediction_records[prediction_position][0])

        return ChainAlignment(
            native_chain=native_chain,
            prediction_chain=prediction_chain,
            matched_residues=len(matched_native_indices),
            native_indices=tuple(matched_native_indices),
            prediction_indices=tuple(matched_prediction_indices),
            native_sequence=native_sequence,
            prediction_sequence=prediction_sequence,
        )


def infer_structure_alignment(
    native: PDBStructure,
    prediction: PDBStructure,
) -> StructureAlignment:
    return StructureMatcher().align(native, prediction)
