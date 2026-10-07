from pathlib import Path

import pandas as pd
import pytest

from sma_xspecies.splicing import (
    significant_events,
    standardize_rmats,
)


def make_explicit_rmats_fixture(
    path: Path,
    reported_differences: list[float],
) -> None:
    pd.DataFrame({
        "ID": [1, 2],
        "GeneID": [
            '"ENSG1.2"',
            '"ENSG2.1"',
        ],
        "geneSymbol": [
            "A",
            "B",
        ],
        "chr": [
            "1",
            "2",
        ],
        "strand": [
            "+",
            "-",
        ],
        "exonStart_0base": [
            100,
            200,
        ],
        "exonEnd": [
            150,
            250,
        ],
        "upstreamES": [
            50,
            150,
        ],
        "upstreamEE": [
            80,
            180,
        ],
        "downstreamES": [
            170,
            270,
        ],
        "downstreamEE": [
            200,
            300,
        ],
        "FDR": [
            0.01,
            0.20,
        ],
        "PValue": [
            0.001,
            0.10,
        ],
        # Deposited convention: CTL minus SMA.
        "IncLevelDifference":
            reported_differences,
        "IncLevel_CTL": [
            "0.4,0.4",
            "0.8,0.8",
        ],
        "IncLevel_SMA": [
            "0.2,0.2",
            "0.5,0.5",
        ],
    }).to_csv(
        path,
        sep="\t",
        index=False,
    )


def test_rmats_is_event_level_and_thresholded(
    tmp_path: Path,
):
    path = tmp_path / "SE.MATS.txt"

    make_explicit_rmats_fixture(
        path,
        reported_differences=[
            0.2,
            0.3,
        ],
    )

    events = standardize_rmats(
        path,
        "TEST",
    )

    selected = significant_events(
        events
    )

    assert len(events) == 2
    assert len(selected) == 1

    assert events.loc[
        0,
        "gene_id",
    ] == "ENSG1"

    assert events.loc[
        0,
        "gene_symbol",
    ] == "A"

    assert events.loc[
        0,
        "mean_psi_control",
    ] == pytest.approx(0.4)

    assert events.loc[
        0,
        "mean_psi_disease",
    ] == pytest.approx(0.2)

    # Project convention: SMA minus control.
    assert events.loc[
        0,
        "delta_psi",
    ] == pytest.approx(-0.2)

    assert events.loc[
        0,
        "reported_ctl_minus_sma",
    ] == pytest.approx(0.2)

    assert events.loc[
        0,
        "direction_definition",
    ] == "SMA_minus_control"

    assert selected[
        "gene_symbol"
    ].tolist() == ["A"]

    expected_event_id = (
        "1:+:100:150:50:80:170:200"
    )

    assert events.loc[
        0,
        "event_id",
    ] == expected_event_id

    assert events[
        "event_id"
    ].is_unique

    assert events[
        "event_id"
    ].str.split(":").map(len).eq(8).all()


def test_rmats_rejects_incorrect_orientation(
    tmp_path: Path,
):
    path = tmp_path / "SE_wrong_orientation.MATS.txt"

    # These values incorrectly represent SMA minus CTL
    # while the deposited explicit schema is required
    # to report CTL minus SMA.
    make_explicit_rmats_fixture(
        path,
        reported_differences=[
            -0.2,
            -0.3,
        ],
    )

    with pytest.raises(
        ValueError,
        match=(
            "does not match "
            "CTL-minus-SMA orientation"
        ),
    ):
        standardize_rmats(
            path,
            "TEST",
        )
