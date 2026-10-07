#!/usr/bin/env python3

# ============================================================
# PHASE 6 — PRODUCTION GENOME ACQUISITION
#
# Klebsiella Global Genomic Epidemiology
#
# Version: Phase6_Production_v1.1
#
# Workflow:
# Phase5G manifest
#       ↓
# fixed 1,000-genome batches
#       ↓
# local NVMe staging
#       ↓
# NCBI dehydrated package
#       ↓
# rehydrate --gzip
#       ↓
# gzip + FASTA validation
#       ↓
# SHA256 source
#       ↓
# copy to CIFS /mnt/genome
#       ↓
# SHA256 destination verification
#       ↓
# persistent manifest update
#       ↓
# local staging cleanup
#
# Restart-safe.
# ============================================================

from pathlib import Path
from datetime import datetime, timezone

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import time

import pandas as pd


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(
    os.environ.get(
        "PHASE6_ROOT",
        "/mnt/genome/"
        "Klebsiella_Global_Genomic_Epidemiology/"
        "Phase6_Genome_Acquisition",
    )
)

SEED_MANIFEST = (
    ROOT
    / "manifests"
    / "05_Phase6_Download_Manifest_Seed.csv"
)

MASTER_MANIFEST = (
    ROOT
    / "manifests"
    / "Phase6_Production_Master_Manifest.csv"
)

BATCH_TABLE = (
    ROOT
    / "manifests"
    / "Phase6_Batch_Status.csv"
)


GENOME_ROOT = (
    ROOT
    / "genomes"
    / "production"
)

LOG_ROOT = (
    ROOT
    / "logs"
    / "production"
)

FAILED_ROOT = (
    ROOT
    / "failed"
)

QC_ROOT = (
    ROOT
    / "qc"
    / "production"
)

FREEZE_ROOT = (
    ROOT
    / "freeze"
)


# ------------------------------------------------------------
# Local NVMe staging
# ------------------------------------------------------------

LOCAL_ROOT = Path(
    "/var/tmp/"
    "Klebsiella_Phase6_Production"
)

LOCAL_PACKAGE_ROOT = (
    LOCAL_ROOT
    / "packages"
)

LOCAL_REHYDRATE_ROOT = (
    LOCAL_ROOT
    / "rehydrated"
)

LOCAL_ACCESSION_ROOT = (
    LOCAL_ROOT
    / "accessions"
)


# ------------------------------------------------------------
# NCBI datasets executable
# ------------------------------------------------------------

DATASETS_BIN = Path(
    shutil.which("datasets")
    or "/home/romen/bin/datasets"
)


# ------------------------------------------------------------
# Production parameters
# ------------------------------------------------------------

BATCH_SIZE = 1000

MAX_WORKERS = 12

MAX_RETRIES = 3

RETRY_WAIT_SECONDS = 30


# Keep substantial safety margins.
MIN_FINAL_FREE_GB = 300

MIN_LOCAL_FREE_GB = 50


ASSEMBLY_RE = re.compile(
    r"^(?:GCA|GCF)_\d+\.\d+$",
    re.I,
)


PIPELINE_VERSION = (
    "Phase6_Production_v1.1"
)


# ============================================================
# DATA TYPES
# ============================================================

STRING_COLUMNS = [

    "Production_Status",

    "Genome_FASTA_Path",

    "SHA256",

    "MD5",

    "Failure_Reason",

    "Last_Attempt_UTC",

    "Completed_UTC",

    "Acquisition_Source",

    "Pipeline_Version",
]


INTEGER_COLUMNS = [

    "Download_Attempts",

    "File_Size_Bytes",

    "Genome_Length_bp",

    "Contig_Count",
]


FAILURE_STATUSES = {

    "Failed",

    "Missing_From_NCBI",

    "Invalid_Gzip",

    "Invalid_FASTA",

    "Checksum_Failed",
}


# ============================================================
# BASIC HELPERS
# ============================================================

def utc_now():

    return datetime.now(
        timezone.utc
    ).isoformat()


def log(message):

    timestamp = (
        datetime.now()
        .strftime(
            "%Y-%m-%d %H:%M:%S"
        )
    )

    print(
        f"[{timestamp}] {message}",
        flush=True,
    )


def ensure_directories():

    directories = [

        ROOT,

        GENOME_ROOT,

        LOG_ROOT,

        FAILED_ROOT,

        QC_ROOT,

        FREEZE_ROOT,

        LOCAL_ROOT,

        LOCAL_PACKAGE_ROOT,

        LOCAL_REHYDRATE_ROOT,

        LOCAL_ACCESSION_ROOT,
    ]


    for directory in directories:

        directory.mkdir(
            parents=True,
            exist_ok=True,
        )


def free_space_gb(path):

    usage = shutil.disk_usage(
        path
    )

    return (
        usage.free
        / 1024**3
    )


# ============================================================
# CHECKSUM HELPERS
# ============================================================

def sha256_file(
    path,
    chunk_size=4 * 1024 * 1024,
):

    digest = hashlib.sha256()


    with open(
        path,
        "rb",
    ) as handle:

        while True:

            block = handle.read(
                chunk_size
            )

            if not block:
                break

            digest.update(
                block
            )


    return digest.hexdigest()


def md5_file(
    path,
    chunk_size=4 * 1024 * 1024,
):

    digest = hashlib.md5()


    with open(
        path,
        "rb",
    ) as handle:

        while True:

            block = handle.read(
                chunk_size
            )

            if not block:
                break

            digest.update(
                block
            )


    return digest.hexdigest()


# ============================================================
# FASTA VALIDATION
# ============================================================

def validate_gzip(path):

    try:

        with gzip.open(
            path,
            "rb",
        ) as handle:

            while handle.read(
                4 * 1024 * 1024
            ):
                pass


        return True


    except Exception:

        return False


def validate_fasta_header(path):

    try:

        with gzip.open(
            path,
            "rt",
            errors="strict",
        ) as handle:

            first_line = (
                handle.readline()
            )


        return (
            first_line.startswith(
                ">"
            )
        )


    except Exception:

        return False


def sequence_stats(path):

    genome_length = 0

    contigs = 0


    with gzip.open(
        path,
        "rt",
        errors="strict",
    ) as handle:

        for line in handle:

            if line.startswith(
                ">"
            ):

                contigs += 1

            else:

                genome_length += len(
                    line.strip()
                )


    return (
        genome_length,
        contigs,
    )


# ============================================================
# COMMAND EXECUTION
# ============================================================

def run_command(
    command,
    log_path,
):

    started = time.time()


    with open(
        log_path,
        "a",
        encoding="utf-8",
    ) as handle:

        handle.write(
            "\n"
            + "=" * 80
            + "\n"
        )

        handle.write(
            f"UTC: {utc_now()}\n"
        )

        handle.write(
            "COMMAND:\n"
        )

        handle.write(
            " ".join(
                map(
                    str,
                    command,
                )
            )
            + "\n"
        )


    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )


    elapsed = (
        time.time()
        - started
    )


    with open(
        log_path,
        "a",
        encoding="utf-8",
    ) as handle:

        handle.write(
            f"RETURN_CODE: "
            f"{result.returncode}\n"
        )

        handle.write(
            f"ELAPSED_SECONDS: "
            f"{elapsed:.2f}\n"
        )

        handle.write(
            "\nSTDOUT:\n"
        )

        handle.write(
            result.stdout
        )

        handle.write(
            "\nSTDERR:\n"
        )

        handle.write(
            result.stderr
        )


    return result


# ============================================================
# PANDAS DTYPE RESTORATION
# ============================================================

def normalize_manifest_dtypes(df):

    """
    CSV does not preserve pandas nullable/string dtypes.

    Empty string-oriented columns may otherwise be interpreted
    as float64 on reload. This caused the FutureWarnings seen
    during pilot registration.
    """


    for column in STRING_COLUMNS:

        if column not in df.columns:

            df[column] = pd.Series(
                pd.NA,
                index=df.index,
                dtype="string",
            )

        else:

            df[column] = (
                df[column]
                .astype(
                    "string"
                )
            )


    for column in INTEGER_COLUMNS:

        if column not in df.columns:

            df[column] = pd.Series(
                pd.NA,
                index=df.index,
                dtype="Int64",
            )

        else:

            df[column] = (
                pd.to_numeric(
                    df[column],
                    errors="coerce",
                )
                .astype(
                    "Int64"
                )
            )


    if (
        "Production_Batch_ID"
        in df.columns
    ):

        df[
            "Production_Batch_ID"
        ] = (
            df[
                "Production_Batch_ID"
            ]
            .astype(
                "string"
            )
        )


    if "Genome_ID" in df.columns:

        df[
            "Genome_ID"
        ] = (
            df[
                "Genome_ID"
            ]
            .astype(
                "string"
            )
        )


    if (
        "Assembly_Accession"
        in df.columns
    ):

        df[
            "Assembly_Accession"
        ] = (
            df[
                "Assembly_Accession"
            ]
            .astype(
                "string"
            )
            .str.strip()
            .str.upper()
        )


    return df


# ============================================================
# MANIFEST SAVE
# ============================================================

def save_manifest(df):

    df = normalize_manifest_dtypes(
        df.copy()
    )


    temporary_path = (
        MASTER_MANIFEST
        .with_suffix(
            ".csv.tmp"
        )
    )


    df.to_csv(
        temporary_path,
        index=False,
    )


    os.replace(
        temporary_path,
        MASTER_MANIFEST,
    )


# ============================================================
# INITIALIZE / LOAD MASTER MANIFEST
# ============================================================

def initialize_master_manifest():

    # --------------------------------------------------------
    # Existing production manifest
    # --------------------------------------------------------

    if MASTER_MANIFEST.exists():

        log(
            "Loading existing production manifest."
        )


        df = pd.read_csv(
            MASTER_MANIFEST,
            low_memory=False,
        )


        df = normalize_manifest_dtypes(
            df
        )


        required_existing = {

            "Genome_ID",

            "Assembly_Accession",

            "Production_Batch_ID",

            "Production_Status",
        }


        missing = (
            required_existing
            .difference(
                df.columns
            )
        )


        if missing:

            raise ValueError(
                "Existing production manifest "
                "is missing columns: "
                + ", ".join(
                    sorted(
                        missing
                    )
                )
            )


        if len(df) != 124081:

            raise ValueError(
                "Existing production manifest "
                "has unexpected row count: "
                f"{len(df)}"
            )


        return df


    # --------------------------------------------------------
    # First production initialization
    # --------------------------------------------------------

    if not SEED_MANIFEST.exists():

        raise FileNotFoundError(
            "Seed manifest missing:\n"
            f"{SEED_MANIFEST}"
        )


    log(
        "Creating production master manifest."
    )


    df = pd.read_csv(
        SEED_MANIFEST,
        low_memory=False,
    )


    if len(df) != 124081:

        raise ValueError(
            "Unexpected Phase5G genome count: "
            f"{len(df)}"
        )


    required = [

        "Genome_ID",

        "Assembly_Accession",
    ]


    for column in required:

        if column not in df.columns:

            raise ValueError(
                "Missing required column: "
                f"{column}"
            )


    df[
        "Genome_ID"
    ] = (
        df[
            "Genome_ID"
        ]
        .astype(
            "string"
        )
    )


    df[
        "Assembly_Accession"
    ] = (
        df[
            "Assembly_Accession"
        ]
        .astype(
            "string"
        )
        .str.strip()
        .str.upper()
    )


    invalid = (
        ~df[
            "Assembly_Accession"
        ]
        .map(
            lambda value:
                bool(
                    ASSEMBLY_RE.fullmatch(
                        str(
                            value
                        )
                    )
                )
        )
    )


    if invalid.any():

        raise ValueError(
            "Invalid assembly accession count: "
            f"{int(invalid.sum())}"
        )


    if not df[
        "Genome_ID"
    ].is_unique:

        raise ValueError(
            "Genome_ID is not unique."
        )


    if not df[
        "Assembly_Accession"
    ].is_unique:

        raise ValueError(
            "Assembly_Accession is not unique."
        )


    # --------------------------------------------------------
    # Immutable fixed batch assignment
    # --------------------------------------------------------

    df[
        "Production_Batch_ID"
    ] = pd.Series(
        [

            f"BATCH_"
            f"{((index // BATCH_SIZE) + 1):04d}"

            for index
            in range(
                len(df)
            )
        ],
        dtype="string",
    )


    # --------------------------------------------------------
    # Production columns
    # --------------------------------------------------------

    for column in STRING_COLUMNS:

        df[column] = pd.Series(
            pd.NA,
            index=df.index,
            dtype="string",
        )


    for column in INTEGER_COLUMNS:

        df[column] = pd.Series(
            pd.NA,
            index=df.index,
            dtype="Int64",
        )


    df[
        "Production_Status"
    ] = "Pending"


    df[
        "Download_Attempts"
    ] = 0


    df[
        "Pipeline_Version"
    ] = PIPELINE_VERSION


    save_manifest(
        df
    )


    return normalize_manifest_dtypes(
        df
    )


# ============================================================
# PILOT REGISTRATION
# ============================================================

PILOT_DIRS = [

    ROOT
    / "genomes"
    / "TEST_10",

    ROOT
    / "genomes"
    / "TEST_100",

    ROOT
    / "genomes"
    / "TEST_1000",
]


def register_existing_pilots(df):

    df = normalize_manifest_dtypes(
        df
    )


    log(
        "Checking previously downloaded "
        "pilot genomes."
    )


    accession_to_index = dict(
        zip(

            df[
                "Assembly_Accession"
            ]
            .astype(str),

            df.index,
        )
    )


    registered = 0

    already_completed = 0

    rejected = 0


    for pilot_dir in PILOT_DIRS:

        if not pilot_dir.exists():

            log(
                "Pilot directory absent; "
                f"skipping: {pilot_dir}"
            )

            continue


        for source_path in sorted(
            pilot_dir.glob(
                "*.fna.gz"
            )
        ):

            accession = (
                source_path.name
                .removesuffix(
                    ".fna.gz"
                )
            )


            if accession not in (
                accession_to_index
            ):

                continue


            index = (
                accession_to_index[
                    accession
                ]
            )


            if (
                df.at[
                    index,
                    "Production_Status",
                ]
                == "Completed"
            ):

                already_completed += 1

                continue


            # ------------------------------------------------
            # Validate pilot source
            # ------------------------------------------------

            if not validate_gzip(
                source_path
            ):

                log(
                    "Pilot invalid gzip: "
                    f"{accession}"
                )

                rejected += 1

                continue


            if not validate_fasta_header(
                source_path
            ):

                log(
                    "Pilot invalid FASTA: "
                    f"{accession}"
                )

                rejected += 1

                continue


            source_sha256 = (
                sha256_file(
                    source_path
                )
            )


            source_md5 = (
                md5_file(
                    source_path
                )
            )


            (
                genome_length,
                contigs,
            ) = sequence_stats(
                source_path
            )


            batch_id = str(
                df.at[
                    index,
                    "Production_Batch_ID",
                ]
            )


            destination_directory = (
                GENOME_ROOT
                / batch_id
            )


            destination_directory.mkdir(
                parents=True,
                exist_ok=True,
            )


            destination = (
                destination_directory
                / f"{accession}.fna.gz"
            )


            # ------------------------------------------------
            # Copy into production hierarchy if absent
            # ------------------------------------------------

            if not destination.exists():

                partial = (
                    destination_directory
                    / (
                        "."
                        + accession
                        + ".fna.gz.partial"
                    )
                )


                partial.unlink(
                    missing_ok=True
                )


                shutil.copy2(
                    source_path,
                    partial,
                )


                copied_sha256 = (
                    sha256_file(
                        partial
                    )
                )


                if (
                    copied_sha256
                    != source_sha256
                ):

                    partial.unlink(
                        missing_ok=True
                    )


                    log(
                        "Pilot copy checksum "
                        "mismatch: "
                        f"{accession}"
                    )


                    rejected += 1

                    continue


                os.replace(
                    partial,
                    destination,
                )


            # ------------------------------------------------
            # Validate final copy
            # ------------------------------------------------

            destination_sha256 = (
                sha256_file(
                    destination
                )
            )


            if (
                destination_sha256
                != source_sha256
            ):

                log(
                    "Existing production pilot "
                    "checksum mismatch: "
                    f"{accession}"
                )


                rejected += 1

                continue


            if not validate_gzip(
                destination
            ):

                log(
                    "Pilot destination gzip "
                    "invalid: "
                    f"{accession}"
                )


                rejected += 1

                continue


            # ------------------------------------------------
            # Register
            # ------------------------------------------------

            df.at[
                index,
                "Production_Status",
            ] = "Completed"


            df.at[
                index,
                "Genome_FASTA_Path",
            ] = str(
                destination
            )


            df.at[
                index,
                "SHA256",
            ] = destination_sha256


            df.at[
                index,
                "MD5",
            ] = source_md5


            df.at[
                index,
                "File_Size_Bytes",
            ] = int(
                destination.stat().st_size
            )


            df.at[
                index,
                "Genome_Length_bp",
            ] = int(
                genome_length
            )


            df.at[
                index,
                "Contig_Count",
            ] = int(
                contigs
            )


            df.at[
                index,
                "Completed_UTC",
            ] = utc_now()


            df.at[
                index,
                "Acquisition_Source",
            ] = "Validated_Pilot"


            df.at[
                index,
                "Failure_Reason",
            ] = pd.NA


            df.at[
                index,
                "Pipeline_Version",
            ] = PIPELINE_VERSION


            registered += 1


            if (
                registered % 25
                == 0
            ):

                save_manifest(
                    df
                )


                log(
                    "Pilot registration "
                    "progress: "
                    f"{registered}"
                )


    save_manifest(
        df
    )


    log(
        "Pilot genomes newly registered: "
        f"{registered}"
    )


    log(
        "Pilot genomes already completed: "
        f"{already_completed}"
    )


    log(
        "Pilot genomes rejected during "
        "validation: "
        f"{rejected}"
    )


    return normalize_manifest_dtypes(
        df
    )


# ============================================================
# RECOVER EXISTING PRODUCTION FILES
# ============================================================

def recover_existing_production_files(df):

    df = normalize_manifest_dtypes(
        df
    )


    log(
        "Checking existing production files."
    )


    recovered = 0


    pending_indices = (
        df.index[
            ~df[
                "Production_Status"
            ].eq(
                "Completed"
            )
        ]
    )


    for index in pending_indices:

        accession = str(
            df.at[
                index,
                "Assembly_Accession",
            ]
        )


        batch_id = str(
            df.at[
                index,
                "Production_Batch_ID",
            ]
        )


        path = (
            GENOME_ROOT
            / batch_id
            / f"{accession}.fna.gz"
        )


        if not path.exists():

            continue


        if not validate_gzip(
            path
        ):

            continue


        if not validate_fasta_header(
            path
        ):

            continue


        sha256 = sha256_file(
            path
        )


        md5 = md5_file(
            path
        )


        (
            genome_length,
            contigs,
        ) = sequence_stats(
            path
        )


        df.at[
            index,
            "Production_Status",
        ] = "Completed"


        df.at[
            index,
            "Genome_FASTA_Path",
        ] = str(
            path
        )


        df.at[
            index,
            "SHA256",
        ] = sha256


        df.at[
            index,
            "MD5",
        ] = md5


        df.at[
            index,
            "File_Size_Bytes",
        ] = int(
            path.stat().st_size
        )


        df.at[
            index,
            "Genome_Length_bp",
        ] = int(
            genome_length
        )


        df.at[
            index,
            "Contig_Count",
        ] = int(
            contigs
        )


        df.at[
            index,
            "Completed_UTC",
        ] = utc_now()


        df.at[
            index,
            "Acquisition_Source",
        ] = (
            "Recovered_Existing_File"
        )


        df.at[
            index,
            "Failure_Reason",
        ] = pd.NA


        df.at[
            index,
            "Pipeline_Version",
        ] = PIPELINE_VERSION


        recovered += 1


        if (
            recovered % 25
            == 0
        ):

            save_manifest(
                df
            )


            log(
                "Existing-file recovery "
                "progress: "
                f"{recovered}"
            )


    save_manifest(
        df
    )


    log(
        "Existing production files "
        "recovered: "
        f"{recovered}"
    )


    return normalize_manifest_dtypes(
        df
    )


# ============================================================
# BATCH STATUS TABLE
# ============================================================

def create_batch_table(df):

    rows = []


    for batch_id, group in (
        df.groupby(
            "Production_Batch_ID",
            sort=True,
        )
    ):

        requested = len(
            group
        )


        completed = int(
            group[
                "Production_Status"
            ]
            .eq(
                "Completed"
            )
            .sum()
        )


        failed = int(
            group[
                "Production_Status"
            ]
            .isin(
                FAILURE_STATUSES
            )
            .sum()
        )


        pending = (
            requested
            - completed
            - failed
        )


        rows.append({

            "Batch_ID":
                batch_id,

            "Requested_N":
                requested,

            "Completed_N":
                completed,

            "Failed_N":
                failed,

            "Pending_N":
                pending,

            "Completion_Percent":
                (
                    100
                    * completed
                    / requested
                ),

            "Batch_Complete":
                (
                    completed
                    == requested
                ),
        })


    table = pd.DataFrame(
        rows
    )


    table.to_csv(
        BATCH_TABLE,
        index=False,
    )


    return table


# ============================================================
# STORAGE PREFLIGHT
# ============================================================

def storage_preflight():

    final_free = (
        free_space_gb(
            ROOT
        )
    )


    local_free = (
        free_space_gb(
            LOCAL_ROOT
        )
    )


    log(
        "Final storage free: "
        f"{final_free:.1f} GB"
    )


    log(
        "Local NVMe free: "
        f"{local_free:.1f} GB"
    )


    if (
        final_free
        <
        MIN_FINAL_FREE_GB
    ):

        raise RuntimeError(
            "Final storage below "
            "safety threshold."
        )


    if (
        local_free
        <
        MIN_LOCAL_FREE_GB
    ):

        raise RuntimeError(
            "Local staging storage below "
            "safety threshold."
        )


# ============================================================
# WRITE ACCESSION FILE
# ============================================================

def write_batch_accessions(
    group,
    accession_file,
):

    accessions = (
        group[
            "Assembly_Accession"
        ]
        .astype(str)
        .tolist()
    )


    accession_file.write_text(
        "\n".join(
            accessions
        )
        + "\n",
        encoding="utf-8",
    )


# ============================================================
# PREPARE NCBI PACKAGE
# ============================================================

def prepare_batch_package(
    batch_id,
    pending_group,
):

    package_zip = (
        LOCAL_PACKAGE_ROOT
        / f"{batch_id}.zip"
    )


    extracted_dir = (
        LOCAL_REHYDRATE_ROOT
        / batch_id
    )


    accession_file = (
        LOCAL_ACCESSION_ROOT
        / f"{batch_id}.txt"
    )


    log_path = (
        LOG_ROOT
        / f"{batch_id}.log"
    )


    # --------------------------------------------------------
    # Remove incomplete staging from a prior attempt
    # --------------------------------------------------------

    package_zip.unlink(
        missing_ok=True
    )


    if extracted_dir.exists():

        shutil.rmtree(
            extracted_dir
        )


    write_batch_accessions(
        pending_group,
        accession_file,
    )


    # --------------------------------------------------------
    # Dehydrated package
    # --------------------------------------------------------

    command = [

        str(
            DATASETS_BIN
        ),

        "download",

        "genome",

        "accession",

        "--inputfile",

        str(
            accession_file
        ),

        "--include",

        "genome",

        "--dehydrated",

        "--filename",

        str(
            package_zip
        ),

        "--no-progressbar",
    ]


    success = False


    for attempt in range(
        1,
        MAX_RETRIES + 1,
    ):

        log(
            f"{batch_id}: "
            "dehydrated download attempt "
            f"{attempt}"
        )


        result = run_command(
            command,
            log_path,
        )


        if (
            result.returncode
            == 0
            and package_zip.exists()
            and package_zip.stat().st_size
            > 0
        ):

            success = True

            break


        time.sleep(
            RETRY_WAIT_SECONDS
        )


    if not success:

        raise RuntimeError(
            "Dehydrated package "
            "download failed."
        )


    # --------------------------------------------------------
    # Extract package
    # --------------------------------------------------------

    extracted_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    result = subprocess.run(

        [

            "unzip",

            "-q",

            str(
                package_zip
            ),

            "-d",

            str(
                extracted_dir
            ),
        ],

        capture_output=True,

        text=True,
    )


    if result.returncode != 0:

        raise RuntimeError(
            "ZIP extraction failed: "
            + result.stderr
        )


    # --------------------------------------------------------
    # Rehydrate genome files directly as gzip
    # --------------------------------------------------------

    command = [

        str(
            DATASETS_BIN
        ),

        "rehydrate",

        "--directory",

        str(
            extracted_dir
        ),

        "--gzip",

        "--max-workers",

        str(
            MAX_WORKERS
        ),

        "--no-progressbar",
    ]


    success = False


    for attempt in range(
        1,
        MAX_RETRIES + 1,
    ):

        log(
            f"{batch_id}: "
            "rehydrate attempt "
            f"{attempt}"
        )


        result = run_command(
            command,
            log_path,
        )


        if (
            result.returncode
            == 0
        ):

            success = True

            break


        time.sleep(
            RETRY_WAIT_SECONDS
        )


    if not success:

        raise RuntimeError(
            "Rehydration failed."
        )


    return (

        package_zip,

        extracted_dir,

        accession_file,

        log_path,
    )


# ============================================================
# PROCESS ONE GENOME
# ============================================================

def process_one_genome(
    df,
    index,
    extracted_dir,
):

    accession = str(
        df.at[
            index,
            "Assembly_Accession",
        ]
    )


    batch_id = str(
        df.at[
            index,
            "Production_Batch_ID",
        ]
    )


    source_dir = (

        extracted_dir

        / "ncbi_dataset"

        / "data"

        / accession
    )


    candidates = sorted(
        source_dir.glob(
            "*.fna.gz"
        )
    )


    df.at[
        index,
        "Last_Attempt_UTC",
    ] = utc_now()


    attempt_value = (
        pd.to_numeric(

            pd.Series(
                [
                    df.at[
                        index,
                        "Download_Attempts",
                    ]
                ]
            ),

            errors="coerce",
        )
        .fillna(0)
        .iloc[0]
    )


    df.at[
        index,
        "Download_Attempts",
    ] = int(
        attempt_value
    ) + 1


    # --------------------------------------------------------
    # Missing from NCBI package
    # --------------------------------------------------------

    if not candidates:

        df.at[
            index,
            "Production_Status",
        ] = (
            "Missing_From_NCBI"
        )


        df.at[
            index,
            "Failure_Reason",
        ] = (
            "No *.fna.gz file found "
            "after successful rehydration."
        )


        return False


    source = candidates[0]


    # --------------------------------------------------------
    # Source validation
    # --------------------------------------------------------

    if not validate_gzip(
        source
    ):

        df.at[
            index,
            "Production_Status",
        ] = "Invalid_Gzip"


        df.at[
            index,
            "Failure_Reason",
        ] = (
            "Source gzip integrity "
            "test failed."
        )


        return False


    if not validate_fasta_header(
        source
    ):

        df.at[
            index,
            "Production_Status",
        ] = "Invalid_FASTA"


        df.at[
            index,
            "Failure_Reason",
        ] = (
            "Source file did not begin "
            "with a FASTA header."
        )


        return False


    source_sha256 = (
        sha256_file(
            source
        )
    )


    # --------------------------------------------------------
    # Permanent batch directory
    # --------------------------------------------------------

    final_dir = (
        GENOME_ROOT
        / batch_id
    )


    final_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    destination = (

        final_dir

        / f"{accession}.fna.gz"
    )


    temporary_destination = (

        final_dir

        / (
            "."
            + accession
            + ".fna.gz.partial"
        )
    )


    temporary_destination.unlink(
        missing_ok=True
    )


    # --------------------------------------------------------
    # Copy from local NVMe to CIFS share
    # --------------------------------------------------------

    # Copy file contents only.
    #
    # Do not use shutil.copy2() here. copy2() additionally
    # preserves filesystem metadata/timestamps via copystat(),
    # which can trigger CIFS SET_FILE_INFO/utimensat operations.
    # Sequence integrity is independently verified by SHA256.
    shutil.copyfile(
        source,
        temporary_destination,
    )


    copied_sha256 = (
        sha256_file(
            temporary_destination
        )
    )


    # --------------------------------------------------------
    # Source/destination checksum comparison
    # --------------------------------------------------------

    if (
        copied_sha256
        != source_sha256
    ):

        temporary_destination.unlink(
            missing_ok=True
        )


        df.at[
            index,
            "Production_Status",
        ] = "Checksum_Failed"


        df.at[
            index,
            "Failure_Reason",
        ] = (
            "SHA256 mismatch between "
            "local staging and CIFS copy."
        )


        return False


    # Atomic rename within destination filesystem.
    os.replace(
        temporary_destination,
        destination,
    )


    # --------------------------------------------------------
    # Destination validation
    # --------------------------------------------------------

    if not validate_gzip(
        destination
    ):

        destination.unlink(
            missing_ok=True
        )


        df.at[
            index,
            "Production_Status",
        ] = "Invalid_Gzip"


        df.at[
            index,
            "Failure_Reason",
        ] = (
            "Destination gzip "
            "validation failed."
        )


        return False


    if not validate_fasta_header(
        destination
    ):

        destination.unlink(
            missing_ok=True
        )


        df.at[
            index,
            "Production_Status",
        ] = "Invalid_FASTA"


        df.at[
            index,
            "Failure_Reason",
        ] = (
            "Destination FASTA header "
            "validation failed."
        )


        return False


    (
        genome_length,
        contigs,
    ) = sequence_stats(
        destination
    )


    # --------------------------------------------------------
    # Mark completed
    # --------------------------------------------------------

    df.at[
        index,
        "Production_Status",
    ] = "Completed"


    df.at[
        index,
        "Genome_FASTA_Path",
    ] = str(
        destination
    )


    df.at[
        index,
        "SHA256",
    ] = copied_sha256


    df.at[
        index,
        "MD5",
    ] = md5_file(
        destination
    )


    df.at[
        index,
        "File_Size_Bytes",
    ] = int(
        destination.stat().st_size
    )


    df.at[
        index,
        "Genome_Length_bp",
    ] = int(
        genome_length
    )


    df.at[
        index,
        "Contig_Count",
    ] = int(
        contigs
    )


    df.at[
        index,
        "Completed_UTC",
    ] = utc_now()


    df.at[
        index,
        "Acquisition_Source",
    ] = (
        "NCBI_Datasets_Production"
    )


    df.at[
        index,
        "Failure_Reason",
    ] = pd.NA


    df.at[
        index,
        "Pipeline_Version",
    ] = PIPELINE_VERSION


    return True


# ============================================================
# PROCESS ONE BATCH
# ============================================================

def process_batch(
    df,
    batch_id,
    max_genomes=None,
):

    df = normalize_manifest_dtypes(
        df
    )


    storage_preflight()


    batch_mask = (
        df[
            "Production_Batch_ID"
        ].eq(
            batch_id
        )
    )


    pending_mask = (

        batch_mask

        &

        ~df[
            "Production_Status"
        ].eq(
            "Completed"
        )
    )


    pending_group = (
        df.loc[
            pending_mask
        ]
        .copy()
    )


    if (
        max_genomes
        is not None
    ):

        if max_genomes <= 0:

            raise ValueError(
                "--max-genomes must be "
                "greater than zero."
            )

        pending_group = (
            pending_group
            .head(
                max_genomes
            )
            .copy()
        )


    if pending_group.empty:

        log(
            f"{batch_id}: "
            "already complete."
        )


        return df


    log(
        f"{batch_id}: "
        f"{len(pending_group)} "
        "genomes pending."
    )


    # --------------------------------------------------------
    # Package-level acquisition
    # --------------------------------------------------------

    try:

        (
            package_zip,

            extracted_dir,

            accession_file,

            log_path,

        ) = prepare_batch_package(

            batch_id,

            pending_group,
        )


    except Exception as exc:

        log(
            f"{batch_id}: "
            "package-level failure: "
            f"{exc}"
        )


        for index in (
            pending_group.index
        ):

            attempt_value = (
                pd.to_numeric(

                    pd.Series(
                        [
                            df.at[
                                index,
                                "Download_Attempts",
                            ]
                        ]
                    ),

                    errors="coerce",
                )
                .fillna(0)
                .iloc[0]
            )


            df.at[
                index,
                "Download_Attempts",
            ] = int(
                attempt_value
            ) + 1


            df.at[
                index,
                "Production_Status",
            ] = "Failed"


            df.at[
                index,
                "Last_Attempt_UTC",
            ] = utc_now()


            df.at[
                index,
                "Failure_Reason",
            ] = (
                "Batch-level "
                "download/rehydration "
                f"failure: {exc}"
            )


        save_manifest(
            df
        )


        return (
            normalize_manifest_dtypes(
                df
            )
        )


    # --------------------------------------------------------
    # Genome-level processing
    # --------------------------------------------------------

    completed = 0

    failed = 0

    total = len(
        pending_group
    )


    for number, index in enumerate(

        pending_group.index,

        start=1,
    ):

        success = process_one_genome(

            df,

            index,

            extracted_dir,
        )


        if success:

            completed += 1

        else:

            failed += 1


        # Reduced CIFS manifest checkpointing.
        #
        # Completed genome files are durable on final storage
        # and can be reconstructed with --recover-existing.
        # Avoid repeatedly rewriting the full master CSV over CIFS.
        if number == total:

            save_manifest(
                df
            )


            log(
                f"{batch_id}: "
                f"{number}/{total} processed; "
                f"completed={completed}; "
                f"failed={failed}"
            )


    # --------------------------------------------------------
    # Batch QC output
    # --------------------------------------------------------

    batch_qc_columns = [

        "Genome_ID",

        "Assembly_Accession",

        "Production_Batch_ID",

        "Production_Status",

        "Genome_FASTA_Path",

        "File_Size_Bytes",

        "Genome_Length_bp",

        "Contig_Count",

        "SHA256",

        "Download_Attempts",

        "Failure_Reason",
    ]


    batch_qc = (
        df.loc[
            batch_mask,
            batch_qc_columns,
        ]
        .copy()
    )


    batch_qc.to_csv(

        QC_ROOT
        / f"{batch_id}_QC.csv",

        index=False,
    )


    create_batch_table(
        df
    )


    # --------------------------------------------------------
    # Clean local staging
    # --------------------------------------------------------

    try:

        package_zip.unlink(
            missing_ok=True
        )


        accession_file.unlink(
            missing_ok=True
        )


        if extracted_dir.exists():

            shutil.rmtree(
                extracted_dir
            )


    except Exception as exc:

        log(
            f"{batch_id}: "
            "staging cleanup warning: "
            f"{exc}"
        )


    log(
        f"{batch_id}: "
        "batch finished. "
        f"completed={completed}; "
        f"failed={failed}"
    )


    return (
        normalize_manifest_dtypes(
            df
        )
    )


# ============================================================
# RETRY FAILED RECORDS
# ============================================================

def reset_retryable_failures(df):

    df = normalize_manifest_dtypes(
        df
    )


    mask = (

        df[
            "Production_Status"
        ].isin(
            FAILURE_STATUSES
        )

        &

        df[
            "Download_Attempts"
        ]
        .fillna(0)
        .lt(
            MAX_RETRIES
        )
    )


    count = int(
        mask.sum()
    )


    if count:

        df.loc[
            mask,
            "Production_Status",
        ] = "Pending_Retry"


    save_manifest(
        df
    )


    log(
        "Retryable records reset: "
        f"{count}"
    )


    return (
        normalize_manifest_dtypes(
            df
        )
    )


# ============================================================
# STATUS SUMMARY
# ============================================================

def write_summary(df):

    df = normalize_manifest_dtypes(
        df
    )


    status_counts = (

        df[
            "Production_Status"
        ]
        .value_counts(
            dropna=False
        )
        .rename_axis(
            "Status"
        )
        .reset_index(
            name="Genome_N"
        )
    )


    status_counts.to_csv(

        QC_ROOT
        / (
            "Phase6_Production_"
            "Status_Summary.csv"
        ),

        index=False,
    )


    completed = (
        df.loc[
            df[
                "Production_Status"
            ].eq(
                "Completed"
            )
        ]
        .copy()
    )


    total_bytes = (

        pd.to_numeric(

            completed[
                "File_Size_Bytes"
            ],

            errors="coerce",
        )
        .fillna(0)
        .sum()
    )


    summary = {

        "phase":
            "Phase6_Genome_Acquisition",

        "pipeline_version":
            PIPELINE_VERSION,

        "timestamp_utc":
            utc_now(),

        "input_genomes":
            int(
                len(df)
            ),

        "completed_genomes":
            int(
                len(completed)
            ),

        "remaining_genomes":
            int(
                len(df)
                - len(completed)
            ),

        "completion_percent":
            float(
                100
                * len(completed)
                / len(df)
            ),

        "compressed_archive_gib":
            float(
                total_bytes
                / 1024**3
            ),

        "final_storage_free_gb":
            float(
                free_space_gb(
                    ROOT
                )
            ),

        "local_storage_free_gb":
            float(
                free_space_gb(
                    LOCAL_ROOT
                )
            ),
    }


    with open(

        FREEZE_ROOT
        / (
            "Phase6_Production_"
            "Run_Summary.json"
        ),

        "w",

        encoding="utf-8",

    ) as handle:

        json.dump(
            summary,
            handle,
            indent=2,
        )


    log(
        "Current production status:"
    )


    print(
        status_counts.to_string(
            index=False
        )
    )


    log(
        "Completed: "
        f"{len(completed):,}/"
        f"{len(df):,} "
        f"("
        f"{summary['completion_percent']:.2f}"
        f"%)"
    )


    log(
        "Compressed archive: "
        f"{summary['compressed_archive_gib']:.2f} "
        "GiB"
    )


# ============================================================
# FAILURE / PENDING LEDGER
# ============================================================

def write_failure_ledger(df):

    unresolved = (
        df.loc[
            ~df[
                "Production_Status"
            ].eq(
                "Completed"
            )
        ]
        .copy()
    )


    unresolved.to_csv(

        FAILED_ROOT
        / "Phase6_Unresolved_Downloads.csv",

        index=False,
    )


# ============================================================
# ARGUMENT PARSER
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(

        description=(
            "Phase 6 production "
            "genome acquisition"
        )
    )


    parser.add_argument(

        "--start-batch",

        type=int,

        default=1,

        help=(
            "First batch number "
            "to consider."
        ),
    )


    parser.add_argument(

        "--end-batch",

        type=int,

        default=None,

        help=(
            "Last batch number "
            "to consider."
        ),
    )


    parser.add_argument(

        "--max-batches",

        type=int,

        default=None,

        help=(
            "Maximum number of "
            "incomplete batches to "
            "process in this run."
        ),
    )


    parser.add_argument(

        "--register-pilots",

        action="store_true",

        help=(
            "Validate/register existing "
            "TEST_10, TEST_100 and "
            "TEST_1000 genomes."
        ),
    )


    parser.add_argument(

        "--max-genomes",

        type=int,

        default=None,

        help=(
            "Maximum number of pending genomes "
            "to process per selected batch."
        ),
    )


    parser.add_argument(

        "--recover-existing",

        action="store_true",

        help=(
            "Recover valid existing "
            "production files not yet "
            "marked complete."
        ),
    )


    parser.add_argument(

        "--retry-failures",

        action="store_true",

        help=(
            "Reset retryable failed "
            "records."
        ),
    )


    parser.add_argument(

        "--summary-only",

        action="store_true",

        help=(
            "Print/write current status "
            "without downloading."
        ),
    )


    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()


    ensure_directories()


    # --------------------------------------------------------
    # NCBI datasets CLI validation
    # --------------------------------------------------------

    if not DATASETS_BIN.exists():

        raise FileNotFoundError(
            "datasets executable "
            "not found: "
            f"{DATASETS_BIN}"
        )


    version = subprocess.run(

        [
            str(
                DATASETS_BIN
            ),
            "version",
        ],

        capture_output=True,

        text=True,
    )


    if version.returncode != 0:

        raise RuntimeError(
            "datasets version check failed:\n"
            + version.stderr
        )


    log(
        version.stdout.strip()
    )


    # --------------------------------------------------------
    # Storage check
    # --------------------------------------------------------

    storage_preflight()


    # --------------------------------------------------------
    # Master manifest
    # --------------------------------------------------------

    df = initialize_master_manifest()


    # --------------------------------------------------------
    # Optional pilot registration
    # --------------------------------------------------------

    if args.register_pilots:

        df = register_existing_pilots(
            df
        )


    # --------------------------------------------------------
    # Optional existing-file recovery
    # --------------------------------------------------------

    if args.recover_existing:

        df = (
            recover_existing_production_files(
                df
            )
        )


    # --------------------------------------------------------
    # Optional failed-record retry reset
    # --------------------------------------------------------

    if args.retry_failures:

        df = reset_retryable_failures(
            df
        )


    create_batch_table(
        df
    )


    # --------------------------------------------------------
    # Summary-only mode
    # --------------------------------------------------------

    if args.summary_only:

        save_manifest(
            df
        )

        write_summary(
            df
        )

        write_failure_ledger(
            df
        )

        return


    # ========================================================
    # BUILD BATCH LIST
    # ========================================================

    batch_numbers = sorted(
        {

            int(
                str(
                    value
                )
                .split(
                    "_"
                )[1]
            )

            for value
            in df[
                "Production_Batch_ID"
            ]
            .dropna()
            .unique()
        }
    )


    batch_numbers = [

        number

        for number
        in batch_numbers

        if (
            number
            >= args.start_batch
        )

        and

        (
            args.end_batch
            is None

            or

            number
            <= args.end_batch
        )
    ]


    selected_batches = []


    for number in batch_numbers:

        batch_id = (
            f"BATCH_{number:04d}"
        )


        batch = (
            df.loc[
                df[
                    "Production_Batch_ID"
                ].eq(
                    batch_id
                )
            ]
        )


        if (
            batch[
                "Production_Status"
            ]
            .eq(
                "Completed"
            )
            .all()
        ):

            continue


        selected_batches.append(
            batch_id
        )


        if (
            args.max_batches
            is not None

            and

            len(
                selected_batches
            )
            >= args.max_batches
        ):

            break


    log(
        "Selected incomplete batches: "
        f"{len(selected_batches)}"
    )


    if selected_batches:

        log(
            "First selected batch: "
            f"{selected_batches[0]}"
        )


        log(
            "Last selected batch: "
            f"{selected_batches[-1]}"
        )


    # ========================================================
    # PROCESS BATCHES
    # ========================================================

    for batch_id in selected_batches:

        try:

            df = process_batch(
                df,
                batch_id,
                max_genomes=args.max_genomes,
            )


        except KeyboardInterrupt:

            log(
                "Keyboard interrupt received."
            )


            save_manifest(
                df
            )


            create_batch_table(
                df
            )


            write_summary(
                df
            )


            write_failure_ledger(
                df
            )


            raise


        except Exception as exc:

            log(
                "Unexpected batch exception "
                f"{batch_id}: "
                f"{exc}"
            )


            save_manifest(
                df
            )


    # ========================================================
    # FINAL CHECKPOINT
    # ========================================================

    save_manifest(
        df
    )


    create_batch_table(
        df
    )


    write_summary(
        df
    )


    write_failure_ledger(
        df
    )


    log(
        "Phase 6 production run finished."
    )


if __name__ == "__main__":

    main()
