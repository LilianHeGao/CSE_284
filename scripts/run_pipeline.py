#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_env_file(path: Path) -> dict[str, str]:
    env: dict[str, str] = {}
    if not path.exists():
        return env

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        env[key.strip()] = normalize_config_value(value)
    return env


def normalize_config_value(value: str) -> str:
    if os.name != "nt":
        return value

    match = re.match(r"^/([a-zA-Z])/(.*)$", value)
    if not match:
        return value

    drive = match.group(1).upper()
    rest = match.group(2)
    return f"{drive}:/{rest}"


def load_config(config_arg: str) -> tuple[dict[str, str], Path]:
    config_path = Path(config_arg)
    if not config_path.is_absolute():
        config_path = REPO_ROOT / config_path
    config_path = config_path.resolve()

    cfg = parse_env_file(config_path)
    local_cfg = parse_env_file(Path(f"{config_path}.local"))
    cfg.update(local_cfg)
    return cfg, config_path


def is_windows_abs(path_str: str) -> bool:
    return bool(re.match(r"^[a-zA-Z]:[/\\]", path_str))


def looks_like_path(value: str) -> bool:
    return (
        "/" in value
        or "\\" in value
        or value.startswith(".")
        or is_windows_abs(value)
    )


def repo_path(path_str: str) -> Path:
    if is_windows_abs(path_str):
        return Path(path_str)
    path = Path(path_str)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def command_exists(command: str) -> bool:
    if looks_like_path(command):
        return repo_path(command).exists()
    return shutil.which(command) is not None


def ensure_dir(path_str: str) -> None:
    repo_path(path_str).mkdir(parents=True, exist_ok=True)


def ensure_file(path_str: str, hint: str | None = None) -> None:
    path = repo_path(path_str)
    if path.exists():
        return

    msg = f"Required file not found: {path_str}"
    if hint:
        msg = f"{msg}. {hint}"
    raise SystemExit(msg)


def ensure_bfile(prefix: str, hint: str | None = None) -> None:
    missing = []
    for ext in (".bed", ".bim", ".fam"):
        path = repo_path(f"{prefix}{ext}")
        if not path.exists():
            missing.append(path)

    if not missing:
        return

    msg = f"Required PLINK bfile not found for prefix '{prefix}'"
    if hint:
        msg = f"{msg}. {hint}"
    raise SystemExit(msg)


def run_command(args: list[str]) -> None:
    print(">", " ".join(args))
    try:
        subprocess.run(args, cwd=REPO_ROOT, check=True)
    except FileNotFoundError as exc:
        raise SystemExit(
            f"Command not found: {args[0]}. Set its path in config/project.env.local."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"Command failed with exit code {exc.returncode}: {' '.join(args)}") from exc


def move_output(src_rel: str, dest_rel: str) -> None:
    src = repo_path(src_rel)
    dest = repo_path(dest_rel)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not src.exists():
        raise FileNotFoundError(f"Expected output file not found: {src}")
    if dest.exists():
        dest.unlink()
    shutil.move(str(src), str(dest))


def build_lmm_covars(eigenvec_rel: str, covars_rel: str) -> None:
    eigenvec = repo_path(eigenvec_rel)
    covars = repo_path(covars_rel)
    covars.parent.mkdir(parents=True, exist_ok=True)

    with eigenvec.open("r", encoding="utf-8") as src, covars.open("w", encoding="utf-8") as dst:
        for idx, line in enumerate(src):
            fields = line.strip().split()
            if not fields:
                continue

            if idx == 0 and len(fields) >= 3 and fields[2].upper().startswith("PC"):
                continue

            dst.write("\t".join(fields[2:]) + "\n")


def build_causal_snplist(bfile_prefix: str, out_rel: str, count: int = 50) -> None:
    bim_path = repo_path(f"{bfile_prefix}.bim")
    out_path = repo_path(out_rel)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    snps: list[str] = []
    with bim_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().split()
            if len(fields) >= 2:
                snps.append(fields[1])

    chosen = snps[:count]
    with out_path.open("w", encoding="utf-8") as handle:
        for snp in chosen:
            handle.write(f"{snp}\n")


def context_from_config(cfg: dict[str, str]) -> dict[str, str]:
    ctx = dict(cfg)

    ctx.setdefault("RAW_DIR", "data/raw")
    ctx.setdefault("PROC_DIR", "data/processed")
    ctx.setdefault("RESULTS_DIR", "results")
    ctx.setdefault("PHENO_SOURCE_CSV", f"{ctx['RAW_DIR']}/values.csv")
    ctx.setdefault("PHENO_NAME", "PHENO")
    ctx.setdefault("SUBSET_NAME", "all")
    ctx.setdefault("RUN_LABEL", ctx["SUBSET_NAME"])

    run_label = ctx["RUN_LABEL"]
    chr_name = ctx["CHR"]

    ctx.setdefault("PHENO_FILE", f"{ctx['PROC_DIR']}/{run_label}.pheno.txt")
    ctx.setdefault("KEEP_IDS_FILE", f"{ctx['PROC_DIR']}/{run_label}.keep_ids.txt")
    ctx.setdefault("PHENO_SAMPLES_FILE", f"{ctx['PROC_DIR']}/{run_label}.samples.tsv")
    ctx.setdefault("PFILE_PREFIX", f"{ctx['PROC_DIR']}/chr{chr_name}.{run_label}")
    ctx.setdefault("BFILE_PREFIX", f"{ctx['PROC_DIR']}/chr{chr_name}.{run_label}.qc")
    ctx.setdefault("CLEAN_VCF_GZ", f"{ctx['RAW_DIR']}/1001G.Chr{chr_name}_clean.vcf.gz")
    ctx.setdefault("BCFTOOLS", "bcftools")
    ctx.setdefault("PLINK2", "plink2")
    ctx.setdefault("GEMMA", "gemma")
    ctx.setdefault("GCTA", "gcta64")
    ctx.setdefault("NEW_ID_MAX_ALLELE_LEN", "40")

    ctx.setdefault("LR_PREFIX", f"{ctx['RESULTS_DIR']}/lr/{run_label}_lr")
    ctx.setdefault("PCA_PREFIX", f"{ctx['RESULTS_DIR']}/lr_pcs/{run_label}_pca")
    ctx.setdefault("LR_PCS_PREFIX", f"{ctx['RESULTS_DIR']}/lr_pcs/{run_label}_lr_pcs")
    ctx.setdefault("LMM_PCA_PREFIX", f"{ctx['RESULTS_DIR']}/lmm/{run_label}_pca")
    ctx.setdefault("LMM_COVARS", f"{ctx['RESULTS_DIR']}/lmm/{run_label}_lmm_covars.txt")
    ctx.setdefault("KINSHIP_PREFIX", f"{run_label}_kinship")
    ctx.setdefault("KINSHIP_OUT", f"{ctx['RESULTS_DIR']}/lmm/{run_label}_kinship.cXX.txt")
    ctx.setdefault("LMM_PREFIX", f"{run_label}_lmm")
    ctx.setdefault("LMM_OUT", f"{ctx['RESULTS_DIR']}/lmm/{run_label}_lmm.assoc.txt")
    ctx.setdefault("EVAL_OUT_PREFIX", f"{ctx['RESULTS_DIR']}/plots/{run_label}_benchmark")

    return ctx


def cmd_check_tools(ctx: dict[str, str]) -> None:
    required = {
        "plink2": ctx["PLINK2"],
        "python": sys.executable,
    }
    optional = {
        "bcftools (VCF normalization)": ctx["BCFTOOLS"],
        "gemma (LMM step)": ctx["GEMMA"],
        "gcta64 (simulation step)": ctx["GCTA"],
    }

    missing = [label for label, tool in required.items() if not command_exists(tool)]
    for label, tool in optional.items():
        if not command_exists(tool):
            print(f"Optional dependency not found: {label} ({tool})")

    if missing:
        raise SystemExit(f"Missing dependencies: {', '.join(missing)}")

    print("All required commands found.")


def cmd_prepare_data(ctx: dict[str, str]) -> None:
    ensure_dir(ctx["RAW_DIR"])
    ensure_dir(ctx["PROC_DIR"])
    ensure_dir(ctx["RESULTS_DIR"])
    ensure_file(
        ctx["VCF_GZ"],
        "Set VCF_GZ in config/project.env.local or place the VCF at this path.",
    )

    run_command(
        [
            sys.executable,
            "scripts/00_build_subset_inputs.py",
            "--values-csv",
            ctx["PHENO_SOURCE_CSV"],
            "--subset-name",
            ctx["SUBSET_NAME"],
            "--countries",
            ctx.get("SUBSET_COUNTRIES", ""),
            "--pheno-name",
            ctx["PHENO_NAME"],
            "--out-pheno",
            ctx["PHENO_FILE"],
            "--out-keep",
            ctx["KEEP_IDS_FILE"],
            "--out-samples",
            ctx["PHENO_SAMPLES_FILE"],
        ]
    )

    vcf_for_plink = ctx["CLEAN_VCF_GZ"]
    if command_exists(ctx["BCFTOOLS"]):
        run_command(
            [
                ctx["BCFTOOLS"],
                "view",
                "-v",
                "snps",
                "-m2",
                "-M2",
                ctx["VCF_GZ"],
                "-Oz",
                "-o",
                ctx["CLEAN_VCF_GZ"],
            ]
        )
    else:
        print("bcftools not found; sanitizing VCF with Python fallback.")
        run_command(
            [
                sys.executable,
                "scripts/sanitize_vcf.py",
                "--in-vcf",
                ctx["VCF_GZ"],
                "--out-vcf",
                ctx["CLEAN_VCF_GZ"],
            ]
        )

    run_command(
        [
            ctx["PLINK2"],
            "--vcf",
            vcf_for_plink,
            "--keep",
            ctx["KEEP_IDS_FILE"],
            "--chr-set",
            "5",
            "--allow-extra-chr",
            "--new-id-max-allele-len",
            ctx["NEW_ID_MAX_ALLELE_LEN"],
            "--set-all-var-ids",
            "@:#:$r:$a",
            "--make-pgen",
            "--out",
            ctx["PFILE_PREFIX"],
        ]
    )

    run_command(
        [
            ctx["PLINK2"],
            "--pfile",
            ctx["PFILE_PREFIX"],
            "--maf",
            ctx["MAF"],
            "--geno",
            ctx["GENO_MISSING"],
            "--make-bed",
            "--out",
            ctx["BFILE_PREFIX"],
        ]
    )


def cmd_run_lr(ctx: dict[str, str]) -> None:
    ensure_dir(f"{ctx['RESULTS_DIR']}/lr")
    ensure_dir(f"{ctx['RESULTS_DIR']}/lr_pcs")
    ensure_bfile(ctx["BFILE_PREFIX"], "Run 'prepare-data' first.")

    run_command(
        [
            ctx["PLINK2"],
            "--bfile",
            ctx["BFILE_PREFIX"],
            "--pheno",
            ctx["PHENO_FILE"],
            "--pheno-name",
            ctx["PHENO_NAME"],
            "--glm",
            "hide-covar",
            "allow-no-covars",
            "--out",
            ctx["LR_PREFIX"],
        ]
    )

    run_command(
        [
            ctx["PLINK2"],
            "--bfile",
            ctx["BFILE_PREFIX"],
            "--pca",
            ctx["NUM_PCS"],
            "--out",
            ctx["PCA_PREFIX"],
        ]
    )

    run_command(
        [
            ctx["PLINK2"],
            "--bfile",
            ctx["BFILE_PREFIX"],
            "--pheno",
            ctx["PHENO_FILE"],
            "--pheno-name",
            ctx["PHENO_NAME"],
            "--covar",
            f"{ctx['PCA_PREFIX']}.eigenvec",
            "--glm",
            "hide-covar",
            "--out",
            ctx["LR_PCS_PREFIX"],
        ]
    )


def cmd_run_lmm(ctx: dict[str, str]) -> None:
    ensure_dir(f"{ctx['RESULTS_DIR']}/lmm")
    ensure_bfile(ctx["BFILE_PREFIX"], "Run 'prepare-data' first.")
    if not command_exists(ctx["GEMMA"]):
        raise SystemExit(
            "GEMMA not found. Set GEMMA=<path> in config/project.env.local or skip run-lmm."
        )

    run_command(
        [
            ctx["GEMMA"],
            "-bfile",
            ctx["BFILE_PREFIX"],
            "-gk",
            "1",
            "-o",
            ctx["KINSHIP_PREFIX"],
        ]
    )
    move_output(f"output/{ctx['KINSHIP_PREFIX']}.cXX.txt", ctx["KINSHIP_OUT"])

    run_command(
        [
            ctx["PLINK2"],
            "--bfile",
            ctx["BFILE_PREFIX"],
            "--pca",
            ctx["NUM_PCS"],
            "--out",
            ctx["LMM_PCA_PREFIX"],
        ]
    )
    build_lmm_covars(f"{ctx['LMM_PCA_PREFIX']}.eigenvec", ctx["LMM_COVARS"])

    run_command(
        [
            ctx["GEMMA"],
            "-bfile",
            ctx["BFILE_PREFIX"],
            "-k",
            ctx["KINSHIP_OUT"],
            "-lmm",
            "4",
            "-c",
            ctx["LMM_COVARS"],
            "-o",
            ctx["LMM_PREFIX"],
        ]
    )
    move_output(f"output/{ctx['LMM_PREFIX']}.assoc.txt", ctx["LMM_OUT"])


def cmd_evaluate(ctx: dict[str, str]) -> None:
    ensure_dir(f"{ctx['RESULTS_DIR']}/plots")
    ensure_file(
        f"{ctx['LR_PREFIX']}.PHENO1.glm.linear",
        "Run 'run-lr' first.",
    )
    ensure_file(
        f"{ctx['LR_PCS_PREFIX']}.PHENO1.glm.linear",
        "Run 'run-lr' first.",
    )
    eval_cmd = [
        sys.executable,
        "scripts/05_evaluate.py",
        "--lr",
        f"{ctx['LR_PREFIX']}.PHENO1.glm.linear",
        "--lr-pcs",
        f"{ctx['LR_PCS_PREFIX']}.PHENO1.glm.linear",
        "--out-prefix",
        ctx["EVAL_OUT_PREFIX"],
    ]

    if repo_path(ctx["LMM_OUT"]).exists():
        eval_cmd.extend(["--lmm", ctx["LMM_OUT"]])
    else:
        print(f"LMM file not found at {ctx['LMM_OUT']}; evaluating LR and LR+PCs only.")

    run_command(eval_cmd)


def cmd_simulate(ctx: dict[str, str]) -> None:
    ensure_dir(f"{ctx['RESULTS_DIR']}/sim")
    ensure_bfile(ctx["BFILE_PREFIX"], "Run 'prepare-data' first.")

    if not command_exists(ctx["GCTA"]):
        raise SystemExit(f"Optional dependency not found: gcta64 ({ctx['GCTA']})")

    causal_snplist = f"{ctx['RESULTS_DIR']}/sim/causal.snplist"
    build_causal_snplist(ctx["BFILE_PREFIX"], causal_snplist)

    run_command(
        [
            ctx["GCTA"],
            "--bfile",
            ctx["BFILE_PREFIX"],
            "--simu-qt",
            "--simu-causal-loci",
            causal_snplist,
            "--simu-hsq",
            "0.5",
            "--simu-rep",
            "1",
            "--out",
            f"{ctx['RESULTS_DIR']}/sim/sim",
        ]
    )


def cmd_all(ctx: dict[str, str]) -> None:
    cmd_check_tools(ctx)
    cmd_prepare_data(ctx)
    cmd_run_lr(ctx)
    if command_exists(ctx["GEMMA"]):
        cmd_run_lmm(ctx)
    else:
        print("Skipping run-lmm: GEMMA not found.")
    cmd_evaluate(ctx)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Cross-platform pipeline runner.")
    parser.add_argument("--config", default="config/project.env")

    sub = parser.add_subparsers(dest="command", required=True)
    for name in ["check-tools", "prepare-data", "run-lr", "run-lmm", "evaluate", "simulate", "all"]:
        sub.add_parser(name)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    cfg, _ = load_config(args.config)
    ctx = context_from_config(cfg)

    commands = {
        "check-tools": cmd_check_tools,
        "prepare-data": cmd_prepare_data,
        "run-lr": cmd_run_lr,
        "run-lmm": cmd_run_lmm,
        "evaluate": cmd_evaluate,
        "simulate": cmd_simulate,
        "all": cmd_all,
    }
    commands[args.command](ctx)


if __name__ == "__main__":
    main()
