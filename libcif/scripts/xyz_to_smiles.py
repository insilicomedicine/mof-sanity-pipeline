#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import click
import sys
import json
from tqdm import tqdm


def convert_xyz_folder(folder_path: Path, max_time: int):
    log_path = folder_path / "smi_conversion.log"
    json_log_path = folder_path / "smi_conversion.json"
    results = []
    
    # Parse building blocks information if available
    building_blocks_info = {}
    num_building_blocks_file = folder_path / "num_building_blocks.txt"
    if num_building_blocks_file.exists():
        try:
            with num_building_blocks_file.open("r") as f:
                for line in f:
                    if ': ' in line:
                        filename, count = line.strip().split(': ')
                        building_blocks_info[filename] = int(count)
        except Exception as e:
            print(f"Error reading building blocks file: {e}")

    with log_path.open("w") as log_file, contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
        print(f"Processing folder: {folder_path}")
        for xyz_path in folder_path.glob("*.xyz"):
            print(f"Converting: {xyz_path.name}", flush=True)

            entry = {"file": xyz_path.name}
            
            # Add building block information for this specific file
            if xyz_path.name in building_blocks_info:
                entry["building_block_count"] = building_blocks_info[xyz_path.name]
                
            # Determine building block type
            if xyz_path.name.startswith('linker'):
                entry["building_block_type"] = "linker"
            elif xyz_path.name.startswith('node'):
                entry["building_block_type"] = "node"
            elif xyz_path.name.startswith('disconnected'):
                entry["building_block_type"] = "disconnected"
            else:
                entry["building_block_type"] = "unknown"

            if "infinite" in xyz_path.name.lower():
                print(f"Skipping due to 'infinite' in name: {xyz_path.name}")
                entry["status"] = "infinite"
                results.append(entry)
                continue

            try:
                with xyz_path.open("r") as f:
                    lines = f.readlines()

                n_atoms = int(lines[0].strip())
                if len(lines) < n_atoms + 2:
                    msg = "File too short"
                    print(f"{msg}: {xyz_path}")
                    entry.update({"status": "failure", "error": msg})
                    results.append(entry)
                    continue

                cleaned_lines = [str(n_atoms) + "\n", xyz_path.stem + "\n"]
                cleaned_lines.extend(lines[2:2 + n_atoms])
            except Exception as e:
                msg = f"Error processing: {e}"
                print(msg)
                entry.update({"status": "failure", "error": str(e)})
                results.append(entry)
                continue

            tmp_xyz_path = xyz_path.with_suffix(".tmp.xyz")
            try:
                tmp_xyz_path.write_text("".join(cleaned_lines))

                obabel_variants = [
                    {
                        "name": xyz_path.stem + "_plain.smi",
                        "cmd": ["obabel", "-ixyz", str(tmp_xyz_path), "-osmi", "-O", str(folder_path / (xyz_path.stem + "_plain.smi"))]
                    }
                ]

                if "linker" in xyz_path.name.lower() or "disconnected" in xyz_path.name.lower():
                    obabel_variants.append({
                        "name": xyz_path.stem + "_ph7.smi",
                        "cmd": ["obabel", "-ixyz", str(tmp_xyz_path), "-osmi", "-O", str(folder_path / (xyz_path.stem + "_ph7.smi")), "--addh", "--p", "7.4"]
                    })

                entry["variants"] = []
                for variant in obabel_variants:
                    result = subprocess.run(
                        variant["cmd"],
                        capture_output=True, text=True, timeout=max_time
                    )

                    log_file.write(result.stdout)
                    log_file.write(result.stderr)

                    smi_path_variant = folder_path / variant["name"]
                    variant_entry = {
                        "variant": variant["name"]
                    }

                    if result.returncode != 0 or "1 molecule converted" not in result.stderr:
                        msg = f"Conversion failed: {variant['name']}"
                        print(msg)
                        variant_entry.update({
                            "status": "failure",
                            "error": result.stderr.strip()
                        })
                    else:
                        print(f"Success: {variant['name']}", flush=True)
                        try:
                            smi_content = " ".join(smi_path_variant.read_text().strip().split()[:-1])
                        except Exception as e:
                            smi_content = f"<failed to read .smi: {e}>"
                        variant_entry.update({
                            "status": "success",
                            "output_smi": smi_content
                        })

                    entry["variants"].append(variant_entry)

            except subprocess.TimeoutExpired:
                msg = f"Open Babel timeout after {max_time}s"
                print(msg)
                entry.update({"status": "failure", "error": msg})
            except Exception as e:
                msg = f"Open Babel error: {e}"
                print(msg)
                entry.update({"status": "failure", "error": str(e)})
            finally:
                tmp_xyz_path.unlink(missing_ok=True)
                results.append(entry)

    with json_log_path.open("w") as jf:
        json.dump(results, jf, indent=2)


def process_subdir(args):
    subdir, max_time = args
    if subdir.is_dir():
        convert_xyz_folder(subdir, max_time)


@click.command()
@click.argument("metafolder", type=click.Path(exists=True, file_okay=False))
@click.option("--max-workers", default=64, show_default=True, help="Number of parallel workers")
@click.option("--max-time", default=10, show_default=True, help="Maximum time (s) to allow for obabel execution")
@click.option("--output-json", default=None, help="Filename for combined JSON summary (default: <metafolder>.json)")
def main(metafolder, max_workers, max_time, output_json):
    """
    Convert all .xyz files inside each subfolder of METAFOLDER to .smi using Open Babel.
    Logs are saved to smi_conversion.log and smi_conversion.json in each subfolder.
    After processing, a combined JSON summary is written.
    """
    root = Path(metafolder).resolve()
    subdirs = [d for d in root.iterdir() if d.is_dir()]
    if output_json is None:
        output_json = f"{root.name}.json"
    output_json = Path(output_json).resolve()

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_subdir, (d, max_time)): d for d in subdirs
        }
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing folders"):
            d = futures[future]
            try:
                future.result()
            except Exception as exc:
                print(f"Exception in {d}: {exc}", file=sys.stderr)

    # Combine all smi_conversion.json files
    combined = {}
    for subdir in tqdm(subdirs):
        json_path = subdir / "smi_conversion.json"
        if json_path.exists():
            try:
                with json_path.open("r") as jf:
                    data = json.load(jf)
                    combined[subdir.name] = data
            except Exception as e:
                print(f"Error reading {json_path}: {e}", file=sys.stderr)

    with output_json.open("w") as outf:
        json.dump(combined, outf, indent=2)

    print(f"\n✅ Combined summary written to: {output_json}")


if __name__ == "__main__":
    main()
