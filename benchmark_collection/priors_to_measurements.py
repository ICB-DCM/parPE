#!/usr/bin/env python3
"""For all benchmark collection problems, convert objectivePriors to
measurements."""
import os
import sys
from pathlib import Path
from petab.v1 import (Problem, write_observable_df, write_measurement_df,
                      write_parameter_df)
from petab.v1.priors import priors_to_measurements
from petab.v1.yaml import load_yaml
from petab.v1.C import (MEASUREMENT_FILES, OBSERVABLE_FILES, PARAMETER_FILE,
                        PROBLEMS)


def get_problems_dir() -> Path:
    """Get the directory containing the benchmark collection problems."""
    if len(sys.argv) == 2:
        problems_dir = sys.argv[1]
    elif "BENCHMARK_COLLECTION" in os.environ:
        problems_dir = os.environ["BENCHMARK_COLLECTION"]
    else:
        print("Usage: prior_to_measurements.py <benchmark_collection_dir>")
        sys.exit(1)

    problems_dir = Path(problems_dir)
    if not problems_dir.is_dir():
        print(f"Directory {problems_dir} does not exist.")
        sys.exit(1)

    return problems_dir


def save_problem(yaml_path: Path, petab_problem: Problem):
    """Save the updated PEtab problem."""
    # only meausrements, observables, and parameters are changed
    config = load_yaml(yaml_path)
    assert len(config[PROBLEMS]) == 1
    problem_config = config[PROBLEMS][0]
    assert len(problem_config[MEASUREMENT_FILES]) == 1
    assert len(problem_config[OBSERVABLE_FILES]) == 1
    parameter_file = yaml_path.parent / config[PARAMETER_FILE]
    measurement_file = yaml_path.parent / problem_config[MEASUREMENT_FILES][0]
    observable_file = yaml_path.parent / problem_config[OBSERVABLE_FILES][0]
    # back up the original files
    parameter_file.rename(parameter_file.with_suffix(".tsv.bak"))
    measurement_file.rename(measurement_file.with_suffix(".tsv.bak"))
    observable_file.rename(observable_file.with_suffix(".tsv.bak"))
    # write the updated files
    write_parameter_df(petab_problem.parameter_df, parameter_file)
    write_observable_df(petab_problem.observable_df, observable_file)
    write_measurement_df(petab_problem.measurement_df, measurement_file)


def main():
    problems_dir = get_problems_dir()

    for problem_dir in sorted(problems_dir.iterdir()):
        if not problem_dir.is_dir():
            continue

        problem_id = problem_dir.name
        yaml_path = problem_dir / f"{problem_id}.yaml"

        if not yaml_path.is_file():
            continue
        print(f"{problem_id}...")
        petab_problem = Problem.from_yaml(yaml_path)
        num_measurements_old = len(petab_problem.measurement_df)

        try:
            petab_problem = priors_to_measurements(petab_problem)
            num_measurements_new = len(petab_problem.measurement_df)
            if num_measurements_new != num_measurements_old:
                print(f"\tConverting {problem_id}...")
                save_problem(yaml_path, petab_problem)
            else:
                print(f"\tNothing to do for {problem_id}.")
        except NotImplementedError as e:
            print(f"\t{problem_id}: {e}")


if __name__ == "__main__":
    main()
