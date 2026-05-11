#!/usr/bin/env python3
import requests
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, help="例：ENCSR176KEW、ENCFF896WDT")
args = parser.parse_args()


# Experiment IDから（e.g, ENCSR176KEW）
exp_id = args.input
url = f"https://www.encodeproject.org/experiments/{exp_id}/?format=json"

response = requests.get(url).json()

for rep in response.get('replicates', []):
    acc = rep.get('accession', {})
    lib = rep.get('library', {})
    strand = lib.get('strand_specificity')
    print(f"Exp: {exp_id}, Library: {acc}, Strand: {strand}")
