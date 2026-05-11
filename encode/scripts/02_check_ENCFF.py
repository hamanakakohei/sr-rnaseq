#!/usr/bin/env python3
import requests
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, help="例：ENCSR176KEW、ENCFF896WDT")
args = parser.parse_args()


# File IDから（e.g., ENCFF896WDT）
file_id = args.input
url = f"https://www.encodeproject.org/files/{file_id}/?format=json"

response = requests.get(url).json()

ass = response.get('assembly')
type = response.get('output_type')
print(f"Exp: {file_id}, Assembly: {ass}, Type: {type}")
