"""
Script to find taxa with "Methano" in their names from GTDB taxonomy table "gtdb-taxonomy-table.tsv" downloaded from https://gtdb.ecogenomic.org/browsers.
Collects taxa at each level (Phylum, Class, Order, Family) that start with "Methano"
Outputs the results to a text file and prints them to the console.

$ python find_Methano_in_GTDB.py
"""

import csv
from pathlib import Path


def collect_methano_taxa(input_path: Path):
	phylum_tuple = set()
	class_tuple = set()
	order_tuple = set()
	family_tuple = set()

	with input_path.open("r", encoding="utf-8", newline="") as f:
		reader = csv.DictReader(f, delimiter="\t")

		required_cols = {"Domain", "Phylum", "Class", "Order", "Family"}
		missing = required_cols - set(reader.fieldnames or [])
		if missing:
			raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

		for row in reader:
			domain = (row.get("Domain") or "").strip()
			phylum = (row.get("Phylum") or "").strip()
			class_ = (row.get("Class") or "").strip()
			order = (row.get("Order") or "").strip()
			family = (row.get("Family") or "").strip()

			if phylum.startswith("Methano"):
				phylum_tuple.add(f"d__{domain}; p__{phylum}")

			if class_.startswith("Methano"):
				class_tuple.add(f"d__{domain}; p__{phylum}; c__{class_}")

			if order.startswith("Methano"):
				order_tuple.add(f"d__{domain}; p__{phylum}; c__{class_}; o__{order}")

			if family.startswith("Methano"):
				family_tuple.add(f"d__{domain}; p__{phylum}; c__{class_}; o__{order}; f__{family}")

	return tuple(sorted(phylum_tuple)), tuple(sorted(class_tuple)), tuple(sorted(order_tuple)), tuple(sorted(family_tuple))


def print_tuple_items(items: tuple[str, ...], title: str):
	print(title)
	for item in items:
		print(item)
	print()


def write_tuple_items(output_path: Path, groups: list[tuple[str, tuple[str, ...]]]):
	with output_path.open("w", encoding="utf-8", newline="") as f:
		for title, items in groups:
			f.write(f"{title}\n")
			for item in items:
				f.write(f"{item}\n")
			f.write("\n")


def main():
	input_file = Path("gtdb-taxonomy-table.tsv")
	output_file = Path("methano_taxonomy_output.txt")

	if not input_file.exists():
		raise FileNotFoundError(f"Input file not found: {input_file}")

	phylum_items, class_items, order_items, family_items = collect_methano_taxa(input_file)

	groups = [
		("Phylum tuple:", phylum_items),
		("Class tuple:", class_items),
		("Order tuple:", order_items),
		("Family tuple:", family_items),
	]

	write_tuple_items(output_file, groups)

	for title, items in groups:
		print_tuple_items(items, title)


if __name__ == "__main__":
	main()
