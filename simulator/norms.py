import json
from typing import List
import os


class CircuitNorms:
    """
    A class to read and parse circuit norm data from JSON files.

    This class reads JSON files in the format of final_bits_norm.json, which contains
    a list of lists of string values representing integers, and converts them to
    a list of lists of integers.
    """

    def __init__(self, h_norms: List[List[int]], log_base_q: int):
        """
        Initialize a CircuitNorms object.

        Args:
            h_norms: List of lists of integers representing the norms.
        """
        self.h_norms = h_norms
        self.log_base_q = log_base_q

    @classmethod
    def from_json(cls, json_data: dict, log_base_q: int) -> "CircuitNorms":
        """
        Create a CircuitNorms object from JSON data.

        Args:
            json_data: A dictionary containing the JSON data with an 'h_norms' key
                      that maps to a list of lists of string values.

        Returns:
            A CircuitNorms object with the parsed data.
        """
        if "h_norms" not in json_data:
            raise ValueError("JSON data must contain an 'h_norms' key")

        # Convert string values to integers
        h_norms = []
        for norm_list in json_data["h_norms"]:
            h_norms.append([int(value) for value in norm_list])

        return cls(h_norms, log_base_q)

    @classmethod
    def load_from_file(cls, file_path: str, log_base_q: int) -> "CircuitNorms":
        """
        Load a CircuitNorms object from a JSON file.

        Args:
            file_path: Path to the JSON file containing the norm data.

        Returns:
            A CircuitNorms object with the parsed data.
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        with open(file_path, "r") as f:
            json_data = json.load(f)

        return cls.from_json(json_data, log_base_q)

    def get_h_norms(self) -> List[List[int]]:
        """
        Get the h_norms data.

        Returns:
            A list of lists of integers representing the norms.
        """
        return self.h_norms

    def compute_norms(self, m: int, n: int, base: int) -> List[int]:
        bit_norms = [0 for _ in range(len(self.h_norms))]
        max_deg = max([len(norm) for norm in self.h_norms])
        power_ms = [m**i for i in range(max_deg)]
        for coeffs in self.h_norms:
            for i, coeff in enumerate(coeffs):
                bit_norms[i] += coeff * power_ms[i]
        assert len(bit_norms) % self.log_base_q == 0
        norms = []
        for i in range(0, len(bit_norms), self.log_base_q):
            sum = 0
            for j in range(self.log_base_q):
                sum += bit_norms[self.log_base_q * i + j] * (base - 1) * n * m
            norms.append(sum)
        return norms

    def __len__(self) -> int:
        """
        Get the number of norm lists.

        Returns:
            The number of norm lists.
        """
        return len(self.h_norms)

    def __getitem__(self, index: int) -> List[int]:
        """
        Get a specific norm list by index.

        Args:
            index: The index of the norm list to get.

        Returns:
            A list of integers representing the norms at the specified index.
        """
        return self.h_norms[index]
