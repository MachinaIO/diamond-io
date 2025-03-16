from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from typing import List


# Helper function for parallel processing - must be at module level for pickling
def _process_coefficient(args):
    """Process a single coefficient from its byte representation."""
    i, bytes_data, byte_size, bit_vector_size, bit_vector, modulus = args

    # Calculate start and end indices for this coefficient's bytes
    start = 1 + bit_vector_size + (i * byte_size)
    end = start + byte_size
    value_bytes = bytes_data[start:end]

    # Convert bytes to integer (little-endian)
    value = int.from_bytes(value_bytes, byteorder="little")

    # Check if this coefficient is negative
    byte_idx = i // 8
    bit_idx = i % 8
    is_negative = (bit_vector[byte_idx] & (1 << bit_idx)) != 0

    # Convert back from centered representation
    if is_negative:
        # If negative flag is set, compute modulus - value
        final_value = modulus - value
    else:
        # Otherwise, use value as is
        final_value = value

    return final_value


def from_compact_bytes(
    bytes_data: bytes, ring_dimension: int, modulus: int
) -> List[int]:
    """
    Convert a compact byte representation back to a polynomial.

    Args:
        bytes_data (bytes): The compact byte representation
        ring_dimension (int): The ring dimension
        modulus (int): The modulus

    Returns:
        list: The coefficients of the polynomial
    """
    # First byte contains the byte size per coefficient
    byte_size = bytes_data[0]

    # Next n/8 bytes contain the bit vector
    bit_vector_size = ring_dimension // 8
    bit_vector = bytes_data[1 : 1 + bit_vector_size]

    # Prepare arguments for parallel processing
    args_list = [
        (i, bytes_data, byte_size, bit_vector_size, bit_vector, modulus)
        for i in range(ring_dimension)
    ]

    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        coeffs = list(executor.map(_process_coefficient, args_list))

    return coeffs


import unittest


class TestFromCompactBytes(unittest.TestCase):
    def test_basic_example(self):
        """Test the from_compact_bytes function with a simple example."""
        # Example parameters
        ring_dimension = 8  # Small example for demonstration
        modulus = 257  # Small prime for demonstration

        # Create a sample byte array
        # Format: [byte_size, bit_vector, coefficient bytes...]
        # For this example:
        # - byte_size = 1 (each coefficient takes 1 byte)
        # - bit_vector = [0b00000101] (coefficients 0 and 2 are negative)
        # - coefficients = [10, 20, 30, 40, 50, 60, 70, 80]
        sample_bytes = bytes(
            [
                1,  # byte_size
                0b00000101,  # bit_vector
                10,
                20,
                30,
                40,
                50,
                60,
                70,
                80,  # coefficient values
            ]
        )

        # Decode using the parallel version
        coeffs_parallel = from_compact_bytes(sample_bytes, ring_dimension, modulus)

        # Expected result:
        # For coefficients 0 and 2 (which are marked as negative in the bit vector),
        # we compute modulus - value, so 257 - 10 = 247 and 257 - 30 = 227
        # For the rest, we use the value as is
        expected = [257 - 10, 20, 257 - 30, 40, 50, 60, 70, 80]

        # Assert that the result matches the expected output
        self.assertEqual(coeffs_parallel, expected)

    def test_larger_ring_dimension(self):
        """Test with a larger ring dimension."""
        # Example parameters
        ring_dimension = 16  # Larger ring dimension
        modulus = 257  # Small prime for demonstration

        # Create a sample byte array
        # Format: [byte_size, bit_vector, coefficient bytes...]
        # For this example:
        # - byte_size = 1 (each coefficient takes 1 byte)
        # - bit_vector = [0b00010001, 0b00010001] (coefficients 0, 4, 8, 12 are negative)
        # - coefficients = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160]
        sample_bytes = bytes(
            [
                1,  # byte_size
                0b00010001,
                0b00010001,  # bit_vector
                10,
                20,
                30,
                40,
                50,
                60,
                70,
                80,  # coefficient values (first 8)
                90,
                100,
                110,
                120,
                130,
                140,
                150,
                160,  # coefficient values (next 8)
            ]
        )

        # Decode using the parallel version
        coeffs_parallel = from_compact_bytes(sample_bytes, ring_dimension, modulus)

        # Expected result:
        # For coefficients 0, 4, 8, 12 (which are marked as negative in the bit vector),
        # we compute modulus - value
        # For the rest, we use the value as is
        expected = [
            257 - 10,
            20,
            30,
            40,
            257 - 50,
            60,
            70,
            80,
            257 - 90,
            100,
            110,
            120,
            257 - 130,
            140,
            150,
            160,
        ]

        # Assert that the result matches the expected output
        self.assertEqual(coeffs_parallel, expected)

    def test_different_modulus(self):
        """Test with a different modulus."""
        # Example parameters
        ring_dimension = 8  # Small example for demonstration
        modulus = 1031  # Different prime modulus

        # Create a sample byte array
        # Format: [byte_size, bit_vector, coefficient bytes...]
        # For this example:
        # - byte_size = 2 (each coefficient takes 2 bytes)
        # - bit_vector = [0b00000101] (coefficients 0 and 2 are negative)
        # - coefficients = [300, 400, 500, 600, 700, 800, 900, 1000]
        sample_bytes = bytes(
            [
                2,  # byte_size
                0b00000101,  # bit_vector
                # coefficient values (each 2 bytes, little-endian)
                44,
                1,  # 300 (44 + 1*256)
                144,
                1,  # 400 (144 + 1*256)
                244,
                1,  # 500 (244 + 1*256)
                88,
                2,  # 600 (88 + 2*256)
                188,
                2,  # 700 (188 + 2*256)
                32,
                3,  # 800 (32 + 3*256)
                132,
                3,  # 900 (132 + 3*256)
                232,
                3,  # 1000 (232 + 3*256)
            ]
        )

        # Decode using the parallel version
        coeffs_parallel = from_compact_bytes(sample_bytes, ring_dimension, modulus)

        # Expected result:
        # For coefficients 0 and 2 (which are marked as negative in the bit vector),
        # we compute modulus - value, so 1031 - 300 = 731 and 1031 - 500 = 531
        # For the rest, we use the value as is
        expected = [1031 - 300, 400, 1031 - 500, 600, 700, 800, 900, 1000]

        # Assert that the result matches the expected output
        self.assertEqual(coeffs_parallel, expected)


if __name__ == "__main__":
    unittest.main()
