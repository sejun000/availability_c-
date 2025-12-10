#!/usr/bin/env python3
"""
Extract SSD failure times from Alibaba trace CSV and output as simple list.
Output format: one failure time per line (in hours)
"""

import csv
import sys
import argparse

def extract_failure_times(input_csv, output_file, column_name='S1', time_unit='days'):
    """
    Extract failure times from CSV and write to output file.

    Args:
        input_csv: Path to input CSV file
        output_file: Path to output file
        column_name: Column name containing failure time (default: S1)
        time_unit: Unit of time in CSV ('days' or 'hours')
    """
    failure_times = []

    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)

        for row in reader:
            value = row.get(column_name, '').strip()

            # Skip NA or empty values
            if value == 'NA' or value == '' or value is None:
                continue

            try:
                time_val = float(value)

                # Convert to hours if needed
                if time_unit == 'days':
                    time_val *= 24.0

                if time_val >= 0:
                    failure_times.append(time_val)
            except ValueError:
                continue

    # Sort failure times
    failure_times.sort()

    # Write to output file
    with open(output_file, 'w') as f:
        for t in failure_times:
            f.write(f"{t:.2f}\n")

    # Print statistics
    if failure_times:
        print(f"Extracted {len(failure_times)} failure times")
        print(f"Min: {min(failure_times):.2f} hours ({min(failure_times)/24:.1f} days)")
        print(f"Max: {max(failure_times):.2f} hours ({max(failure_times)/24:.1f} days)")
        print(f"Mean: {sum(failure_times)/len(failure_times):.2f} hours")
        print(f"Median: {failure_times[len(failure_times)//2]:.2f} hours")
        print(f"Output written to: {output_file}")
    else:
        print("No valid failure times found!")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Extract SSD failure times from CSV')
    parser.add_argument('input', help='Input CSV file')
    parser.add_argument('-o', '--output', default='ssd_failure_times.txt',
                        help='Output file (default: ssd_failure_times.txt)')
    parser.add_argument('-c', '--column', default='S1',
                        help='Column name for failure time (default: S1)')
    parser.add_argument('-u', '--unit', default='days', choices=['days', 'hours'],
                        help='Time unit in CSV (default: days)')

    args = parser.parse_args()

    extract_failure_times(args.input, args.output, args.column, args.unit)

if __name__ == '__main__':
    main()
