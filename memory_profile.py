#!/usr/bin/env python3
"""
Memory Usage Analyzer for Diamond-IO Logs

This script runs a cargo test command, captures the logs, and provides insights
about the most memory-intensive steps in the circuit obfuscation process.
"""

import re
import pandas as pd
from datetime import datetime
from pathlib import Path
import subprocess
import sys
import argparse


def run_cargo_test(command):
    """Run a command and capture the output, ensuring logs are saved even if interrupted."""
    print(f"Running command: {' '.join(command)}")
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True
    )
    
    # Capture output line by line
    all_output = []
    try:
        for line in iter(process.stdout.readline, ''):
            print(line, end='')  # Print in real-time
            all_output.append(line)
            
            # Stop when we see the completion message
            if "OBFUSCATION COMPLETED" in line:
                print("Obfuscation completed, stopping log capture.")
                break
    except KeyboardInterrupt:
        print("\nProcess interrupted! Saving captured logs...\n")
    finally:
        # Ensure process is terminated properly
        process.stdout.close()
        process.terminate()
        process.wait()

    return ''.join(all_output)



def strip_ansi_codes(text):
    """Remove ANSI color codes from text."""
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return ansi_escape.sub('', text)


def analyze_memory_usage_from_string(log_output):
    """
    Parses memory usage logs and returns a DataFrame.
    Skips invalid entries instead of crashing.
    """
    # First strip all ANSI codes from the output
    log_output = strip_ansi_codes(log_output)
    
    # Pattern for lines with message
    pattern_with_message = re.compile(
        r".*?diamond_io::utils\s*:\s*(.*?)\s*\|\|\s*Current physical/virtual memory usage:\s*(\d+)\s*\|\s*(\d+)"
    )
    
    # Pattern for lines with just numbers
    pattern_numbers = re.compile(
        r"(\d+)\s*\|\s*(\d+)"
    )

    data = []
    
    for line in log_output.splitlines():
        # Try to match line with message first
        match = pattern_with_message.search(line)
        if match:
            try:
                message = match.group(1).strip()
                physical_memory = int(match.group(2))
                virtual_memory = int(match.group(3))
                data.append({
                    'message': message,
                    'physical_memory': physical_memory,
                    'virtual_memory': virtual_memory
                })
            except ValueError as e:
                print(f"Skipping invalid log entry: {line} (Error: {e})")
        else:
            # Try to match line with just numbers
            match = pattern_numbers.search(line)
            if match:
                try:
                    physical_memory = int(match.group(1))
                    virtual_memory = int(match.group(2))
                    # Use a generic message for number-only lines
                    data.append({
                        'message': f"Memory measurement {len(data) + 1}",
                        'physical_memory': physical_memory,
                        'virtual_memory': virtual_memory
                    })
                except ValueError as e:
                    print(f"Skipping invalid log entry: {line} (Error: {e})")
    
    if not data:
        return None  # No valid logs

    df = pd.DataFrame(data)
    
    # Calculate changes and percentages
    df['physical_memory_change'] = df['physical_memory'].diff().fillna(0)
    df['virtual_memory_change'] = df['virtual_memory'].diff().fillna(0)
    
    # Calculate percentage changes
    df['physical_percentage_increase'] = (df['physical_memory_change'] / df['physical_memory'].shift(1) * 100).fillna(0)
    df['virtual_percentage_increase'] = (df['virtual_memory_change'] / df['virtual_memory'].shift(1) * 100).fillna(0)
    
    return df


def format_bytes(bytes_value):
    """Format bytes to human-readable format."""
    if bytes_value < 1024:
        return f"{bytes_value} B"
    elif bytes_value < 1024 ** 2:
        return f"{bytes_value / 1024:.2f} KB"
    elif bytes_value < 1024 ** 3:
        return f"{bytes_value / (1024 ** 2):.2f} MB"
    else:
        return f"{bytes_value / (1024 ** 3):.2f} GB"


def generate_physical_memory_table(df):
    """Generate a table for physical memory usage."""
    table = "===== PHYSICAL MEMORY USAGE =====\n\n"
    table += f"{'Step Description':<70} {'Memory Usage':<15} {'Absolute Change':<20} {'% Change':<15}\n"
    table += "-" * 120 + "\n"
    
    for _, row in df.iterrows():
        table += f"{row['message']:<70} {format_bytes(row['physical_memory']):<15} "
        table += f"{format_bytes(row['physical_memory_change']):<20} "
        
        if row['physical_memory_change'] != 0:
            table += f"{row['physical_percentage_increase']:.2f}%\n"
        else:
            table += "0.00%\n"
    
    # Add summary statistics
    initial_physical = df['physical_memory'].iloc[0]
    final_physical = df['physical_memory'].iloc[-1]
    max_physical = df['physical_memory'].max()
    
    table += "\n===== SUMMARY =====\n"
    table += f"Initial physical memory: {format_bytes(initial_physical)}\n"
    table += f"Final physical memory: {format_bytes(final_physical)}\n"
    table += f"Peak physical memory: {format_bytes(max_physical)}\n"
    table += f"Total physical memory increase: {format_bytes(final_physical - initial_physical)}\n"
    
    return table


def generate_virtual_memory_table(df):
    """Generate a table for virtual memory usage."""
    table = "===== VIRTUAL MEMORY USAGE =====\n\n"
    table += f"{'Step Description':<70} {'Memory Usage':<15} {'Absolute Change':<20} {'% Change':<15}\n"
    table += "-" * 120 + "\n"
    
    for _, row in df.iterrows():
        table += f"{row['message']:<70} {format_bytes(row['virtual_memory']):<15} "
        table += f"{format_bytes(row['virtual_memory_change']):<20} "
        
        if row['virtual_memory_change'] != 0:
            table += f"{row['virtual_percentage_increase']:.2f}%\n"
        else:
            table += "0.00%\n"
    
    # Add summary statistics
    initial_virtual = df['virtual_memory'].iloc[0]
    final_virtual = df['virtual_memory'].iloc[-1]
    max_virtual = df['virtual_memory'].max()
    
    table += "\n===== SUMMARY =====\n"
    table += f"Initial virtual memory: {format_bytes(initial_virtual)}\n"
    table += f"Final virtual memory: {format_bytes(final_virtual)}\n"
    table += f"Peak virtual memory: {format_bytes(max_virtual)}\n"
    table += f"Total virtual memory increase: {format_bytes(final_virtual - initial_virtual)}\n"
    
    return table


def generate_combined_memory_table(df):
    """Generate a table for combined physical and virtual memory usage."""
    table = "===== COMBINED MEMORY USAGE =====\n\n"
    table += f"{'Step Description':<70} {'Physical Memory':<15} {'Virtual Memory':<15} {'Physical Change':<15} {'% Change':<10} {'Virtual Change':<15} {'% Change':<10}\n"
    table += "-" * 150 + "\n"
    
    for _, row in df.iterrows():
        table += f"{row['message']:<70} "
        table += f"{format_bytes(row['physical_memory']):<15} "
        table += f"{format_bytes(row['virtual_memory']):<15} "
        table += f"{format_bytes(row['physical_memory_change']):<15} "
        
        # Add physical memory percentage change
        if row['physical_memory_change'] != 0:
            table += f"{row['physical_percentage_increase']:.2f}%{'':<5} "
        else:
            table += f"0.00%{'':<5} "
        
        # Add virtual memory change and percentage
        table += f"{format_bytes(row['virtual_memory_change']):<15} "
        
        if row['virtual_memory_change'] != 0:
            table += f"{row['virtual_percentage_increase']:.2f}%\n"
        else:
            table += f"0.00%\n"
    
    # Add summary statistics
    initial_physical = df['physical_memory'].iloc[0]
    final_physical = df['physical_memory'].iloc[-1]
    max_physical = df['physical_memory'].max()
    
    initial_virtual = df['virtual_memory'].iloc[0]
    final_virtual = df['virtual_memory'].iloc[-1]
    max_virtual = df['virtual_memory'].max()
    
    table += "\n===== SUMMARY =====\n"
    table += f"Initial memory (physical/virtual): {format_bytes(initial_physical)} / {format_bytes(initial_virtual)}\n"
    table += f"Final memory (physical/virtual): {format_bytes(final_physical)} / {format_bytes(final_virtual)}\n"
    table += f"Peak memory (physical/virtual): {format_bytes(max_physical)} / {format_bytes(max_virtual)}\n"
    table += f"Total increase (physical/virtual): {format_bytes(final_physical - initial_physical)} / {format_bytes(final_virtual - initial_virtual)}\n"
    
    return table

def analyze_log_file(log_file_path):
    """Analyze memory usage from an existing log file."""
    try:
        with open(log_file_path, 'r') as f:
            log_content = f.read()
        return analyze_memory_usage_from_string(log_content)
    except Exception as e:
        print(f"Error reading or analyzing log file: {e}")
        return None

def save_analysis_files(df, logs_dir, timestamp):
    """Save analysis files if we have valid data."""
    if df is None or df.empty:
        print("No valid log entries found. Skipping file save.")
        return False

    physical_table = generate_physical_memory_table(df)
    virtual_table = generate_virtual_memory_table(df)
    combined_table = generate_combined_memory_table(df)
    
    print("\n" + combined_table)
    
    physical_table_path = logs_dir / f"physical_memory_analysis_{timestamp}.txt"
    virtual_table_path = logs_dir / f"virtual_memory_analysis_{timestamp}.txt"
    combined_table_path = logs_dir / f"combined_memory_analysis_{timestamp}.txt"
    
    with open(physical_table_path, 'w') as f:
        f.write(physical_table)
    with open(virtual_table_path, 'w') as f:
        f.write(virtual_table)
    with open(combined_table_path, 'w') as f:
        f.write(combined_table)

    print(f"Physical memory analysis saved to {physical_table_path}")
    print(f"Virtual memory analysis saved to {virtual_table_path}")
    print(f"Combined memory analysis saved to {combined_table_path}")
    return True

def main():
    """Run the memory profiler with the specified command or analyze an existing log file."""
    parser = argparse.ArgumentParser(description='Memory Usage Analyzer for Diamond-IO Logs')
    parser.add_argument('--log-file', help='Path to an existing log file to analyze')
    parser.add_argument('command', nargs=argparse.REMAINDER, help='Command to run (if not using --log-file)')
    args = parser.parse_args()

    logs_dir = Path("logs")
    logs_dir.mkdir(exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if args.log_file:
        # Analyze existing log file
        df = analyze_log_file(args.log_file)
        if df is not None and not df.empty:
            save_analysis_files(df, logs_dir, timestamp)
        return

    if not args.command:
        print("Error: No command specified.")
        print("Usage: python memory_profiler.py <command> [args...]")
        print("       python memory_profiler.py --log-file <path>")
        sys.exit(1)

    print("Running memory profiler...")
    log_output = ""
    df = None

    try:
        log_output = run_cargo_test(args.command)
        df = analyze_memory_usage_from_string(log_output)
    except KeyboardInterrupt:
        print("\nProcess interrupted! Checking logs before saving...\n")
        if log_output:  # Only analyze if we have some output
            df = analyze_memory_usage_from_string(log_output)
    except Exception as e:
        print(f"Error during execution: {e}")
        if log_output:  # Try to analyze any output we got
            df = analyze_memory_usage_from_string(log_output)

    if df is not None and not df.empty:
        # Save raw logs only if we have valid data
        log_file_path = logs_dir / f"test_logs_{timestamp}.txt"
        with open(log_file_path, 'w') as f:
            f.write(log_output)
        print(f"Raw logs saved to {log_file_path}")
        
        # Save analysis files
        save_analysis_files(df, logs_dir, timestamp)
    else:
        print("No valid log entries found. No files will be written.")

if __name__ == "__main__":
    main()
