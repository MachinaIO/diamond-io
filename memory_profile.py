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


def run_cargo_test(command):
    """Run a command and capture the output.
    
    Args:
        command: List of command arguments to execute.
    """
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
    for line in iter(process.stdout.readline, ''):
        print(line, end='')  # Print in real-time
        all_output.append(line)
        
        # Stop when we see the completion message
        if "OBFUSCATION COMPLETED" in line:
            print("Obfuscation completed, stopping log capture.")
            break
    
    # Wait for the process to complete
    process.stdout.close()
    process.wait()
    
    return ''.join(all_output)


def strip_ansi_codes(text):
    """Remove ANSI color codes from text."""
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return ansi_escape.sub('', text)


def analyze_memory_usage_from_string(log_string):
    """Analyze memory usage from a log string."""
    # Strip ANSI color codes from the log string
    clean_log = strip_ansi_codes(log_string)
    lines = clean_log.splitlines()
    
    # Parse log lines
    log_entries = []
    
    # New pattern to match the actual log format where memory info is on the same line
    # Example: 2025-03-25T08:35:31.295990Z  INFO diamond_io::utils: Sampled public data || Current physical/virtural memory usage: 7651328 | 420061642752
    log_pattern = r'(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+Z)\s+INFO\s+([\w:]+):\s+(.*?)\s+\|\|\s+Current physical/virtural memory usage:\s+(\d+)\s+\|\s+(\d+)'
    
    for line in lines:
        line = line.strip()
        match = re.search(log_pattern, line)
        
        if match:
            timestamp_str = match.group(1)
            module = match.group(2)
            message = match.group(3).strip()
            physical_memory = int(match.group(4))
            virtual_memory = int(match.group(5))
            
            # Convert timestamp to datetime object
            timestamp = datetime.strptime(timestamp_str, '%Y-%m-%dT%H:%M:%S.%fZ')
            
            log_entries.append({
                'timestamp': timestamp,
                'module': module,
                'message': message,
                'physical_memory': physical_memory,
                'virtual_memory': virtual_memory
            })
    
    if not log_entries:
        print("No valid log entries found.")
        return None
    
    # Convert to DataFrame for easier analysis
    df = pd.DataFrame(log_entries)
    
    # Calculate memory changes between steps
    df['physical_memory_change'] = df['physical_memory'].diff().fillna(0)
    df['virtual_memory_change'] = df['virtual_memory'].diff().fillna(0)
    
    # Calculate percentage increases
    df['physical_percentage_increase'] = (df['physical_memory_change'] / df['physical_memory'].shift(1)) * 100
    df['virtual_percentage_increase'] = (df['virtual_memory_change'] / df['virtual_memory'].shift(1)) * 100
    
    # Calculate elapsed time in seconds from the first log entry
    start_time = df['timestamp'].iloc[0]
    df['elapsed_seconds'] = (df['timestamp'] - start_time).dt.total_seconds()
    
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

def main():
    """Run the memory profiler with the specified command."""
    # Create logs directory if it doesn't exist
    logs_dir = Path("logs")
    logs_dir.mkdir(exist_ok=True)
    
    # Get current timestamp for file naming
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Get command from command line arguments
    if len(sys.argv) > 1:
        command = sys.argv[1:]
    else:
        print("Error: No command specified.")
        print("Usage: python memory_profiler.py <command> [args...]")
        print("Example: python memory_profiler.py cargo test -r --test test_io_dummy_param --no-default-features -- --nocapture")
        sys.exit(1)
    
    # Run the command and capture the output
    print("Running memory profiler...")
    log_output = run_cargo_test(command)
    
    # Save the raw logs to a file
    log_file_path = logs_dir / f"test_logs_{timestamp}.txt"
    with open(log_file_path, 'w') as f:
        f.write(log_output)
    print(f"Raw logs saved to {log_file_path}")
    
    # Analyze the logs
    df = analyze_memory_usage_from_string(log_output)
    
    if df is not None:
        # Generate the tables
        physical_table = generate_physical_memory_table(df)
        virtual_table = generate_virtual_memory_table(df)
        combined_table = generate_combined_memory_table(df)
        
        # Output to console
        print("\n" + combined_table)
        
        # Save to files
        physical_table_path = logs_dir / f"physical_memory_analysis_{timestamp}.txt"
        with open(physical_table_path, 'w') as f:
            f.write(physical_table)
        print(f"Physical memory analysis saved to {physical_table_path}")
        
        virtual_table_path = logs_dir / f"virtual_memory_analysis_{timestamp}.txt"
        with open(virtual_table_path, 'w') as f:
            f.write(virtual_table)
        print(f"Virtual memory analysis saved to {virtual_table_path}")
        
        combined_table_path = logs_dir / f"combined_memory_analysis_{timestamp}.txt"
        with open(combined_table_path, 'w') as f:
            f.write(combined_table)
        print(f"Combined memory analysis saved to {combined_table_path}")


if __name__ == "__main__":
    main()
