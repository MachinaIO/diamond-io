#!/usr/bin/env python3
"""
Memory Usage Analyzer for Diamond-IO Logs

This script runs a cargo test command, captures the logs, and provides insights
about the most memory-intensive steps in the circuit obfuscation process.
"""

import re
import pandas as pd
from datetime import datetime
import argparse
from pathlib import Path
import subprocess
import os


def run_cargo_test():
    """Run the cargo test command and capture the output."""
    command = [
        "cargo", "test", "-r", "--package", "diamond-io", "--lib", 
        "--features=parallel", "--", "io::test::test_io_just_mul_enc_and_bit", 
        "--exact", "--show-output"
    ]
    
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
    
    # Parse log lines in pairs
    log_entries = []
    
    # The log file has pairs of lines:
    # 1. A line with timestamp, module, and step description
    # 2. A line with memory usage information
    
    i = 0
    while i < len(lines) - 1:
        # First line should contain timestamp, module, and message
        line1 = lines[i].strip()
        line2 = lines[i+1].strip()
        
        # Extract step description from the first line
        step_pattern = r'(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+Z)\s+INFO\s+([\w:]+):\s+(.*)'
        step_match = re.search(step_pattern, line1)
        
        # Extract memory usage from the second line
        memory_pattern = r'Current physical/virtural memory usage: (\d+) / (\d+)'
        memory_match = re.search(memory_pattern, line2)
        
        if step_match and memory_match:
            timestamp_str = step_match.group(1)
            module = step_match.group(2)
            message = step_match.group(3).strip()
            physical_memory = int(memory_match.group(1))
            virtual_memory = int(memory_match.group(2))
            
            # Convert timestamp to datetime object
            timestamp = datetime.strptime(timestamp_str, '%Y-%m-%dT%H:%M:%S.%fZ')
            
            log_entries.append({
                'timestamp': timestamp,
                'module': module,
                'message': message,
                'physical_memory': physical_memory,
                'virtual_memory': virtual_memory
            })
            
            # Move to the next pair of lines
            i += 2
        else:
            # If the pattern doesn't match, move to the next line
            i += 1
    
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
    table += f"{'Step Description':<70} {'Physical Memory':<15} {'Virtual Memory':<15} {'Physical Change':<15} {'Virtual Change':<15}\n"
    table += "-" * 130 + "\n"
    
    for _, row in df.iterrows():
        table += f"{row['message']:<70} "
        table += f"{format_bytes(row['physical_memory']):<15} "
        table += f"{format_bytes(row['virtual_memory']):<15} "
        table += f"{format_bytes(row['physical_memory_change']):<15} "
        table += f"{format_bytes(row['virtual_memory_change']):<15}\n"
    
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


def plot_memory_usage(df, output_file=None):
    """Plot memory usage over time."""
    plt.figure(figsize=(12, 6))
    
    # Plot physical memory usage
    plt.plot(df['elapsed_seconds'], df['physical_memory'] / (1024 ** 2), 
             marker='o', linestyle='-', markersize=4, label='Physical Memory (MB)')
    
    # Plot virtual memory usage on a secondary y-axis
    ax2 = plt.gca().twinx()
    ax2.plot(df['elapsed_seconds'], df['virtual_memory'] / (1024 ** 3), 
             marker='s', linestyle='-', color='red', markersize=4, label='Virtual Memory (GB)')
    
    # Add labels and title
    plt.xlabel('Time (seconds)')
    plt.gca().set_ylabel('Physical Memory (MB)')
    ax2.set_ylabel('Virtual Memory (GB)')
    plt.title('Memory Usage Over Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Combine legends
    lines1, labels1 = plt.gca().get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    
    # Add annotations for significant physical memory increases
    significant_changes = df[df['physical_memory_change'] > df['physical_memory_change'].quantile(0.9)]
    for _, row in significant_changes.iterrows():
        plt.annotate(
            row['message'][:20] + '...' if len(row['message']) > 20 else row['message'],
            xy=(row['elapsed_seconds'], row['physical_memory'] / (1024 ** 2)),
            xytext=(10, 10),
            textcoords='offset points',
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=.2')
        )
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def main():
    # Create logs directory if it doesn't exist
    logs_dir = Path("logs")
    logs_dir.mkdir(exist_ok=True)
    
    # Get current timestamp for file naming
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Run the cargo test command and capture the output
    log_output = run_cargo_test()
    
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
        
        # Save the tables to files
        physical_table_path = logs_dir / f"physical_memory_analysis_{timestamp}.txt"
        with open(physical_table_path, 'w') as f:
            f.write(physical_table)
        print(f"Physical memory analysis saved to {physical_table_path}")
        
        virtual_table_path = logs_dir / f"virtual_memory_analysis_{timestamp}.txt"
        with open(virtual_table_path, 'w') as f:
            f.write(virtual_table)
        print(f"Virtual memory analysis saved to {virtual_table_path}")


if __name__ == "__main__":
    main()
