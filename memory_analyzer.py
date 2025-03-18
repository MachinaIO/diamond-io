#!/usr/bin/env python3
"""
Memory Usage Analyzer for Diamond-IO Logs

This script parses a logs.txt file and provides insights about the most 
memory-intensive steps in the circuit obfuscation process.
"""

import re
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
from pathlib import Path


def analyze_memory_usage(log_file):
    """Analyze memory usage from the log file."""
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
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
    
    # Calculate elapsed time in seconds from the first log entry
    start_time = df['timestamp'].iloc[0]
    df['elapsed_seconds'] = (df['timestamp'] - start_time).dt.total_seconds()
    
    return df


def print_memory_insights(df):
    """Print insights about memory usage."""
    print("\n===== MEMORY USAGE INSIGHTS =====\n")
    
    # Overall statistics
    initial_physical = df['physical_memory'].iloc[0]
    final_physical = df['physical_memory'].iloc[-1]
    max_physical = df['physical_memory'].max()
    
    print(f"Initial physical memory: {format_bytes(initial_physical)}")
    print(f"Final physical memory: {format_bytes(final_physical)}")
    print(f"Peak physical memory: {format_bytes(max_physical)}")
    print(f"Total physical memory increase: {format_bytes(final_physical - initial_physical)}")
    
    # Calculate percentage increases if not already calculated
    if 'percentage_increase' not in df.columns:
        df['percentage_increase'] = (df['physical_memory_change'] / df['physical_memory'].shift(1)) * 100
    
    # Print all steps with memory usage details
    print("\n===== ALL STEPS WITH MEMORY USAGE DETAILS =====\n")
    print(f"{'Step Description':<70} {'Memory Usage':<15} {'Absolute Change':<20} {'% Change':<15}")
    print("-" * 120)
    
    for _, row in df.iterrows():
        print(f"{row['message']:<70} {format_bytes(row['physical_memory']):<15} ", end="")
        print(f"{format_bytes(row['physical_memory_change']):<20} ", end="")
        
        if row['physical_memory_change'] != 0:
            print(f"{row['percentage_increase']:.2f}%")
        else:
            print("0.00%")


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


def plot_memory_usage(df, output_file=None):
    """Plot memory usage over time."""
    plt.figure(figsize=(12, 6))
    
    # Plot physical memory usage
    plt.plot(df['elapsed_seconds'], df['physical_memory'] / (1024 ** 2), 
             marker='o', linestyle='-', markersize=4, label='Physical Memory')
    
    # Add labels and title
    plt.xlabel('Time (seconds)')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Memory Usage Over Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # Add annotations for significant memory increases
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
    parser = argparse.ArgumentParser(description='Analyze memory usage from Diamond-IO logs.')
    parser.add_argument('log_file', type=str, help='Path to the logs.txt file')
    parser.add_argument('--plot', action='store_true', help='Generate a memory usage plot')
    parser.add_argument('--output', type=str, help='Output file for the plot (PNG format)')
    
    args = parser.parse_args()
    
    log_file = args.log_file
    if not Path(log_file).exists():
        print(f"Error: Log file '{log_file}' not found.")
        return
    
    df = analyze_memory_usage(log_file)
    if df is not None:
        print_memory_insights(df)
        
        if args.plot:
            plot_memory_usage(df, args.output)


if __name__ == "__main__":
    main()
