#!/usr/bin/env python3
"""
Setup script to ensure Linux executables have execute permissions.
This runs automatically when the app starts to fix permissions on Streamlit Cloud.
"""
import os
import stat
from pathlib import Path

def setup_executables():
    """Set execute permissions on Linux executables."""
    # Find the Files_for_GUI directory
    script_dir = Path(__file__).parent
    files_gui_dir = script_dir / "Files_for_GUI"
    
    if not files_gui_dir.exists():
        return
    
    # List of executables that need execute permissions
    executables = ["vina", "autogrid4", "autodock4"]
    
    for exe_name in executables:
        exe_path = files_gui_dir / exe_name
        if exe_path.exists() and not exe_path.is_dir():
            try:
                # Get current permissions
                current_mode = exe_path.stat().st_mode
                # Add execute permissions for user, group, and others
                new_mode = current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
                os.chmod(exe_path, new_mode)
                print(f"✓ Set execute permissions on {exe_path}")
            except Exception as e:
                print(f"⚠ Could not set permissions on {exe_path}: {e}")

if __name__ == "__main__":
    setup_executables()

