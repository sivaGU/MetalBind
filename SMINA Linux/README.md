SMINA (Linux) Binary
====================

This directory stores the Linux build of **smina** that MetalloDock uses when
deployed on Streamlit Cloud (Ubuntu).

Contents:

- `smina` – the statically-linked Linux binary (renamed from the upstream
  `smina.static` download).  
- `LICENSE` (optional) – include the license file that accompanies the binary.

Usage notes:

1. The detection logic in `MetalloDock.py` looks for `smina` in this folder as
   well as in `Files_for_GUI/`.  No extra configuration is required.
2. After cloning on Linux, ensure the executable bit is set:

   ```bash
   chmod +x "SMINA Linux/smina"
   ```

3. Commit both the binary and this README to the repository so that Streamlit
   Cloud has access to the executable when the app starts.
SMINA (Linux)
=============

This directory is a placeholder for the Linux build of the **smina** docking
executable used by the MetalloDock hybrid workflow on Streamlit Cloud.

To enable SMINA on Linux:

1. Download the latest `smina.linux` binary from the official repository:
   https://github.com/mwojcikowski/smina/releases
2. Extract the archive and copy the executable into this folder, renaming it to
   simply `smina`.
3. Make the file executable:

   ```bash
   chmod +x "SMINA Linux/smina"
   ```

4. Copy the accompanying SMINA license file (e.g. `LICENSE.GNU`) into this same
   folder so the repository remains compliant.

After committing these files, MetalloDock will automatically detect the binary
at runtime (falling back to `Files_for_GUI/smina` if present).

