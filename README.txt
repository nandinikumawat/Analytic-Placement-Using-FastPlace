
# README

## Instructions

1. Place the source (`.cpp`), header (`.h`), and Python files (`plot_placement.py`) in the directory for each benchmark.

2. Run `make` to compile the program.

3. Execute the parser using the command:
   ```
   ./parser "benchmark_name"
   ```

   Replace `"benchmark_name"` with the specific benchmark you want to process.

4. The Makefile will:
   - Generate the plots (`placement_before.png` and `placement_after.png`).
   - Produce the coordinates before and after placement (`coordinates_before_spreading.kiaPad` and `coordinates_after_spreading.kiaPad`).

5. Wirelength values will also be output for each benchmark.

## Output Files

- **Plots**:
  - `placement_before.png`
  - `placement_after.png`
  
- **Coordinates**:
  - `coordinates_before_spreading.kiaPad`
  - `coordinates_after_spreading.kiaPad`

- **Console Output**:
  - Total wirelength before and after placement.

## Notes

- Ensure all required files are in the same directory as the benchmark files.
- Use appropriate permissions for running the `make` and `./parser` commands.
