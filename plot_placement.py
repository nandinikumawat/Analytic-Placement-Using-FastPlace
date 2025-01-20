import matplotlib.pyplot as plt

def read_coordinates(filename):
    """Read coordinates from the file and return cells and pads as separate dictionaries."""
    cells = {}
    pads = {}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            name, x, y = parts[0], float(parts[1]), float(parts[2])
            if name.startswith('Cell'):
                cells[name] = (x, y)
            elif name.startswith('Pad'):
                pads[name] = (x, y)
    return cells, pads

def plot_coordinates(cells, pads, output_file, title):
    """Plot the coordinates for cells and pads."""
    plt.figure(figsize=(8, 8))

    # Plot cells
    for name, (x, y) in cells.items():
        plt.scatter(x, y, color='blue', label='Cells' if name == list(cells.keys())[0] else "")
        plt.text(x, y, name, fontsize=8, color='blue')

    # Plot pads
    for name, (x, y) in pads.items():
        plt.scatter(x, y, color='red', marker='x', label='Pads' if name == list(pads.keys())[0] else "")
        plt.text(x, y, name, fontsize=8, color='red')

    plt.title(title)
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.grid(True)
    plt.legend()
    plt.savefig(output_file)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python3 plot_placement.py <input_file> <output_image>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_image = sys.argv[2]

    cells, pads = read_coordinates(input_file)
    plot_coordinates(cells, pads, output_image, f"Cell and Pad Placement ({input_file})")
