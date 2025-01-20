CFLAGS=-g

all: placer generate_coordinates plot_before plot_after

placer: suraj_parser.o main.o
	g++ $(CFLAGS) -o placer suraj_parser.o main.o

suraj_parser.o: suraj_parser.cpp
	g++ $(CFLAGS) -c suraj_parser.cpp

main.o: main.cpp suraj_parser.h
	g++ $(CFLAGS) -c main.cpp

# Generate the .kiaPad files using the placer program
generate_coordinates: placer
	./placer toy01  # Replace 'toy01' with your input circuit file prefix

# Generate plot for coordinates before spreading
plot_before: coordinates_before_spreading.kiaPad plot_placement.py
	python3 plot_placement.py coordinates_before_spreading.kiaPad placement_before.png

# Generate plot for coordinates after spreading
plot_after: coordinates_after_spreading.kiaPad plot_placement.py
	python3 plot_placement.py coordinates_after_spreading.kiaPad placement_after.png

# Clean up generated files
clean:
	rm -f *.o placer coordinates_before_spreading.kiaPad coordinates_after_spreading.kiaPad placement_before.png placement_after.png

.PHONY: all clean generate_coordinates plot_before plot_after


