### Fortune's sweepline algorithm

The algorithm i have meant to create is supposed to sweep through the sites one by one and create a beachline for where the points are. it is called a beachline beacuse the way the line is calculated by 


# Voronoi Diagram with Fortune's Algorithm
A Voronoi diagram divides a plane into regions where each region belongs to the closest point (called a "site") from a set of input points. The edges of the diagram are the locations that are equally distanced between two sites.

Fortune's Algorithm is a fast and efficient method to compute Voronoi diagrams using a sweepline approach. It is supposed to operate in O(n log n) time by sweeping a horizontal line from top to bottom, handling site events and circle events as it goes. With a site event being when it has reached a input point, and a cricle event being the intersection between three lines (vertex). My algorithm has not been fully optimised yet but it is fairly fast still.

### Visuals and implementation
To make this project i have used (wikipedia)[https://en.wikipedia.org/wiki/Fortune%27s_algorithm], and a blog post by (jacquesheunis)[ https://jacquesheunis.com/post/fortunes-algorithm/]

To visualize the Voronoi diagram, I have used SDL2 for rendering, based on the (SDL2 CMake template)[https://github.com/llanillo/clion-cmake-sdl2-template].

Aside from SDL2 for the graphical interface, the algorithm is implemented purely using standard C++ libraries, ensuring efficiency without relying on external dependencies.

### Showcase
Beachline and upcoming circle-sites (blue-line)
![70beachline](https://github.com/user-attachments/assets/0973ae99-6208-499a-b2a1-490b2aec447e)

Voronoi diagram with 70 points
![voronoi 70](https://github.com/user-attachments/assets/a0ebeb76-f831-4bd0-82c3-2b2267867af7)

Voronoi diagram with 400 points
![Voronoi diagram 400 points](https://github.com/user-attachments/assets/9745d71e-722d-4e1a-bf1d-fcd2ba2e151f)

Voronoi diagram with 5000 points
![Voronoi diagram 5000 points](https://github.com/user-attachments/assets/16eac9de-0610-4f37-bffe-cc9ea8663614)


#### Further goals 
Optimise, weighted voronoi diagram, images


