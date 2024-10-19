#pragma once
//sort events by y-coord in a PriorityQueue

//Sweep status with a balanced binary search tree

/*
calculate the formula for the middle points between a point and the entire sweepline. given variables: Point - x, y Sweepline: y, variable x
formula for points p(x,y) closest to site (s) than sweepline (l):
y <= (x^2 - 2*x * xs + xs^2 + ys^2 - yl^2)/(2*ys - 2yl)

x varies.


*/

#include "utilities.h"

#include <vector>
#include <set>
#include <ostream>

static double sweepline_epsilon = 1e-9;

class voronoi_diagram {
    private:
        std::vector<point> input_points;
        std::size_t num_input_points;
        std::set<site_event> event_queue; //creates a min-heap
        site_event current_event;

        std::set<half_edge> half_edges;
        std::vector<edge> diagram_edges;
        std::vector<point> vertices;

        beachline beachline; //organize breakpoints and active sites from left to right, x=0-> x=100
        sweepline sweepline;
        int display_w = 800;
        int display_h = 600;
    public:
        voronoi_diagram(); //generates a random voronoi_diagram
        explicit voronoi_diagram(std::vector<point> input_points);
        void next_site();
        void add_circle_event(point p1,point p2,point p3, bool ordered);
        void remove_circle_event(point p1,point p2,point p3, bool ordered);
        bool remove_arc_site_at_intersection();
        void generate_half_edges_at_new_site();
        void complete_edges();
        void update_circle_event();
        void update_breakpoints();
        void update_beachline();
        void run_next_event();
        void run_voronoi();
        void display_full();
        void display_end();
};


