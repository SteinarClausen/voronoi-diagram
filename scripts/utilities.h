//
// Created by stein on 10/5/2024.
//

#pragma once

#include <utility>
#include <vector>
#include <set>
#include <ostream>

struct vector2D;

struct point {
    double x;
    double y;

    point(const double x, const double y) : x(x), y(y) {}
    friend std::ostream& operator<<(std::ostream& os, const point& p);
    friend bool operator==(const point& p1, const point& p2);
    bool operator<(const point& other) const;
    point operator+(const vector2D& vec) const;

};

struct circle {
    point center;
    double radius;

    circle(const point c, const double r) : center(c), radius(r) {}
};

struct vector2D {
    double x;
    double y;

    vector2D(const double x, const double y) : x(x), y(y) {};

    vector2D operator-(const vector2D& other) const;
    vector2D operator+(const vector2D& other) const;
    friend std::ostream& operator<<(std::ostream& os, const vector2D& p);
    void normalize();
    vector2D operator*(double scalar) const;
    double cross(const vector2D& other) const;
    double dot(const vector2D& other) const;
};

struct edge {
    point start;
    point end;

    edge(const point start, const point end) : start(start), end(end) {};
};

struct sweepline {
    double y;
    explicit sweepline(double const y) : y(y) {}
};

struct half_edge { //when halfedges are created an opposite edge is created and is linked with the other edge
    point start;
    double sweepline_start;
    std::pair<point,point> arc_sites;
    vector2D direction;

    half_edge(const point start, const double sweepline, const vector2D& vector_2d, std::pair<point,point> points) : start(start), sweepline_start(sweepline), arc_sites(std::move(points)), direction(vector_2d) {}
    point get_breakpoint_position(double sweepline) const; //TODO: REMOVE?
    bool operator<(const half_edge& other) const;
};

class beachline {
    public:
        struct CompareByX
        {
            bool operator()(const point& lhs, const point& rhs) const {return lhs.x < rhs.x;}
        };
        std::vector<point> active_arc_sites; //arc growing from corresponding site
        std::set<point, CompareByX> breakpoints; //splits the beachline up by x-value
        std::vector<vector2D> breakpoint_vectors;
        int new_arc_site_index = 0;
        int getBreakpointPlacementIndex(const point& p) const; //returns the correct index for the active arc-site to be placed.
        friend std::ostream& operator<<(std::ostream& os, const std::set<point, CompareByX>& breakpoints);
};

std::ostream& operator<<(std::ostream& os, const std::vector<point>& active_arc_sites);

class site_event {
public:
        point site;
        bool isCircleEvent; //when three site-lines intersect
        double y;
        std::vector<point> circlePoints;
        double radius=0;
    public:
        site_event() : site(0,0), isCircleEvent(false) , y(0.0), circlePoints({}){}
        site_event(const point p, const bool isSiteEvent, const double y, const std::vector<point>& points={}): site(p), isCircleEvent(isSiteEvent), y(y), circlePoints(points) {};
        bool operator<(const site_event& other) const;
        site_event& operator=(const site_event& other);
        double getY() const {return y;}
        point getSite() const {return site;}
        bool getIsCircleEvent() const {return isCircleEvent;}
        point getCirclePoints(const int it) const {return circlePoints.at(it);}
};

double calculate_y_parabola(double x_parabola, double x_site, double y_site, double y_sweepline);

double calculate_y_parabola_derivative(double x_parabola,double x_site,double y_site,double y_sweepline);


double calculate_parabola_intersection(point a, point b, double y_sweepline);

//mirror a point on the line AB - useful for getting "sister" circle-sites incase it is
point mirror_point(point mirror_point, point A, point B);

circle circumcircle(point A, point B, point C);
