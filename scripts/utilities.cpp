//
// Created by stein on 10/5/2024.
//

#include "utilities.h"

#include <vector>
#include <iostream>
#include <cmath>

std::ostream& operator<<(std::ostream& os, const point& p) {
    os << "(" << p.x << "," << p.y << ")";
    return os;
}

bool operator==(const point& p1, const point& p2)
{
    if (p1.x == p2.x && p1.y == p2.y)
    {
        return true;
    }
    return false;
}

bool point::operator<(const point& other) const
{
    return x < other.x;
}


point point::operator+(const vector2D& vec) const
{
    return {x + vec.x, y + vec.y};
}

vector2D vector2D::operator*(const double scalar) const
{
    return {x * scalar, y * scalar};
}

vector2D vector2D::operator+(const vector2D& other) const
{
    return {x+other.x,y+other.x};
}

vector2D vector2D::operator-(const vector2D& other) const
{
    return {x-other.x,y-other.x};
}

std::ostream& operator<<(std::ostream& os, const vector2D& p)
{
    os << "(" << p.x << "," << p.y << ")";
    return os;
}

double vector2D::cross(const vector2D& other) const
{
    return x * other.y - y * other.x;
}

double vector2D::dot(const vector2D& other) const
{
    return x*other.x + y*other.y;
}


void vector2D::normalize()
{
    double const length = sqrt(x*x+y*y);
    if (length == 0)
    {
        return;
    }
    this->x = x / length;
    this->y = y / length;
}

point half_edge::get_breakpoint_position(const double sweepline) const
{
    const double scalar = (sweepline-sweepline_start)*0.5;
    const point position = start+(direction*scalar);
    return position;
}

bool half_edge::operator<(const half_edge& other) const
{
    if (start.x == other.start.x)
    {
        return direction.x < other.direction.x;
    }
    return start.x < other.start.x;
}



int beachline::getBreakpointPlacementIndex(const point& p) const {
    int index = 0;

    for (auto it = breakpoints.begin(); it != breakpoints.end(); ++it, ++index) { // Iterate through breakpoints and find the correct index
        if (it->x > p.x) {
            return index; // Return the index of the first point with x > pt.x
        }
    }
    return index;     // If no point is found with x > pt.x, return the index of the last element + 1
}

std::ostream& operator<<(std::ostream& os, const std::vector<point>& active_arc_sites) {
    os << "[";
    for (const point& p : active_arc_sites)
    {
        os << p;
    }
    os << "]" << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::set<point, beachline::CompareByX>& breakpoints) {
    os << "{";
    for (auto it = breakpoints.begin(); it != breakpoints.end(); ++it) {
        if (it != breakpoints.begin()) os << ", ";
        os << *it;
    }
    os << "}";
    return os;
}

site_event& site_event::operator=(const site_event& other) {
    if (this != &other) { // Prevent self-assignment
        site = other.site;
        isCircleEvent = other.isCircleEvent;
        y = other.y;
        circlePoints = other.circlePoints;
    }
    return *this;
}

bool site_event::operator<(const site_event& other) const
{
    if (y != other.y)
        return y < other.y;
    return site.x<other.site.x;
}


double calculate_y_parabola(const double x_parabola, const double x_site, const double y_site, const double y_sweepline) {
    return (x_parabola*x_parabola - 2*x_parabola * x_site + x_site*x_site + y_site*y_site - y_sweepline*y_sweepline)/(2*y_site - 2*y_sweepline);
}

double calculate_y_parabola_derivative(const double x_parabola, const double x_site, const double y_site, const double y_sweepline) {
    return (2*x_parabola-2*x_site)/(2*y_site - 2*y_sweepline);
}

double calculate_parabola_intersection(const point a, const point b, const double y_sweepline)
{
    const double a_x_site = a.x;
    const double a_y_site = a.y;
    const double b_x_site = b.x;
    const double b_y_site = b.y;

    const double A = 2*(b_y_site-a_y_site);
    if (A == 0)
    {
        return std::min(a_x_site,b_x_site) + (std::max(a_x_site, b_x_site) - std::min(a_x_site, b_x_site))/2; //avoid division by zero, take middle point between the two
    }
    const double B = 4*(a_y_site*b_x_site-b_y_site*a_x_site+y_sweepline*a_x_site-y_sweepline*b_x_site);
    const double C = (a_x_site*a_x_site + a_y_site*a_y_site - y_sweepline*y_sweepline)*(2*b_y_site - 2*y_sweepline) -
               (b_x_site*b_x_site + b_y_site*b_y_site - y_sweepline*y_sweepline)*(2*a_y_site - 2*y_sweepline);

    const double discriminant = std::abs(B*B-4*A*C); //to make sure there is no error with the sqrt of a negative number, only happens with very small numbers so it doesn't matter

    const double sqrt_discriminant = std::sqrt(discriminant);
    const double intersection_x1 = (-B + sqrt_discriminant) / (2 * A)+1e-5;
    const double intersection_x2 = (-B - sqrt_discriminant) / (2 * A)-1e-5;

    if((a_x_site < b_x_site && a_y_site > b_y_site) || (a_x_site > b_x_site && a_y_site > b_y_site)) //TODO: bare a_y > b_y
    {
        return std::max(intersection_x1,intersection_x2);
    } else
    {
        return std::min(intersection_x1,intersection_x2);
    }
}

point mirror_point(const point mirror_point, const point A, const point B)
{
    point projection(0,0);

    const double ABx = B.x - A.x;
    const double ABy = B.y - A.y;

    const double APx = mirror_point.x - A.x;
    const double APy = mirror_point.y - A.y;

    const double dot_product_AP_AB = APx * ABx + APy * ABy;
    const double dot_product_AB_AB = ABx * ABx + ABy * ABy;

    // Projection coordinates
    projection.x = A.x + (dot_product_AP_AB / dot_product_AB_AB) * ABx;
    projection.y = A.y + (dot_product_AP_AB / dot_product_AB_AB) * ABy;

    const point mirrored = {(2*projection.x-mirror_point.x),(2*projection.y-mirror_point.y)};
    return mirrored;
}

circle circumcircle(const point A, const point B, const point C) {
    const double determinant = (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));

    if (determinant == 0) {
        throw std::invalid_argument("The points are colinear, no circumcircle exists");
    }

    const double circumcenter_x = ((A.x*A.x + A.y*A.y)*(B.y-C.y) + (B.x*B.x + B.y*B.y)*(C.y-A.y) + (C.x*C.x + C.y*C.y)*(A.y-B.y)) / (2*determinant);
    const double circumcenter_y = ((A.x*A.x + A.y*A.y)*(C.x-B.x) + (B.x*B.x + B.y*B.y)*(A.x-C.x) + (C.x*C.x + C.y*C.y)*(B.x-A.x)) / (2*determinant);

    const point circumcenter(circumcenter_x,circumcenter_y);

    const double radius = std::sqrt((A.x-circumcenter_x) * (A.x-circumcenter_x) + (A.y-circumcenter_y)*(A.y-circumcenter_y));
    const circle circumcircle(circumcenter,radius);

    return circumcircle;
}