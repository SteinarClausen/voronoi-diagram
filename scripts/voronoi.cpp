//
// Created by stein on 9/2/2024.
//

#include "voronoi.h"
#include "utilities.h"
#include <SDL.h>
#include <cmath>
#include <vector>
#include <random>

#include <iostream>

voronoi_diagram::voronoi_diagram(std::vector<point> input_points) : input_points(std::move(input_points)), num_input_points(input_points.size()), sweepline(0.0){
    for (const point& p : this->input_points)
    {
        event_queue.emplace(p,false,p.y);
    }
}

voronoi_diagram::voronoi_diagram() : voronoi_diagram(std::vector<point>{}) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dist_x(0.0, display_w);
    std::uniform_real_distribution<> dist_y(0.0, display_h);

    constexpr int points_int = 500;
    std::vector<point> points;

    points.reserve(points_int);
    for (int i = 0; i < points_int; i++) {
        points.emplace_back(dist_x(gen), dist_y(gen));
    }
    for (int i = 0; i < points_int; i++)
    {
        std::cout << points[i] << " ";
    }
    std::cout << std::endl;
    *this = voronoi_diagram(points);

}

void voronoi_diagram::next_site() {
    if (event_queue.empty())
    {
        std::cerr << "Error in voronoi_diagram::next_site(), event_queue empty" << std::endl;
    }
    if (event_queue.begin()->getSite().y>display_h+1)
    {
        event_queue.erase(event_queue.begin());
        if (!event_queue.empty())
        {
            next_site();
        }
        return;
    }
    if (event_queue.begin()->getSite().x>display_w+1 || event_queue.begin()->getSite().x<-1)
    {
        event_queue.erase(event_queue.begin());
        if (!event_queue.empty())
        {
            next_site();
        }
        return;
    }
    current_event = *event_queue.begin();
    sweepline.y = current_event.y + sweepline_epsilon; //TODO: Add sweepline_epsilon when using it
    event_queue.erase(event_queue.begin());
}

//TODO: Make a smaller function that contains the point sorting.
void voronoi_diagram::add_circle_event(point p1, point p2, point p3, const bool ordered = false) { //order matters - i,i+1,i+2
    if (p1==p2 || p2==p3 || p3==p1)
    {
        return;
    }
    const circle c = circumcircle(p1, p2, p3);
    const point p = {c.center.x, c.center.y};
    if (!ordered)
    {
        //need to get the points in the right order, to do this one can use the derivative of the function at the point.x and with a sweepline at point.y.
        //the arcsite to remove is the median derivative.
        const double p1_derivative = calculate_y_parabola_derivative(p.x,p1.x,p1.y,p.y + c.radius);
        const double p2_derivative = calculate_y_parabola_derivative(p.x,p2.x,p2.y,p.y + c.radius);
        const double p3_derivative = calculate_y_parabola_derivative(p.x,p3.x,p3.y,p.y + c.radius);

        //sort the point so it goes in the correct order, lowest positive number to the highest number, middle value is the one removed
        if(p1_derivative>p2_derivative){
            std::swap(p1,p2);
        }
        if(p1_derivative>p3_derivative){
            std::swap(p1,p3);
        }
        if(p2_derivative>p3_derivative){
            std::swap(p2,p3);
        }
    }

    if(p.y+c.radius > sweepline.y) { //TODO: optimise
        for (auto it = event_queue.begin(); it != event_queue.end();)
        {
            if (it->circlePoints.empty()) {
                ;
            } else if(it->circlePoints.at(0)==p1 && it->circlePoints.at(1)==p2 && it->circlePoints.at(2)==p3) {
                return;
            }
            ++it;
        }
        site_event site_event{p, true,p.y + c.radius, {p1,p2,p3}};
        site_event.radius = c.radius;
        event_queue.insert(site_event); //order matters here
    }
}

void voronoi_diagram::remove_circle_event(point p1, point p2, point p3, bool ordered = false) {
    if (p1==p2 || p2==p3 || p3==p1)
    {
        return;
    }
    const circle c = circumcircle(p1, p2, p3);
    const point p = {c.center.x, c.center.y + c.radius};

    if (!ordered)
    {
        //need to get the points in the right order, to do this one can use the derivative of the function at the point.x and with a sweepline at point.y.
        //the arcsite to remove is the median derivative.
        const double p1_derivative = calculate_y_parabola_derivative(p.x,p1.x,p1.y,p.y + c.radius);
        const double p2_derivative = calculate_y_parabola_derivative(p.x,p2.x,p2.y,p.y + c.radius);
        const double p3_derivative = calculate_y_parabola_derivative(p.x,p3.x,p3.y,p.y + c.radius);

        //sort the point so it goes in the correct order, lowest positive number to the highest number, middle value is the one removed
        if(p1_derivative>p2_derivative){
            std::swap(p1,p2);
        }
        if(p1_derivative>p3_derivative){
            std::swap(p1,p3);
        }
        if(p2_derivative>p3_derivative){
            std::swap(p2,p3);
        }
    }
    if(p.y > sweepline.y) { //TODO: optimise
        for (auto it = event_queue.begin(); it != event_queue.end();)
        {
            if (it->circlePoints.empty()) {
                ;
            } else if(it->circlePoints.at(0)==p1 && it->circlePoints.at(1)==p2 && it->circlePoints.at(2)==p3) {
                event_queue.erase(it);
                return;
            }
            ++it;
        }
    }
}

void voronoi_diagram::update_circle_event() //TODO: optimise
{
    if (beachline.active_arc_sites.size() < 5)
    {
        return;
    }
    size_t const active_size = beachline.active_arc_sites.size();
    int const index = beachline.new_arc_site_index;
    if (index == 1) {
        add_circle_event(beachline.active_arc_sites.at(index),beachline.active_arc_sites.at(index+1),beachline.active_arc_sites.at(index+2));
    }
    else if (beachline.new_arc_site_index == active_size-2)
    {
        add_circle_event(beachline.active_arc_sites.at(index-2),beachline.active_arc_sites.at(index-1),beachline.active_arc_sites.at(index));
    }
    else
    {
        remove_circle_event(beachline.active_arc_sites.at(index-2), beachline.active_arc_sites.at(index+1), beachline.active_arc_sites.at(index+2), true);
        add_circle_event(beachline.active_arc_sites.at(index-2), beachline.active_arc_sites.at(index-1), beachline.active_arc_sites.at(index));
        add_circle_event(beachline.active_arc_sites.at(index), beachline.active_arc_sites.at(index+1), beachline.active_arc_sites.at(index+2));
    }
}

//if you take the derivative of the function in all the active_arcsite points in the circle_event point
//you get three different derivatives, the middle value of these is the arc_site to remove
//need to be done before breakpoints clear
bool voronoi_diagram::remove_arc_site_at_intersection() {
    const point p1 = current_event.getCirclePoints(0);
    const point p2 = current_event.getCirclePoints(1);
    const point p3 = current_event.getCirclePoints(2);

    // Reverse iterate through active_arc_sites
    for (auto it = beachline.active_arc_sites.rbegin(); it != beachline.active_arc_sites.rend() - 2;) {
        if (*it == p3 && *(it + 1) == p2 && *(it + 2) == p1) {
            // Erase the middle site (p2)
            auto const forward_it = beachline.active_arc_sites.erase((it+2).base());

            if (forward_it-1 != beachline.active_arc_sites.begin()) {
                remove_circle_event(*(forward_it - 2), *(forward_it-1), p2, true);
                add_circle_event(*(forward_it - 2), *(forward_it-1), *(forward_it));
            } //pain to debug
            if (forward_it != beachline.active_arc_sites.end()) {
                remove_circle_event(p2, *forward_it, *(forward_it + 1),true);
                add_circle_event(*(forward_it-1), *(forward_it), *(forward_it + 1));
            }
            return true;
        } else {
            ++it;
        }
    }
    return false;
}

void voronoi_diagram::generate_half_edges_at_new_site()
{
    const point p1 = current_event.getCirclePoints(0);
    const point p2 = current_event.getCirclePoints(1); // the one removed
    const point p3 = current_event.getCirclePoints(2);

    const std::pair<point, point> arc1(std::min(p1, p2), std::max(p1, p2)); //old breakpoint
    const std::pair<point, point> arc2(std::min(p2, p3), std::max(p2, p3)); //old breakpoint
    const std::pair<point, point> arc3(std::min(p1, p3), std::max(p1, p3)); //new breakpoint

    bool arc1_removed = false;
    bool arc2_removed = false;

    const double temp = current_event.getSite().y;
    current_event.site.y = current_event.site.y - current_event.radius;
    for (auto it = half_edges.begin(); it != half_edges.end(); )
    {
        if (it->arc_sites == arc1) {
            diagram_edges.emplace_back(it->start, current_event.site);
            it = half_edges.erase(it); // Safely erase and get the next iterator
            arc1_removed = true;
        } else if (it->arc_sites == arc2) {
            diagram_edges.emplace_back(it->start, current_event.site);
            it = half_edges.erase(it); // Safely erase and get the next iterator
            arc2_removed = true;
        } else {
            ++it;
        }
    }
    if (!arc1_removed) {
        const double x = calculate_parabola_intersection(p2, p1, sweepline.y);
        const double y = calculate_y_parabola(x, p1.x, p1.y, sweepline.y);
        vector2D direction(x-current_event.getSite().x, y-current_event.getSite().y);
        direction.normalize();
        half_edges.emplace(current_event.getSite(), sweepline.y, direction, arc1);
    }

    if (!arc2_removed) {
        const double x = calculate_parabola_intersection(p3, p2, sweepline.y);
        const double y = calculate_y_parabola(x, p3.x, p2.y, sweepline.y);
        vector2D direction(x-current_event.getSite().x, y-current_event.getSite().y);
        direction.normalize();
        half_edges.emplace(current_event.getSite(), sweepline.y, direction, arc2);
    }

    // Create a new half-edge for the third intersection
    const double x = calculate_parabola_intersection(p3, p1, sweepline.y);
    const double y = calculate_y_parabola(x, p1.x, p1.y, sweepline.y);
    vector2D direction(current_event.getSite().x-x, current_event.getSite().y-y);
    direction.normalize();
    half_edges.emplace(current_event.getSite(), sweepline.y, direction, arc3);
    current_event.site.y = temp;
}

void voronoi_diagram::complete_edges()
{
    for (auto it = half_edges.begin(); it != half_edges.end(); )
    {
        const point start = it->start;
        const vector2D direction = it->direction; //multiply the direction with t1 or t2 where it is out of frame and finish the lines
        double t1;
        double t2;
        if (direction.x < 0)
        {
            t1 = (-1-start.x)/direction.x;
        } else
        {
            t1 = (display_w+1 - start.x)/direction.x;
        }

        if (direction.y < 0)
        {
            t2 = (-1-start.y)/direction.y;
        } else
        {
            t2 = (display_h+1 - start.y)/direction.y;
        }
        point end = start + direction * std::min(t1,t2);
        if (start.x < display_w && start.y < display_h && start.x > 0 && start.y > 0)
        {
            diagram_edges.emplace_back(start, end);
        }
        it = half_edges.erase(it);
    }
}

void voronoi_diagram::update_beachline() {
    if (current_event.getIsCircleEvent())
    {
        if (remove_arc_site_at_intersection())
        {
            generate_half_edges_at_new_site();
        }
    }
    else if (!beachline.active_arc_sites.empty()) //if it is a regular site event
    {
        //need to place the new site in active_arc_sites at the correct index given the breakline x-values
        const int Index = beachline.getBreakpointPlacementIndex(current_event.getSite());
        //split old active_site_beachline in two
        if (Index >= 0 && Index <= beachline.active_arc_sites.size()) {
            // Get the iterator to the position where 'Index' points
            auto iter = beachline.active_arc_sites.begin() + Index;

            beachline.new_arc_site_index=Index+1; //need to give an index to where the new_arc_site will land
            // Insert the new site at the position right after the site at 'Index'
            iter = beachline.active_arc_sites.insert(iter + 1, current_event.getSite());

            // Insert a duplicate of the site at 'Index' right after the new site
            beachline.active_arc_sites.insert(iter + 1, beachline.active_arc_sites.at(Index));
            update_circle_event();
        } else {
            std::cerr << "Index out of bounds: " << Index << " Size: " << beachline.active_arc_sites.size() << std::endl;
        }
    }
    else if (beachline.active_arc_sites.empty()) //if the arc_sites vector is empty
    {
        beachline.active_arc_sites.push_back(current_event.getSite());
        beachline.breakpoints.emplace(display_w,0);
    } else
    {
        std::cerr << "Beachline encountered undefined behaviour" << std::endl;
    }
    update_breakpoints();
    //update_vertices

}

void voronoi_diagram::update_breakpoints() {
    beachline.breakpoints.clear();

    for (int i = 0; i < beachline.active_arc_sites.size()-1; i++) {
        double breakpoint_x = calculate_parabola_intersection(beachline.active_arc_sites.at(i), beachline.active_arc_sites.at(i+1), sweepline.y);
        double breakpoint_y = calculate_y_parabola(breakpoint_x, beachline.active_arc_sites.at(i).x, beachline.active_arc_sites.at(i).y, sweepline.y);
        beachline.breakpoints.emplace(breakpoint_x,breakpoint_y); //TODO: finne ut om x-verdier må legges inn i på en annen måte for å gjøre opp for at det skal være et balanced binary tree..

    }
}

void voronoi_diagram::run_next_event()
{
    if(!event_queue.empty())
    {
        next_site();
        if (!beachline.active_arc_sites.empty())
        {
            update_breakpoints(); //TODO: Do this in a smarter way, possibly with vectors with len 1 and multiply by height difference or something
        }
        update_beachline();
    }
    else
    {
        std::cerr << "Event queue is empty" << std::endl;
    }
}

void voronoi_diagram::run_voronoi() {
    while (!event_queue.empty())
    {
        run_next_event();
    }
    complete_edges();
}

void voronoi_diagram::display_full() {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow("Blank Window", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, display_w, display_h, SDL_WINDOW_SHOWN);
    if (window == nullptr) {
        SDL_Log("Failed to create window: %s", SDL_GetError());
        return;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == nullptr) {
        SDL_Log("Failed to create renderer: %s", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return;
    }

    bool running = true;
    bool next_event = false;
    bool key_c_pressed = false;
    SDL_Event event;

    while (running) {
        // Close window with any input
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT || event.type == SDL_MOUSEBUTTONDOWN) {
                running = false;
                break;
            } if (event.key.keysym.sym == SDLK_c)
            {
                next_event = true;
                key_c_pressed = true;
            }
            if (event.type == SDL_KEYUP)
            {
                if (event.key.keysym.sym == SDLK_c)
                {
                    key_c_pressed = false;
                }
            }
        }
        if (next_event && !key_c_pressed)
        {
            if(!event_queue.empty())
            {
                run_next_event();
            }
            next_event = false;
        }
        if (!event_queue.empty())
        {
            run_next_event();
        }
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // Set the background color to purple
        SDL_RenderClear(renderer);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // Set the color to white
        const int height = static_cast<int>(sweepline.y);
        SDL_RenderDrawLine(renderer, 1, static_cast<int>(height), 800, static_cast<int>(height)); // Draws a line from (x1, y1) to (x2, y2)
        //int index = 0;
        //int iteratorIndex = 0;
        // Loop over breakpoints in the beachline set using iterators
        //Skal ha en løkke som går gjennom alle breakpointsene, basert på hvilken breakpoint det er skal det bestemme hvilken active arc site som gjelder
        //Skal starte på 0 og slutte på 800
        auto it = beachline.breakpoints.begin();
        int index = 0;
        int iteratorIndex = 0;
        double previous_y;
        double current_y;
        while (index<this->display_w)
        {
            if(beachline.breakpoints.empty()) {
                break;
            }
            if (index>=it->x)
            {

                previous_y = calculate_y_parabola(static_cast<double>(index-1),beachline.active_arc_sites[iteratorIndex].x,beachline.active_arc_sites[iteratorIndex].y, height);
                while(index>=it->x != (it == beachline.breakpoints.end()))
                {
                    ++it;
                    iteratorIndex++;
                    if (it->x < index && it->x > 0)
                    {
                        SDL_RenderDrawLine(renderer, index, static_cast<int>(previous_y), index, height);
                    }
                }
                current_y = calculate_y_parabola(static_cast<double>(index),beachline.active_arc_sites[iteratorIndex].x,beachline.active_arc_sites[iteratorIndex].y, height);
                SDL_RenderDrawLine(renderer, index-1, static_cast<int>(previous_y), index, static_cast<int>(current_y));
                index++;
            } else
            {
                previous_y = calculate_y_parabola(static_cast<double>(index-1),beachline.active_arc_sites[iteratorIndex].x,beachline.active_arc_sites[iteratorIndex].y, height);
                current_y = calculate_y_parabola(static_cast<double>(index),beachline.active_arc_sites[iteratorIndex].x,beachline.active_arc_sites[iteratorIndex].y, height);
                if (!(current_y < 0 && previous_y <0) && (current_y < display_h && previous_y < display_h))
                {
                    SDL_RenderDrawLine(renderer, index-1, static_cast<int>(previous_y), index, static_cast<int>(current_y));
                }
                index++;
            }
        }
        for (const point p : input_points)
        {
            SDL_RenderDrawPoint(renderer, static_cast<int>(p.x),static_cast<int>(p.y));
            SDL_RenderDrawPoint(renderer, static_cast<int>(p.x),static_cast<int>(p.y)+1);
            SDL_RenderDrawPoint(renderer, static_cast<int>(p.x),static_cast<int>(p.y)-1);
            SDL_RenderDrawPoint(renderer, static_cast<int>(p.x)+1,static_cast<int>(p.y));
            SDL_RenderDrawPoint(renderer, static_cast<int>(p.x)-1,static_cast<int>(p.y));
        }
        if (!event_queue.empty())
        {
            SDL_SetRenderDrawColor(renderer, 30, 40, 255, 255);
            auto iter = event_queue.begin();
            while (iter !=event_queue.end())
            {
                if (iter->getIsCircleEvent())
                {
                    SDL_RenderDrawLine(renderer, 0,static_cast<int>(iter->y),display_w,static_cast<int>(iter->y));
                }
                ++iter;
            }
        }
        SDL_SetRenderDrawColor(renderer, 255, 40, 255, 255);
        if (!half_edges.empty())
        {
            for (auto iterer = half_edges.begin(); iterer != half_edges.end(); ++it)
            {
                SDL_RenderDrawLine(renderer, static_cast<int>(iterer->start.x),static_cast<int>(iterer->start.y),static_cast<int>(iterer->start.x + iterer->direction.x*10),static_cast<int>(iterer->start.y + iterer->direction.y*10));
                ++iterer;
            }
        }

        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); // Set the color to red
        if(!diagram_edges.empty())
        {
            if (event_queue.empty() && !half_edges.empty())
            {
                complete_edges();
            }
            for (const auto & diagram_edge : diagram_edges)
            {
                SDL_RenderDrawLine(renderer, static_cast<int>(diagram_edge.start.x),static_cast<int>(diagram_edge.start.y),static_cast<int>(diagram_edge.end.x),static_cast<int>(diagram_edge.end.y));
            }
        }
        SDL_RenderPresent(renderer);
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void voronoi_diagram::display_end()
{
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow("Voronoi Diagram", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, display_w, display_h, SDL_WINDOW_SHOWN);
    if (window == nullptr) {
        SDL_Log("Failed to create window: %s", SDL_GetError());
        return;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == nullptr) {
        SDL_Log("Failed to create renderer: %s", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return;
    }

    bool running = true;
    SDL_Event event;

    bool voronoi_drawn = false;

    SDL_Texture* voronoi_texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, display_w, display_h);

    run_voronoi();

    SDL_SetRenderTarget(renderer, voronoi_texture);

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    for (const point& p : input_points)
    {
        SDL_RenderDrawPoint(renderer, static_cast<int>(p.x), static_cast<int>(p.y));
        SDL_RenderDrawPoint(renderer, static_cast<int>(p.x), static_cast<int>(p.y) + 1);
        SDL_RenderDrawPoint(renderer, static_cast<int>(p.x), static_cast<int>(p.y) - 1);
        SDL_RenderDrawPoint(renderer, static_cast<int>(p.x) + 1, static_cast<int>(p.y));
        SDL_RenderDrawPoint(renderer, static_cast<int>(p.x) - 1, static_cast<int>(p.y));
    }

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    for (const auto& diagram_edge : diagram_edges)
    {
        SDL_RenderDrawLine(renderer, static_cast<int>(diagram_edge.start.x), static_cast<int>(diagram_edge.start.y),
                           static_cast<int>(diagram_edge.end.x), static_cast<int>(diagram_edge.end.y));
    }

    SDL_SetRenderTarget(renderer, nullptr);

    while (running)
    {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT || event.type == SDL_MOUSEBUTTONDOWN) {
                running = false;
                break;
            }
        }

        SDL_RenderCopy(renderer, voronoi_texture, nullptr, nullptr);
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyTexture(voronoi_texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}


