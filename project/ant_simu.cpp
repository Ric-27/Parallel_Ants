#include <vector>
#include <iostream>
#include <random>
#include "labyrinthe.hpp"
#include "ant.hpp"
#include "pheromone.hpp"
# include "gui/context.hpp"
# include "gui/colors.hpp"
# include "gui/point.hpp"
# include "gui/segment.hpp"
# include "gui/triangle.hpp"
# include "gui/quad.hpp"
# include "gui/event_manager.hpp"
# include "display.hpp"

#include <mpi.h>
#include <chrono>

std::chrono::time_point<std::chrono::system_clock> start[5], end[5];
std::chrono::duration<double> elapsed_seconds[5];

size_t food_quantity = 0;

void print_food_quantity(){
    std::cout << "there is " << food_quantity << " food in the nest" << std::endl;
}

void debug(std::string c = "®"){
    std::cout << c << std::endl;
}

void print_function_time(){
    std::cout << "function advance has been calculated in: "<< elapsed_seconds[1].count() << " seg" << std::endl;
    std::cout << "function evaporation has been calculated in: "<< elapsed_seconds[2].count() << " seg" << std::endl;
    std::cout << "function update has been calculated in: "<< elapsed_seconds[3].count() << " seg" << std::endl;
    std::cout << "function display has been calculated in: "<< elapsed_seconds[4].count() << " seg" << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
}

void advance_time( const labyrinthe& land, pheromone& phen, 
                   const position_t& pos_nest, const position_t& pos_food,
                   std::vector<ant>& ants, std::size_t& cpteur )
{
    start[1] = std::chrono::system_clock::now();

    for ( size_t i = 0; i < ants.size(); ++i ){
        ants[i].advance(phen, land, pos_food, pos_nest, cpteur);    
    }

    end[1] = std::chrono::system_clock::now();
    elapsed_seconds[1] = end[1] - start[1];
    //----------------------------------------------------------    
    start[2] = std::chrono::system_clock::now();

    phen.do_evaporation();

    end[2] = std::chrono::system_clock::now();
    elapsed_seconds[2] = end[2] - start[2];    
    //----------------------------------------------------------
    start[3] = std::chrono::system_clock::now();

    phen.update();

    end[3] = std::chrono::system_clock::now();
    elapsed_seconds[3] = end[3] - start[3];
}

int main(int nargs, char* argv[])
{
    bool validation = false;
    bool debug_validation = false;

    int core_number, total_of_cores;

    const dimension_t dims{32,64};// Dimension du labyrinthe
    const std::size_t life = int(dims.first*dims.second);
    const int nb_ants = 2*dims.first*dims.second; // Nombre de fourmis
    
    const double alpha=0.97; // Coefficient de chaos
    //const double beta=0.9999; // Coefficient d'évaporation
    const double beta=0.999; // Coefficient d'évaporation
                             // 
    labyrinthe laby(dims);
    // Location du nid
    position_t pos_nest{dims.first/2,dims.second/2};
    // Location de la nourriture
    position_t pos_food{dims.first-1,dims.second-1};

    MPI_Init(&nargs, &argv);   
    
    MPI_Comm_rank(MPI_COMM_WORLD, &core_number);
    MPI_Comm_size(MPI_COMM_WORLD, &total_of_cores);
    //
    MPI_Group world_group, display_group, calcul_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Comm display_comm, calcul_comm;

    const int core_0_1[2] = {0, 1};
    int core_1_n[total_of_cores-1];
    for (int i = 1; i < total_of_cores; i++)
    {
        core_1_n[i - 1] = i;
    }
    if(core_number <= 1) {
        MPI_Group_incl(world_group, 2, core_0_1, &display_group);
        MPI_Comm_create_group(MPI_COMM_WORLD, display_group, 0, &display_comm);
    }
    if(core_number > 0) {
        MPI_Group_incl(world_group, total_of_cores - 1, core_1_n, &calcul_group);
        MPI_Comm_create_group(MPI_COMM_WORLD, calcul_group, 1, &calcul_comm);
    }
    //
    int nb_ants_loc = nb_ants / (total_of_cores - 1);
    int buffer_size = 1 + 2 * (nb_ants_loc * (total_of_cores - 1)) + laby.dimensions().first*laby.dimensions().second;
    
    //--------------------------------------------------------------------
    if (core_number >= 1)
    { 
        const double eps = 0.75;  // Coefficient d'exploration
        
        // Définition du coefficient d'exploration de toutes les fourmis.
        ant::set_exploration_coef(eps);
        std::vector<ant> ants;

        ants.reserve(nb_ants_loc);  

        for ( size_t i = 0; i < nb_ants_loc; ++i ) {
            ants.emplace_back(pos_nest, life);
        }
        // On crée toutes les fourmis dans la fourmilière.
        pheromone phen(laby.dimensions(), pos_food, pos_nest, alpha, beta);
        size_t food_loc = 0;    
        while(true) {
            int food_total = 0;
    
            std::vector<double> buffer_ants(2*nb_ants_loc*(total_of_cores - 1));
            std::vector<double> buffer_ants_loc;

            std::vector<double> buffer_phen(laby.dimensions().first*laby.dimensions().second);
            std::vector<double> buffer_phen_loc;

            for (size_t i = 0; i < nb_ants_loc; i++)
            {
                position_t pos_ant = ants[i].get_position();
                buffer_ants_loc.emplace_back((double)pos_ant.first);
                buffer_ants_loc.emplace_back((double)pos_ant.second);
            }

            for (std::size_t i = 0; i < laby.dimensions().first; i++){
                for (std::size_t j = 0; j < laby.dimensions().second; j++)
                {
                    buffer_phen_loc.emplace_back((double)phen(i,j));
                }                 
            }
            MPI_Reduce(&food_loc,&food_total,1,MPI_INT,MPI_SUM,0,calcul_comm);
            //debug("reduce");
            MPI_Gather(buffer_ants_loc.data(),buffer_ants_loc.size(),MPI_DOUBLE,buffer_ants.data(),buffer_ants_loc.size(),MPI_DOUBLE,0,calcul_comm);
            //debug("gather");
            MPI_Allreduce(buffer_phen_loc.data(),buffer_phen.data(),buffer_phen.size(),MPI_DOUBLE,MPI_MAX,calcul_comm);
            //debug("all reduce");
            phen.swap_map(buffer_phen);
            if(core_number==1){
                std::vector<double> buffer;
                buffer.emplace_back((double)food_total);
                buffer.insert(buffer.end(),buffer_ants.begin(),buffer_ants.end());
                buffer.insert(buffer.end(),buffer_phen.begin(),buffer_phen.end());
                
                /*if(!debug_validation){
                    std::cout << food_total << std::endl;
                    debug_validation = false;
                }*/

                MPI_Request request;
                MPI_Status status;
                MPI_Isend(buffer.data(),buffer.size(),MPI_DOUBLE,0,101,display_comm,&request);
                advance_time(laby, phen, pos_nest, pos_food, ants, food_loc); 
                MPI_Wait(&request,&status);
            } else {
                advance_time(laby, phen, pos_nest, pos_food, ants, food_loc);
            }
            
        }
    }
    //------------------------------------------------------------

    if(core_number == 0){
        start[0] = std::chrono::system_clock::now();

        int pher_start = 2 * (nb_ants_loc * (total_of_cores - 1)) + 1;

        MPI_Status(status);

        std::vector<ant> ants;
        ants.reserve(nb_ants);

        pheromone phen(laby.dimensions(), pos_food, pos_nest, alpha, beta);

        std::vector<double> buffer(buffer_size);
        
        MPI_Recv(buffer.data(),buffer.size(),MPI_DOUBLE,1,101,display_comm,&status);
        //debug("recibi la primera");        
        food_quantity = buffer[0];
        for (size_t i = 1; i < pher_start; i += 2)
        {
            ants.emplace_back(position_t(buffer[i],buffer[i+1]),life);
            //std::cout << buffer[i] << " " << buffer[i+1] << std::endl;
        }

        phen.swap_map(std::vector<double>(buffer.begin() + pher_start,buffer.end()));        

        gui::context graphic_context(nargs, argv);
        gui::window& win =  graphic_context.new_window(h_scal*laby.dimensions().second,h_scal*laby.dimensions().first+266);

        display_t displayer( laby, phen, pos_nest, pos_food, ants, win );
        
        gui::event_manager manager;
        manager.on_key_event(int('q'), [] (int code) { MPI_Abort(MPI_COMM_WORLD,MPI_SUCCESS); });

        manager.on_key_event(int('f'), [] (int code) { print_food_quantity(); });

        manager.on_key_event(int('t'), [] (int code) { print_function_time(); });

        manager.on_display([&] { displayer.display(food_quantity); win.blit(); });

        manager.on_idle([&] () {
                   
            //----------------------------------------------------------
            start[4] = std::chrono::system_clock::now();

            displayer.display(food_quantity);

            end[4] = std::chrono::system_clock::now();
            elapsed_seconds[4] = end[4] - start[4];
            //----------------------------------------------------------        
            win.blit();
            //----------------------------------------------------------

            if(food_quantity >= 10000 && !validation){
                end[0] = std::chrono::system_clock::now();
                elapsed_seconds[0] = end[0] - start[0];
                std::cout << ants.size() <<" ants found " << food_quantity << " pieces of food in: "<< elapsed_seconds[0].count() << " seg" << std::endl;
                validation = true; 
            }

            MPI_Recv(buffer.data(),buffer.size(),MPI_DOUBLE,1,101,display_comm,&status);
            //debug("recibi");
            //std::cout << buffer.size() << std::endl;
            food_quantity = buffer[0];
            for (size_t i = 1, j = 0; i < pher_start; i += 2, j++)
            {
                ants[j].set_position(position_t(buffer[i],buffer[i+1]));
            }
            phen.swap_map(std::vector<double> (buffer.begin() + pher_start,buffer.end()));
        });
        manager.loop();
    }
    MPI_Finalize();
    return 0;
}