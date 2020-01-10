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
bool validation = false;
size_t food_quantity = 0;

void print_food_quantity(){
    std::cout << "there is " << food_quantity << " food in the nest" << std::endl;
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
    start[0] = std::chrono::system_clock::now();

    MPI_Init(&nargs, &argv);
    
    int core_number, total_of_cores;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &core_number);
    MPI_Comm_size(MPI_COMM_WORLD, &total_of_cores);


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

    int buffer_size = 1 + 2 * nb_ants + laby.dimensions().first*laby.dimensions().second;

    //--------------------------------------------------------------------

    if(core_number == 1){
        const double eps = 0.75;  // Coefficient d'exploration
        size_t food_quantity = 0;

        // Définition du coefficient d'exploration de toutes les fourmis.
        ant::set_exploration_coef(eps);

        // On va créer toutes les fourmis dans le nid :
        std::vector<ant> ants;
        ants.reserve(nb_ants);
        for ( size_t i = 0; i < nb_ants; ++i ) ants.emplace_back(pos_nest, life);

        // On crée toutes les fourmis dans la fourmilière.
        pheromone phen(laby.dimensions(), pos_food, pos_nest, alpha, beta);
        
        while(true) {
            std::vector<double> buffer;

            buffer.emplace_back(food_quantity);

            for (int i = 0; i < nb_ants; i++)
            {
                position_t pos_ant = ants[i].get_position();
                buffer.emplace_back(pos_ant.first);
                buffer.emplace_back(pos_ant.second);
            }

            for (int i = 0; i < laby.dimensions().first; i++){
                for (int j = 0; j < laby.dimensions().second; j++)
                {
                    buffer.emplace_back((double)phen(i,j));
                }                
            }
            
            MPI_Request request;
            MPI_Status status;
            MPI_Isend(buffer.data(),buffer_size,MPI_DOUBLE,0,101,MPI_COMM_WORLD,&request);

            advance_time(laby, phen, pos_nest, pos_food, ants, food_quantity); 
            MPI_Wait(&request,&status);
        }
    }

    //------------------------------------------------------------

    if(core_number == 0){
        
        MPI_Status(status);

        std::vector<ant> ants;
        ants.reserve(nb_ants);

        //for ( size_t i = 0; i < nb_ants; ++i ) ants.emplace_back(pos_nest, life);

        pheromone phen(laby.dimensions(), pos_food, pos_nest, alpha, beta);

        std::vector<double> buffer(buffer_size);
        
        MPI_Recv(buffer.data(),buffer_size,MPI_DOUBLE,1,101,MPI_COMM_WORLD,&status);
        
        food_quantity = buffer[0];
        for (int i = 1; i < nb_ants*2 + 1; i += 2)
        {
            ants.emplace_back(position_t(buffer[i],buffer[i+1]),life);
        }
        phen.swap_map(std::vector<double> (buffer.begin() + nb_ants * 2 + 1,buffer.end()));        

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
                std::cout << "Ants found " << food_quantity << " pieces of food in: "<< elapsed_seconds[0].count() << " seg" << std::endl;
                validation = true; 
            }

            MPI_Recv(buffer.data(),buffer_size,MPI_DOUBLE,1,101,MPI_COMM_WORLD,&status);

            food_quantity = buffer[0];
            for (int i = 1, j = 0; i < nb_ants * 2 + 1; i += 2, j++)
            {
                ants[j].set_position(position_t(buffer[i],buffer[i+1]));
                //std::cout << buffer[i] <<" "<< buffer[i+1] << std::endl;
            }
            phen.swap_map(std::vector<double> (buffer.begin() + nb_ants * 2 + 1,buffer.end()));
        });
        manager.loop();
    }
    MPI_Finalize();
    return 0;
}