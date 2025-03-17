#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <string>
#include <algorithm>
#include <Python.h>
#include <fstream>
#include <filesystem>
#include <omp.h>

#include "matplotlibcpp.h"
#include "Vecteur.hpp"
#include "Simulation.hpp"
#include "MakeSim.hpp"
#include "Admixture.hpp"


namespace plt = matplotlibcpp;
using namespace std;

int main()
{
    

    Py_Initialize();        //Initialize python for using it to plot histograms
    
    //===========================================================================
    //                        Panmixie and Selection
    //===========================================================================

    int nThreads = 8;
    int nGen = 16;                                                                      // Number of generation to compute
    int nPop = 2;                                                                       // Number of population
    int nStep = pow(2,12)+1;                                                            // Number of discretization steps
    //Vecteur<Vecteur<Vecteur<float>>> Phi(nGen, Vecteur<Vecteur<float>> ({Vecteur<float>({1.,1.}),Vecteur<float>({1.,1.})}));     // Selection matrix Phi
    //string fileName = "./output/simulation_rapport/ex_0.5_0.5_0-0.0001_0.2_0.7999_Phi_2_1_1_1/";                                                 // Name of the output file
    //string fileName = "./output/simulation_rapport/ex_0.5_0.5_0-0._0._1./";                                                 // Name of the output file

    /*
    // Contribution vector
    Vecteur<Vecteur<float>> contribution(1,{0.5,0.5,0.});                       
    for (int i=0; i<nGen; i++){contribution.push_back({0.,0.,1.});}
    */
    /*
    //Eur P1
    contribution[2][1] = 0.4f; 
    contribution[2][2] = 0.6f;

    //Afr P1 - Eur P2
    contribution[5][1] = 0.7f;
    contribution[5][0] = 0.2f;
    contribution[5][2] = 0.1f;

    //Afr P2 
    contribution[9][0] = 0.5f;
    contribution[9][2] = 0.5f;
    */

    /*
    //Contribution vector
    Vecteur<Vecteur<float>> contribution(1,{0.,0.,0.});                       
    for (int i=0; i<nGen; i++){contribution.push_back({0.,0.,0.});}

    //Get simulation coefficients in scenario/simu_i/simu_i.txt
    float gen, Ne, s1_i, s2_i;
    fstream file;
    string param_file = "./output/parameters/Afr2P-Eur2P/simu_1/simu_1.par";
    file.open(param_file, ios::in);
    string line;
    getline(file, line);
    while(file >> gen >> Ne >> s1_i >> s2_i)
    {   
        if (int(gen)<=nGen)
        {
            contribution[int(gen)][0] = s1_i;
            contribution[int(gen)][1] = s2_i;
            contribution[int(gen)][2] = 1 - s1_i - s2_i;
        }
    }
    file.close();
    */
    //Start time recording
    auto start = chrono::steady_clock::now();

    /*
    //Discretization error on moments of the 16th first marginal distribution
    Vecteur<int> nStep_list({int(pow(2,5)+1), int(pow(2,6)+1), int(pow(2,7)+1), int(pow(2,8)+1), int(pow(2,9)+1), int(pow(2,10)+1), int(pow(2,11)+1), int(pow(2,12)+1), int(pow(2,13)+1), int(pow(2,14)+1)}); //, int(pow(2,15)+1)});
    Vecteur<string> scenarios({"Afr2P-Eur2P", "Afr2P-EurIN", "Afr2P-EurDE", "AfrIN-Eur2P", "AfrIN-EurDE", "AfrDE-Eur2P", "AfrDE-EurIN", "AfrDE-EurDE", "AfrIN-EurIN"});
    int nFirsts = 100;
    int nMoments = 40;

    #pragma omp parallel for //collapse(2) schedule(dynamic) 
    for (size_t s = 0; s < scenarios.size(); s++) {
        for (int i = 1; i <= nFirsts; i++) {
            const string& scenario = scenarios[s];

            // Contribution vector
            Vecteur<Vecteur<float>> contribution(1,{0.,0.,0.});                       
            for (int i=0; i<nGen; i++){contribution.push_back({0.,0.,0.});}

            // Get simulation coefficients in scenario/simu_i/simu_i.txt
            float gen, Ne, s1_i, s2_i;
            {
                //lock_guard<mutex> lock(fileMutex);
                fstream file;
                string param_file = "./output/parameters/"+scenario+"/simu_"+to_string(i)+"/simu_"+to_string(i)+".par";
                file.open(param_file, ios::in);
                string line;
                getline(file, line);
                while(file >> gen >> Ne >> s1_i >> s2_i)
                {   
                    if (int(gen)<=nGen)
                    {
                        contribution[int(gen)][0] = s1_i;
                        contribution[int(gen)][1] = s2_i;
                        contribution[int(gen)][2] = 1 - s1_i - s2_i;
                    }
                }
                file.close();
            }   

            //Get simualtion id
            string id = scenario+"_simu_"+to_string(i);

            //Run the simulations
            Simulation sim(contribution, nGen, fileName);

            //Check simulations run well
            cout << "end sim" + to_string(i) << endl;

            //Create the directory if it don't exist
            string path = "./output/comparaison_moment/" + scenario + "/simu_" + std::to_string(i);
            filesystem::create_directories(path);

            //Save the last simulation in a table under scenario/simu_i
            {
                //lock_guard<mutex> lock(fileMutex); //Prevent confusion during parallel computing and ensure saving i right files
                fstream outputFile;
                string outputPrefix = "./output/comparaison_moment/"+scenario+"/simu_"+to_string(i)+"/endAdmixture_simu_"+to_string(i)+".txt";
                outputFile.open(outputPrefix, ios::app);
                for (const auto& H_1 : sim.final_admixture_proportion.fraction_H1)
                {
                    outputFile << H_1.first << ";" << H_1.second << "\n";
                }
                outputFile.close();
            }

            //Calculate the nMoments first moments
            Vecteur<float> moments = sim.final_admixture_proportion.firstMoments(nMoments);

            //Save the nMoments first moments
            {
                //lock_guard<mutex> lock(fileMutex);
                ofstream newfile;
                newfile.open("./output/comparaison_moment/"+scenario+"/simu_"+to_string(i)+"/firstMoments_simu_"+to_string(i)+".txt", ios::app);
                newfile << id << ";";
                for (int ii=0; ii<nMoments; ii++)
                {
                    if (ii==(nMoments-1))
                    {
                        newfile << moments[ii];
                    }
                    else
                    {
                        newfile << moments[ii] << ";" ;
                    }
                }
                newfile << "\n";
                newfile.close();
            }

            for (int ii = 0; ii<nStep_list.size(); ii++)
            {
                //cout << nStep_list[ii] << endl;
                Simulation dicrete_sim(contribution, nGen, fileName, nStep_list[ii]);

                Vecteur<float> moments = dicrete_sim.final_admixture_proportion.firstMoments(50);
                
                {
                    //lock_guard<mutex> lock(fileMutex);
                    ofstream newfile;
                    newfile.open("./output/comparaison_moment/"+scenario+"/simu_"+to_string(i)+"/firstMoments_simu_"+to_string(i)+".txt", ios::app);
                    newfile << id << "_discrete_" << ii << ";";
                    
                    for (int iii=0; iii<nMoments; iii++)
                    {
                        if (iii==(nMoments-1))
                        {
                            newfile << moments[iii];
                        }
                        else
                        {
                            newfile << moments[iii] << ";" ;
                        }
                    }
                    
                    newfile << "\n";
                    newfile.close();
                }
            }
        }
    }

    */
    //Simulation test(contribution, Phi, nGen, fileName, nStep);
    //Simulation test(contribution, Phi, nGen, fileName, true, false, false, false);
    //Simulation sim(contribution, nGen, fileName, nStep);

    /*
    //Convert to H1 to H to use other functions of Admixture
    Admixture converted_pool = sim.final_admixture_proportion.convert_fraction_H1_to_fraction();

    map<string, float> sumStat = converted_pool.summaryStat(10);

    ofstream newfile;
    newfile.open("./output/find_param/sumStat.txt", ios::app);
    newfile << "test1" << ";" << sumStat["mean"] << ";" 
            << sumStat["variance"] << ";" 
            << sumStat["skewness"] << ";" 
            << sumStat["kurtosis"] << ";" 
            << sumStat["mode"] << ";" 
            << sumStat["min"] << ";" 
            << sumStat["max"] << ";" 
            << sumStat["Q1"] << ";" 
            << sumStat["Q2"] << ";" 
            << sumStat["Q3"] << ";" 
            << sumStat["Q4"] << ";" << sumStat["Q5"] 
            << ";" << sumStat["Q6"] << ";" 
            << sumStat["Q7"] << ";" 
            << sumStat["Q8"] << ";" 
            << sumStat["Q9"] << "\n";
    newfile.close();
    */

    //cout << test.final_admixture_proportion << endl;

    //Reproduce the fortes-lima results
    Vecteur<string> scenarios({"Afr2P-Eur2P", "Afr2P-EurIN", "Afr2P-EurDE", "AfrIN-Eur2P", "AfrIN-EurDE", "AfrDE-Eur2P", "AfrDE-EurIN", "AfrDE-EurDE", "AfrIN-EurIN"});
    MakeSim test_makesim(scenarios, 10000, nStep, nThreads);

    //Our simulations with selection
    //Vecteur<string> scenarios({"AfrDE-EurDE", "AfrDE-EurDE-SelIN", "AfrDE-EurDE-SelDE", "AfrDE-EurDE-SelPC"});
    //MakeSim test_makesim(scenarios, 20000, nStep, nThreads);

    //Finish time recording & display time
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed = end-start;
    cout << nGen << " Generations, Runtime panmixie+selection: " << elapsed.count() << " seconds" << endl;

    Py_Finalize(); //Close python
    cout << "Ended :)\n";
}